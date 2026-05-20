from __future__ import annotations

import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.data import covalent_radii
from ase.data.colors import jmol_colors
from ase.visualize import view
from PySide6.QtCore import QObject, QPointF, Qt, QSignalBlocker, QThread, QTimer, Signal
from PySide6.QtGui import QColor, QMouseEvent, QPainter, QPen, QBrush, QWheelEvent
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSlider,
    QSpinBox,
    QDoubleSpinBox,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from dft_workflow.assembly import AdsorbatePlacement, Assembly, assemble_adsorbates_on_slab
from dft_workflow.config import load_xtb_command
from dft_workflow.molecule_io import load_molecule
from dft_workflow.relaxation import run_gfnff_relaxation, run_mmff_relaxation, run_uff_relaxation, xtb_is_available
from dft_workflow.surfaces import build_fcc_slab
from dft_workflow.validation import validate_export
from dft_workflow.vasp_export import write_job_folder


PACKAGE_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = Path.cwd()
DEFAULT_TEMPLATE_DIR = PACKAGE_ROOT / "Templates"


ASE_BACKGROUND = "#000000"
ASE_FOREGROUND = (0.82, 0.82, 0.82, 1.0)
ASE_RADII_SCALE = 0.89


@dataclass
class AdsorbateModel:
    name: str
    atoms: Atoms
    translation: tuple[float, float, float] = (0.0, 0.0, 2.5)
    rotation_deg: tuple[float, float, float] = (0.0, 0.0, 0.0)


class RelaxationWorker(QObject):
    finished = Signal(object)

    def __init__(self, method: str, placements: list[AdsorbatePlacement], slab: Atoms, parameters: dict) -> None:
        super().__init__()
        self.method = method
        self.placements = placements
        self.slab = slab.copy()
        self.parameters = parameters

    def run(self) -> None:
        if self.method == "GFN-FF via xTB":
            result = run_gfnff_relaxation(
                self.placements,
                self.slab,
                max_atoms=self.parameters["max_atoms"],
                xtb_command=self.parameters["xtb_command"],
                timeout_seconds=self.parameters["timeout_seconds"],
                stages=self.parameters["xtb_stages"],
                steps_per_stage=self.parameters["xtb_steps_per_stage"],
            )
        elif self.method.startswith("MMFF94"):
            result = run_mmff_relaxation(
                self.placements,
                self.slab if self.parameters["openbabel_include_surface"] else None,
                max_atoms=self.parameters["max_atoms"],
                steps=self.parameters["openbabel_steps"],
                batch_size=self.parameters["openbabel_batch_steps"],
            )
        else:
            result = run_uff_relaxation(
                self.placements,
                self.slab if self.parameters["openbabel_include_surface"] else None,
                max_atoms=self.parameters["max_atoms"],
                steps=self.parameters["openbabel_steps"],
                batch_size=self.parameters["openbabel_batch_steps"],
            )
        self.finished.emit(result)


class StructureCanvas(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self.setMinimumSize(500, 380)
        self.setMouseTracking(True)
        self.atoms: Atoms | None = None
        self.rotation = _rotation_matrix(22.0, 0.0, -35.0)
        self.scale = 1.0
        self.pan = np.zeros(2)
        self.view_center = np.zeros(3)
        self.has_focused_once = False
        self.pending_initial_focus = False
        self.focus_request_id = 0
        self.display_kind = "empty"
        self.last_mouse_position: QPointF | None = None
        self.draw_atoms(None)

    def draw_atoms(self, atoms: Atoms | None, refocus: bool = False, display_kind: str = "structure") -> None:
        previous_display_kind = self.display_kind
        self.display_kind = display_kind
        self.atoms = atoms.copy() if atoms is not None else None
        if self.atoms is None or len(self.atoms) == 0:
            self.has_focused_once = False
            self.pending_initial_focus = False
            self.display_kind = "empty"
        if refocus or (self.atoms is not None and len(self.atoms) > 0 and not self.has_focused_once):
            self.pending_initial_focus = True
            self.focus_request_id += 1
            request_id = self.focus_request_id
            QTimer.singleShot(0, lambda: self._apply_pending_focus(request_id))
        elif previous_display_kind == "adsorbate_only" and self.display_kind == "assembly":
            self._rebase_view_center()
        elif self.pending_initial_focus:
            self.pending_initial_focus = False
            self.focus_request_id += 1
        self.update()

    def _rebase_view_center(self) -> None:
        if self.atoms is not None and len(self.atoms) > 0:
            self.view_center = self.atoms.get_positions().mean(axis=0)

    def resizeEvent(self, event) -> None:  # noqa: N802 - Qt API
        super().resizeEvent(event)
        if self.pending_initial_focus:
            self.focus_request_id += 1
            request_id = self.focus_request_id
            QTimer.singleShot(0, lambda: self._apply_pending_focus(request_id))

    def _apply_pending_focus(self, request_id: int) -> None:
        if request_id != self.focus_request_id:
            return
        if self.atoms is None or len(self.atoms) == 0:
            return
        self._focus()
        self.has_focused_once = True
        self.pending_initial_focus = False
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt API
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, True)
        painter.fillRect(self.rect(), QColor(ASE_BACKGROUND))

        if self.atoms is None or len(self.atoms) == 0:
            painter.setPen(QPen(QColor("#dddddd"), 1))
            painter.drawText(self.rect(), Qt.AlignCenter, "No structure loaded")
            return

        positions, radii = self._projected_atom_data()
        self._draw_cell(painter)
        self._draw_axes(painter)
        self._draw_atoms(painter, positions, radii)

    def mousePressEvent(self, event: QMouseEvent) -> None:  # noqa: N802 - Qt API
        self.last_mouse_position = event.position()

    def mouseMoveEvent(self, event: QMouseEvent) -> None:  # noqa: N802 - Qt API
        if self.last_mouse_position is None:
            self.last_mouse_position = event.position()
            return

        delta = event.position() - self.last_mouse_position
        self.last_mouse_position = event.position()
        if event.buttons() & Qt.LeftButton:
            # Apply rotations in screen coordinates so dragging feels like ASE GUI.
            self.rotation = self.rotation @ _rotation_matrix(float(-delta.y()) * 0.5, float(-delta.x()) * 0.5, 0.0)
            self.update()
        elif event.buttons() & (Qt.MiddleButton | Qt.RightButton):
            self.pan += np.array([float(delta.x()), float(delta.y())])
            self.update()

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:  # noqa: N802 - Qt API
        self.last_mouse_position = None

    def wheelEvent(self, event: QWheelEvent) -> None:  # noqa: N802 - Qt API
        factor = 1.15 if event.angleDelta().y() > 0 else 1 / 1.15
        self.scale *= factor
        self.update()

    def mouseDoubleClickEvent(self, event: QMouseEvent) -> None:  # noqa: N802 - Qt API
        self._focus()
        self.update()

    def reset_view(self) -> None:
        self.atoms = None
        self.rotation = _rotation_matrix(22.0, 0.0, -35.0)
        self.scale = 1.0
        self.pan = np.zeros(2)
        self.view_center = np.zeros(3)
        self.has_focused_once = False
        self.pending_initial_focus = False
        self.focus_request_id += 1
        self.display_kind = "empty"
        self.update()

    def _focus(self) -> None:
        self.pan = np.zeros(2)
        if self.atoms is None or len(self.atoms) == 0:
            self.scale = 1.0
            self.view_center = np.zeros(3)
            return
        positions = self.atoms.get_positions()
        radii = covalent_radii[self.atoms.numbers] * ASE_RADII_SCALE
        self.view_center = positions.mean(axis=0)
        projected = positions @ self.rotation
        lower = (projected[:, :2] - radii[:, None]).min(axis=0)
        upper = (projected[:, :2] + radii[:, None]).max(axis=0)
        span = np.maximum(upper - lower, 1e-6)
        margin = 0.78
        self.scale = min(self.width() / span[0], self.height() / span[1]) * margin

    def _projected_atom_data(self) -> tuple[np.ndarray, np.ndarray]:
        assert self.atoms is not None
        positions = self.atoms.get_positions()
        radii = covalent_radii[self.atoms.numbers] * ASE_RADII_SCALE
        projected = positions @ self.rotation
        center = (self.view_center @ self.rotation)[:2]
        screen = (projected[:, :2] - center) * self.scale
        screen[:, 0] += self.width() / 2 + self.pan[0]
        screen[:, 1] = self.height() / 2 - screen[:, 1] + self.pan[1]
        return np.column_stack([screen, projected[:, 2]]), radii * self.scale

    def _draw_atoms(self, painter: QPainter, positions: np.ndarray, radii: np.ndarray) -> None:
        assert self.atoms is not None
        order = np.argsort(positions[:, 2])
        for index in order:
            x, y, _ = positions[index]
            radius = max(float(radii[index]), 2.0)
            color = QColor.fromRgbF(*jmol_colors[self.atoms.numbers[index]])
            painter.setPen(QPen(QColor("#000000"), 1.0))
            painter.setBrush(QBrush(color))
            painter.drawEllipse(QPointF(float(x), float(y)), radius, radius)

    def _draw_cell(self, painter: QPainter) -> None:
        assert self.atoms is not None
        cell = self.atoms.cell.array
        if abs(np.linalg.det(cell)) < 1e-8:
            return
        corners = np.array(
            [
                [0, 0, 0],
                cell[0],
                cell[1],
                cell[2],
                cell[0] + cell[1],
                cell[0] + cell[2],
                cell[1] + cell[2],
                cell[0] + cell[1] + cell[2],
            ]
        )
        corners = self._project_points(corners)
        edges = (
            (0, 1),
            (0, 2),
            (0, 3),
            (1, 4),
            (1, 5),
            (2, 4),
            (2, 6),
            (3, 5),
            (3, 6),
            (4, 7),
            (5, 7),
            (6, 7),
        )
        painter.setPen(QPen(QColor.fromRgbF(*ASE_FOREGROUND), 1.0))
        for start, end in edges:
            painter.drawLine(QPointF(*corners[start]), QPointF(*corners[end]))

    def _draw_axes(self, painter: QPainter) -> None:
        if self.atoms is None or len(self.atoms) == 0:
            origin = np.zeros(3)
            length = 5.0
        else:
            origin = self.atoms.get_positions().mean(axis=0)
            length = max(float(np.ptp(self.atoms.get_positions(), axis=0).max()) * 0.18, 2.0)
        axes = (
            (np.array([origin, origin + [length, 0, 0]], dtype=float), QColor("#cc0000")),
            (np.array([origin, origin + [0, length, 0]], dtype=float), QColor("#008800")),
            (np.array([origin, origin + [0, 0, length]], dtype=float), QColor("#0040dd")),
        )
        for segment, color in axes:
            points = self._project_points(segment)
            painter.setPen(QPen(color, 2.0))
            painter.drawLine(QPointF(*points[0]), QPointF(*points[1]))

    def _project_points(self, points: np.ndarray) -> np.ndarray:
        assert self.atoms is not None
        center = (self.view_center @ self.rotation)[:2]
        projected = points @ self.rotation
        screen = (projected[:, :2] - center) * self.scale
        screen[:, 0] += self.width() / 2 + self.pan[0]
        screen[:, 1] = self.height() / 2 - screen[:, 1] + self.pan[1]
        return screen


def _rotation_matrix(x_deg: float, y_deg: float, z_deg: float) -> np.ndarray:
    x, y, z = np.deg2rad([x_deg, y_deg, z_deg])
    rx = np.array([[1, 0, 0], [0, np.cos(x), -np.sin(x)], [0, np.sin(x), np.cos(x)]])
    ry = np.array([[np.cos(y), 0, np.sin(y)], [0, 1, 0], [-np.sin(y), 0, np.cos(y)]])
    rz = np.array([[np.cos(z), -np.sin(z), 0], [np.sin(z), np.cos(z), 0], [0, 0, 1]])
    return rz @ ry @ rx


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("DFT Job Submission Workflow")
        self.resize(980, 720)

        self.adsorbates: list[AdsorbateModel] = []
        self.slab = None
        self.assembly: Assembly | None = None
        self._loading_adsorbate_controls = False
        self.relaxation_snapshot: list[tuple[Atoms, tuple[float, float, float], tuple[float, float, float]]] | None = None
        self.relaxation_trajectory: list[Atoms] = []
        self.relaxation_frame_index = 0
        self.relaxation_timer = QTimer(self)
        self.relaxation_timer.setInterval(250)
        self.relaxation_timer.timeout.connect(self._advance_relaxation_frame)
        self.relaxation_thread: QThread | None = None
        self.relaxation_worker: RelaxationWorker | None = None

        self._build_ui()
        self._connect_signals()
        self._update_summary()

    def _build_ui(self) -> None:
        root = QWidget()
        layout = QHBoxLayout(root)
        controls = QVBoxLayout()
        layout.addLayout(controls, 0)

        molecule_box = QGroupBox("Adsorbates")
        molecule_layout = QVBoxLayout(molecule_box)
        self.load_molecule_button = QPushButton("Add XYZ/SDF")
        self.remove_adsorbate_button = QPushButton("Remove Selected")
        self.remove_adsorbate_button.setEnabled(False)
        self.adsorbate_selector = QComboBox()
        self.adsorbate_selector.setEnabled(False)
        self.open_ase_viewer_button = QPushButton("Open Current Structure in ASE Viewer")
        self.open_ase_viewer_button.setEnabled(False)
        molecule_layout.addWidget(self.load_molecule_button)
        molecule_layout.addWidget(self.remove_adsorbate_button)
        molecule_layout.addWidget(self.adsorbate_selector)
        molecule_layout.addWidget(self.open_ase_viewer_button)
        controls.addWidget(molecule_box)

        surface_box = QGroupBox("Surface")
        surface_layout = QFormLayout(surface_box)
        self.metal_input = QComboBox()
        self.metal_input.setEditable(True)
        self.metal_input.addItems(["Cu", "Ag", "Au"])
        self.facet_input = QComboBox()
        self.facet_input.addItems(["111", "100", "110"])
        self.size_x_input = self._spinbox(1, 20, 4)
        self.size_y_input = self._spinbox(1, 20, 4)
        self.layers_input = self._spinbox(1, 20, 3)
        self.vacuum_input = self._double_spinbox(1.0, 80.0, 20.0, 1.0)
        self.bottom_vacuum_input = self._double_spinbox(0.0, 40.0, 0.0, 0.5)
        self.lattice_input = self._double_spinbox(0.0, 10.0, 0.0, 0.01)
        self.build_slab_button = QPushButton("Build Slab")
        surface_layout.addRow("Metal", self.metal_input)
        surface_layout.addRow("Facet", self.facet_input)
        surface_layout.addRow("Size X", self.size_x_input)
        surface_layout.addRow("Size Y", self.size_y_input)
        surface_layout.addRow("Layers", self.layers_input)
        surface_layout.addRow("Top vacuum (A)", self.vacuum_input)
        surface_layout.addRow("Bottom vacuum (A)", self.bottom_vacuum_input)
        surface_layout.addRow("FCC lattice a (A, 0=ASE)", self.lattice_input)
        surface_layout.addRow(self.build_slab_button)
        controls.addWidget(surface_box)

        placement_box = QGroupBox("Molecule Placement")
        placement_layout = QFormLayout(placement_box)
        self.x_input, x_widget = self._double_slider(-30.0, 30.0, 0.0, 0.1, 10)
        self.y_input, y_widget = self._double_slider(-30.0, 30.0, 0.0, 0.1, 10)
        self.z_input, z_widget = self._double_slider(0.1, 30.0, 2.5, 0.1, 10)
        self.rx_input, rx_widget = self._double_slider(-180.0, 180.0, 0.0, 1.0, 1)
        self.ry_input, ry_widget = self._double_slider(-180.0, 180.0, 0.0, 1.0, 1)
        self.rz_input, rz_widget = self._double_slider(-180.0, 180.0, 0.0, 1.0, 1)
        self.freeze_layers_input = self._spinbox(0, 20, 1)
        placement_layout.addRow("X shift (A)", x_widget)
        placement_layout.addRow("Y shift (A)", y_widget)
        placement_layout.addRow("Height (A)", z_widget)
        placement_layout.addRow("Rotate X (deg)", rx_widget)
        placement_layout.addRow("Rotate Y (deg)", ry_widget)
        placement_layout.addRow("Rotate Z (deg)", rz_widget)
        placement_layout.addRow("Frozen bottom layers", self.freeze_layers_input)
        controls.addWidget(placement_box)

        export_box = QGroupBox("Export")
        export_layout = QFormLayout(export_box)
        self.template_dir_input = QLineEdit(str(DEFAULT_TEMPLATE_DIR))
        self.output_dir_input = QLineEdit(str(PROJECT_ROOT / "exported_jobs" / "job_001"))
        self.browse_template_button = QPushButton("Browse")
        self.browse_output_button = QPushButton("Browse")
        self.potcar_command_input = QLineEdit("python potcar_generator.py {elements} > POTCAR")
        self.coordinate_mode_input = QComboBox()
        self.coordinate_mode_input.addItems(["Direct", "Cartesian"])
        self.overwrite_input = QCheckBox("Overwrite non-empty folder")
        self.export_button = QPushButton("Export VASP Job Folder")
        export_layout.addRow("Templates", self._path_picker(self.template_dir_input, self.browse_template_button))
        export_layout.addRow("Output folder", self._path_picker(self.output_dir_input, self.browse_output_button))
        export_layout.addRow("POTCAR command", self.potcar_command_input)
        export_layout.addRow("POSCAR coordinates", self.coordinate_mode_input)
        export_layout.addRow(self.overwrite_input)
        export_layout.addRow(self.export_button)

        relaxation_box = QGroupBox("Relaxation")
        relaxation_layout = QFormLayout(relaxation_box)
        self.relaxation_method_input = QComboBox()
        self.relaxation_method_input.addItems(["GFN-FF via xTB", "MMFF94 (Open Babel)", "UFF (Open Babel)"])
        self.run_relaxation_button = QPushButton("Run Relaxation")
        self.run_relaxation_button.setEnabled(False)
        self.reset_relaxation_button = QPushButton("Reset Relaxation")
        self.reset_relaxation_button.setEnabled(False)
        self.play_relaxation_button = QPushButton("Play Relaxation")
        self.play_relaxation_button.setEnabled(False)
        self.relaxation_max_atoms_input = self._spinbox(1, 5000, 600)
        self.xtb_stages_input = self._spinbox(1, 200, 30)
        self.xtb_steps_input = self._spinbox(1, 100, 5)
        self.xtb_timeout_input = self._spinbox(1, 86400, 3600)
        self.openbabel_steps_input = self._spinbox(1, 10000, 500)
        self.openbabel_batch_input = self._spinbox(1, 1000, 10)
        self.openbabel_surface_input = QCheckBox("Include top surface layer")
        self.openbabel_surface_input.setChecked(False)
        relaxation_layout.addRow("Method", self.relaxation_method_input)
        relaxation_layout.addRow("Max atoms", self.relaxation_max_atoms_input)
        relaxation_layout.addRow("xTB stages", self.xtb_stages_input)
        relaxation_layout.addRow("xTB steps/stage", self.xtb_steps_input)
        relaxation_layout.addRow("xTB timeout (s)", self.xtb_timeout_input)
        relaxation_layout.addRow("OB max steps", self.openbabel_steps_input)
        relaxation_layout.addRow("OB batch steps", self.openbabel_batch_input)
        relaxation_layout.addRow(self.openbabel_surface_input)
        relaxation_layout.addRow(self.run_relaxation_button)
        relaxation_layout.addRow(self.reset_relaxation_button)
        relaxation_layout.addRow(self.play_relaxation_button)

        reset_box = QGroupBox("Reset")
        reset_layout = QVBoxLayout(reset_box)
        self.reset_app_button = QPushButton("Clear All Structures")
        self.reset_app_button.setStyleSheet(
            "QPushButton { background-color: #b00020; color: white; font-weight: bold; padding: 6px; }"
            "QPushButton:hover { background-color: #d00028; }"
        )
        reset_layout.addWidget(self.reset_app_button)
        controls.addWidget(reset_box)
        controls.addStretch(1)

        center = QVBoxLayout()
        layout.addLayout(center, 1)
        self.viewer = StructureCanvas()
        self.summary = QTextEdit()
        self.summary.setReadOnly(True)
        self.relaxation_log = QTextEdit()
        self.relaxation_log.setReadOnly(True)
        self.relaxation_log.setMaximumHeight(170)
        center.addWidget(QLabel("Structure Preview"))
        center.addWidget(self.viewer, 1)

        side = QVBoxLayout()
        layout.addLayout(side, 0)
        side.addWidget(export_box)
        side.addWidget(relaxation_box)
        side.addWidget(QLabel("Structure Summary / Validation"))
        side.addWidget(self.summary, 1)
        side.addWidget(QLabel("Relaxation Log"))
        side.addWidget(self.relaxation_log, 1)
        self.setCentralWidget(root)

    def _connect_signals(self) -> None:
        self.load_molecule_button.clicked.connect(self._load_molecule)
        self.remove_adsorbate_button.clicked.connect(self._remove_selected_adsorbate)
        self.adsorbate_selector.currentIndexChanged.connect(self._select_adsorbate)
        self.open_ase_viewer_button.clicked.connect(self._open_ase_viewer)
        self.build_slab_button.clicked.connect(self._build_slab)
        self.browse_template_button.clicked.connect(self._browse_template_dir)
        self.browse_output_button.clicked.connect(self._browse_output_dir)
        self.export_button.clicked.connect(self._export)
        self.run_relaxation_button.clicked.connect(self._run_relaxation)
        self.reset_relaxation_button.clicked.connect(self._reset_relaxation)
        self.play_relaxation_button.clicked.connect(self._toggle_relaxation_playback)
        self.reset_app_button.clicked.connect(self._reset_app)
        for widget in (
            self.x_input,
            self.y_input,
            self.z_input,
            self.rx_input,
            self.ry_input,
            self.rz_input,
            self.freeze_layers_input,
        ):
            widget.valueChanged.connect(self._rebuild_assembly)

    def _load_molecule(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load molecule",
            str(PROJECT_ROOT),
            "Molecule files (*.xyz *.sdf)",
        )
        if not path:
            return
        try:
            molecule_path = Path(path)
            name = self._unique_adsorbate_name(molecule_path.stem)
            self.adsorbates.append(AdsorbateModel(name=name, atoms=load_molecule(path)))
            self._clear_relaxation_snapshot("Relaxation reset because an adsorbate was added.")
            self._refresh_adsorbate_selector(len(self.adsorbates) - 1)
            self._rebuild_assembly(refocus=False)
        except Exception as exc:
            QMessageBox.critical(self, "Molecule load failed", str(exc))

    def _remove_selected_adsorbate(self) -> None:
        index = self.adsorbate_selector.currentIndex()
        if index < 0:
            return
        del self.adsorbates[index]
        self._clear_relaxation_snapshot("Relaxation reset because an adsorbate was removed.")
        selected = min(index, len(self.adsorbates) - 1)
        self._refresh_adsorbate_selector(selected)
        self._rebuild_assembly(refocus=False)

    def _select_adsorbate(self, index: int) -> None:
        if index < 0 or index >= len(self.adsorbates):
            return
        self._load_adsorbate_transform(index)
        self._update_summary(refocus=False)

    def _refresh_adsorbate_selector(self, selected_index: int = 0) -> None:
        with QSignalBlocker(self.adsorbate_selector):
            self.adsorbate_selector.clear()
            for adsorbate in self.adsorbates:
                self.adsorbate_selector.addItem(f"{adsorbate.name} ({len(adsorbate.atoms)} atoms)")
            if self.adsorbates:
                self.adsorbate_selector.setCurrentIndex(max(0, selected_index))
        has_adsorbates = bool(self.adsorbates)
        self.adsorbate_selector.setEnabled(has_adsorbates)
        self.remove_adsorbate_button.setEnabled(has_adsorbates)
        if has_adsorbates:
            self._load_adsorbate_transform(self.adsorbate_selector.currentIndex())

    def _load_adsorbate_transform(self, index: int) -> None:
        adsorbate = self.adsorbates[index]
        self._loading_adsorbate_controls = True
        try:
            controls = (
                (self.x_input, adsorbate.translation[0]),
                (self.y_input, adsorbate.translation[1]),
                (self.z_input, adsorbate.translation[2]),
                (self.rx_input, adsorbate.rotation_deg[0]),
                (self.ry_input, adsorbate.rotation_deg[1]),
                (self.rz_input, adsorbate.rotation_deg[2]),
            )
            for control, value in controls:
                control.setValue(value)
        finally:
            self._loading_adsorbate_controls = False

    def _store_selected_adsorbate_transform(self) -> None:
        if self._loading_adsorbate_controls:
            return
        index = self.adsorbate_selector.currentIndex()
        if index < 0 or index >= len(self.adsorbates):
            return
        self.adsorbates[index].translation = (self.x_input.value(), self.y_input.value(), self.z_input.value())
        self.adsorbates[index].rotation_deg = (self.rx_input.value(), self.ry_input.value(), self.rz_input.value())

    def _unique_adsorbate_name(self, base_name: str) -> str:
        existing = {adsorbate.name for adsorbate in self.adsorbates}
        if base_name not in existing:
            return base_name
        counter = 2
        while f"{base_name}_{counter}" in existing:
            counter += 1
        return f"{base_name}_{counter}"

    def _open_ase_viewer(self) -> None:
        atoms = self._current_atoms()
        if atoms is None:
            QMessageBox.warning(self, "No structure", "Load a molecule or build a slab first.")
            return
        try:
            view(atoms)
        except Exception as exc:
            QMessageBox.critical(self, "ASE Viewer failed", str(exc))

    def _build_slab(self) -> None:
        try:
            lattice = self.lattice_input.value() or None
            self.slab = build_fcc_slab(
                self.metal_input.currentText().strip(),
                self.facet_input.currentText(),
                self.size_x_input.value(),
                self.size_y_input.value(),
                self.layers_input.value(),
                self.vacuum_input.value(),
                lattice,
                self.bottom_vacuum_input.value(),
            )
            self._clear_relaxation_snapshot("Relaxation reset because the slab was rebuilt.")
            self._rebuild_assembly(refocus=False)
        except Exception as exc:
            QMessageBox.critical(self, "Slab build failed", str(exc))

    def _browse_template_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select template folder", self.template_dir_input.text())
        if path:
            self.template_dir_input.setText(path)

    def _browse_output_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select output folder", str(PROJECT_ROOT))
        if path:
            self.output_dir_input.setText(path)

    def _reset_app(self) -> None:
        self.adsorbates.clear()
        self.slab = None
        self.assembly = None
        self.relaxation_snapshot = None
        self._loading_adsorbate_controls = True
        try:
            with QSignalBlocker(self.adsorbate_selector):
                self.adsorbate_selector.clear()
            self.adsorbate_selector.setEnabled(False)
            self.remove_adsorbate_button.setEnabled(False)
            for control, value in (
                (self.x_input, 0.0),
                (self.y_input, 0.0),
                (self.z_input, 2.5),
                (self.rx_input, 0.0),
                (self.ry_input, 0.0),
                (self.rz_input, 0.0),
            ):
                control.setValue(value)
            self.freeze_layers_input.setValue(1)
        finally:
            self._loading_adsorbate_controls = False
        self.viewer.reset_view()
        self.relaxation_log.clear()
        self._clear_relaxation_trajectory()
        self._update_summary(refocus=False)

    def _run_relaxation(self) -> None:
        if self.assembly is None or self.slab is None or not self.adsorbates:
            self._append_relaxation_log("Relaxation requires at least one adsorbate and a built slab.")
            return
        if self.relaxation_thread is not None:
            self._append_relaxation_log("Relaxation is already running.")
            return
        self._store_selected_adsorbate_transform()
        placements = self._current_adsorbate_placements()
        method = self.relaxation_method_input.currentText()
        self._clear_relaxation_trajectory()
        self._append_relaxation_log(f"Starting relaxation: {method}")
        parameters = self._relaxation_parameters()
        if method == "GFN-FF via xTB" and not xtb_is_available(parameters["xtb_command"]):
            self._append_relaxation_log(
                "xTB executable not found. Set xtb in dft_workflow_config.ini, add it to PATH, "
                "or use Open Babel UFF/MMFF94."
            )
            return

        self._set_relaxation_running(True)
        self.relaxation_thread = QThread(self)
        self.relaxation_worker = RelaxationWorker(method, placements, self.slab, parameters)
        self.relaxation_worker.moveToThread(self.relaxation_thread)
        self.relaxation_thread.started.connect(self.relaxation_worker.run)
        self.relaxation_worker.finished.connect(self._finish_relaxation)
        self.relaxation_worker.finished.connect(self.relaxation_thread.quit)
        self.relaxation_worker.finished.connect(self.relaxation_worker.deleteLater)
        self.relaxation_thread.finished.connect(self._cleanup_relaxation_thread)
        self.relaxation_thread.start()

    def _finish_relaxation(self, result) -> None:
        for line in result.log:
            self._append_relaxation_log(line)
        if not result.success:
            self._append_relaxation_log(f"Relaxation failed: {result.message}")
            self._set_relaxation_running(False)
            return

        if self.relaxation_snapshot is None:
            self.relaxation_snapshot = [
                (adsorbate.atoms.copy(), adsorbate.translation, adsorbate.rotation_deg)
                for adsorbate in self.adsorbates
            ]
        for model, relaxed in zip(self.adsorbates, result.adsorbates):
            model.atoms = relaxed.atoms.copy()
            model.translation = relaxed.translation
            model.rotation_deg = relaxed.rotation_deg
        self.relaxation_trajectory = [frame.copy() for frame in result.trajectory]
        self.relaxation_frame_index = 0
        self.play_relaxation_button.setEnabled(bool(self.relaxation_trajectory))

        self._refresh_adsorbate_selector(self.adsorbate_selector.currentIndex())
        self._rebuild_assembly(refocus=False)
        self._append_relaxation_log(result.message)
        self._set_relaxation_running(False)

    def _cleanup_relaxation_thread(self) -> None:
        if self.relaxation_thread is not None:
            self.relaxation_thread.deleteLater()
        self.relaxation_thread = None
        self.relaxation_worker = None

    def _relaxation_parameters(self) -> dict:
        return {
            "max_atoms": self.relaxation_max_atoms_input.value(),
            "xtb_command": load_xtb_command(),
            "xtb_stages": self.xtb_stages_input.value(),
            "xtb_steps_per_stage": self.xtb_steps_input.value(),
            "timeout_seconds": self.xtb_timeout_input.value(),
            "openbabel_steps": self.openbabel_steps_input.value(),
            "openbabel_batch_steps": self.openbabel_batch_input.value(),
            "openbabel_include_surface": self.openbabel_surface_input.isChecked(),
        }

    def _set_relaxation_running(self, running: bool) -> None:
        self.run_relaxation_button.setEnabled(False if running else self.assembly is not None and self.slab is not None and bool(self.adsorbates))
        self.reset_relaxation_button.setEnabled(False if running else self.relaxation_snapshot is not None)
        self.export_button.setEnabled(not running)
        self.reset_app_button.setEnabled(not running)
        self.play_relaxation_button.setEnabled(False if running else bool(self.relaxation_trajectory))

    def _reset_relaxation(self) -> None:
        if self.relaxation_snapshot is None:
            self._append_relaxation_log("No relaxation snapshot to reset.")
            return
        for adsorbate, snapshot in zip(self.adsorbates, self.relaxation_snapshot):
            atoms, translation, rotation_deg = snapshot
            adsorbate.atoms = atoms.copy()
            adsorbate.translation = translation
            adsorbate.rotation_deg = rotation_deg
        self.relaxation_snapshot = None
        self._clear_relaxation_trajectory()
        self._refresh_adsorbate_selector(self.adsorbate_selector.currentIndex())
        self._rebuild_assembly(refocus=False)
        self._append_relaxation_log("Relaxation reset to pre-relaxation structure.")

    def _current_adsorbate_placements(self) -> list[AdsorbatePlacement]:
        return [
            AdsorbatePlacement(adsorbate.atoms, adsorbate.translation, adsorbate.rotation_deg)
            for adsorbate in self.adsorbates
        ]

    def _append_relaxation_log(self, message: str) -> None:
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.relaxation_log.append(f"[{timestamp}] {message}")

    def _toggle_relaxation_playback(self) -> None:
        if not self.relaxation_trajectory:
            self._append_relaxation_log("No relaxation trajectory available for playback.")
            return
        if self.relaxation_timer.isActive():
            self.relaxation_timer.stop()
            self.play_relaxation_button.setText("Play Relaxation")
            return
        self.relaxation_frame_index = 0
        self.play_relaxation_button.setText("Stop Playback")
        self.relaxation_timer.start()

    def _advance_relaxation_frame(self) -> None:
        if not self.relaxation_trajectory:
            self._clear_relaxation_trajectory()
            return
        frame = self.relaxation_trajectory[self.relaxation_frame_index]
        self.viewer.draw_atoms(frame, refocus=False, display_kind="trajectory")
        self.relaxation_frame_index += 1
        if self.relaxation_frame_index >= len(self.relaxation_trajectory):
            self.relaxation_timer.stop()
            self.play_relaxation_button.setText("Play Relaxation")
            self.relaxation_frame_index = 0
            self._update_summary(refocus=False)

    def _clear_relaxation_trajectory(self) -> None:
        self.relaxation_timer.stop()
        self.relaxation_trajectory = []
        self.relaxation_frame_index = 0
        self.play_relaxation_button.setText("Play Relaxation")
        self.play_relaxation_button.setEnabled(False)

    def _clear_relaxation_snapshot(self, message: str | None = None) -> None:
        if self.relaxation_snapshot is not None:
            self.relaxation_snapshot = None
            self._clear_relaxation_trajectory()
            if message is not None:
                self._append_relaxation_log(message)

    def _rebuild_assembly(self, *args, refocus: bool = False) -> None:
        self._store_selected_adsorbate_transform()
        if not self.adsorbates or self.slab is None:
            self.assembly = None
            self._update_summary(refocus=refocus)
            return
        try:
            placements = [
                AdsorbatePlacement(adsorbate.atoms, adsorbate.translation, adsorbate.rotation_deg)
                for adsorbate in self.adsorbates
            ]
            self.assembly = assemble_adsorbates_on_slab(
                placements,
                self.slab,
                self.metal_input.currentText().strip(),
                self.freeze_layers_input.value(),
            )
            self._update_summary(refocus=refocus)
        except Exception as exc:
            self.summary.setPlainText(str(exc))

    def _export(self) -> None:
        if self.assembly is None:
            QMessageBox.warning(self, "Nothing to export", "Load a molecule and build a slab first.")
            return
        warnings = validate_export(
            self.adsorbates[0].atoms if self.adsorbates else None,
            self.slab,
            self.assembly.atoms,
            self.assembly.mobile_flags,
            self.vacuum_input.value(),
        )
        if warnings:
            message = "Warnings:\n" + "\n".join(f"- {warning}" for warning in warnings)
            if QMessageBox.question(self, "Export warnings", message + "\n\nContinue?") != QMessageBox.Yes:
                return
        try:
            write_job_folder(
                self.output_dir_input.text(),
                self.assembly.atoms,
                self.assembly.mobile_flags,
                self.assembly.element_order,
                self.template_dir_input.text(),
                self.potcar_command_input.text(),
                self.overwrite_input.isChecked(),
                self.coordinate_mode_input.currentText(),
            )
            QMessageBox.information(self, "Export complete", f"Wrote files to {self.output_dir_input.text()}")
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))

    def _update_summary(self, refocus: bool = False) -> None:
        lines = []
        self.open_ase_viewer_button.setEnabled(self._current_atoms() is not None)
        can_relax = self.assembly is not None and self.slab is not None and bool(self.adsorbates)
        relaxation_running = self.relaxation_thread is not None
        self.run_relaxation_button.setEnabled(can_relax and not relaxation_running)
        self.reset_relaxation_button.setEnabled(self.relaxation_snapshot is not None and not relaxation_running)
        self.play_relaxation_button.setEnabled(bool(self.relaxation_trajectory) and not relaxation_running)
        if not self.adsorbates:
            lines.append("Adsorbates: none loaded")
        else:
            total_adsorbate_atoms = sum(len(adsorbate.atoms) for adsorbate in self.adsorbates)
            lines.append(f"Adsorbates: {len(self.adsorbates)} structures, {total_adsorbate_atoms} atoms")
            for index, adsorbate in enumerate(self.adsorbates, start=1):
                marker = "*" if index - 1 == self.adsorbate_selector.currentIndex() else " "
                lines.append(f"{marker} {index}. {adsorbate.name}: {len(adsorbate.atoms)} atoms")
        if self.slab is None:
            lines.append("Slab: not built")
        else:
            lines.append(f"Slab: {len(self.slab)} atoms")
            if self.assembly is None:
                cell_lengths = self.slab.cell.lengths()
                lines.append(
                    "Cell (A): "
                    f"a={cell_lengths[0]:.3f}, b={cell_lengths[1]:.3f}, c={cell_lengths[2]:.3f}"
                )
        if self.assembly is not None:
            self.viewer.draw_atoms(self.assembly.atoms, refocus=refocus, display_kind="assembly")
            lines.append(f"Combined: {len(self.assembly.atoms)} atoms")
            lines.append(f"POSCAR element order: {' '.join(self.assembly.element_order)}")
            cell_lengths = self.assembly.atoms.cell.lengths()
            lines.append(
                "Final cell (A): "
                f"a={cell_lengths[0]:.3f}, b={cell_lengths[1]:.3f}, c={cell_lengths[2]:.3f}"
            )
            if self.assembly.adsorbate_surface_distances:
                lines.append("Closest adsorbate-surface atom distances:")
                for adsorbate, distance in zip(self.adsorbates, self.assembly.adsorbate_surface_distances):
                    lines.append(f"- {adsorbate.name}: {distance:.3f} A")
            frozen = sum(1 for flag in self.assembly.mobile_flags if not any(flag))
            lines.append(f"Frozen atoms: {frozen}")
            warnings = validate_export(
                self.adsorbates[0].atoms if self.adsorbates else None,
                self.slab,
                self.assembly.atoms,
                self.assembly.mobile_flags,
                self.vacuum_input.value(),
            )
            if warnings:
                lines.append("")
                lines.append("Warnings:")
                lines.extend(f"- {warning}" for warning in warnings)
        elif self.slab is not None:
            self.viewer.draw_atoms(self.slab, refocus=refocus, display_kind="slab_only")
        elif self.adsorbates:
            self.viewer.draw_atoms(
                self.adsorbates[self.adsorbate_selector.currentIndex()].atoms,
                refocus=refocus,
                display_kind="adsorbate_only",
            )
        else:
            self.viewer.draw_atoms(None, refocus=True, display_kind="empty")
        self.summary.setPlainText("\n".join(lines))

    def _current_atoms(self) -> Atoms | None:
        if self.assembly is not None:
            return self.assembly.atoms
        if self.slab is not None:
            return self.slab
        if self.adsorbates:
            return self.adsorbates[self.adsorbate_selector.currentIndex()].atoms
        return None

    @staticmethod
    def _spinbox(minimum: int, maximum: int, value: int) -> QSpinBox:
        box = QSpinBox()
        box.setRange(minimum, maximum)
        box.setValue(value)
        return box

    @staticmethod
    def _double_spinbox(minimum: float, maximum: float, value: float, step: float) -> QDoubleSpinBox:
        box = QDoubleSpinBox()
        box.setRange(minimum, maximum)
        box.setDecimals(3)
        box.setSingleStep(step)
        box.setValue(value)
        box.setAlignment(Qt.AlignRight)
        box.setMinimumWidth(90)
        return box

    @staticmethod
    def _double_slider(
        minimum: float,
        maximum: float,
        value: float,
        step: float,
        scale: int,
    ) -> tuple[QDoubleSpinBox, QWidget]:
        spinbox = MainWindow._double_spinbox(minimum, maximum, value, step)
        slider = QSlider(Qt.Horizontal)
        slider.setRange(round(minimum * scale), round(maximum * scale))
        slider.setValue(round(value * scale))

        def set_spinbox(raw_value: int) -> None:
            with QSignalBlocker(spinbox):
                spinbox.setValue(raw_value / scale)
            spinbox.valueChanged.emit(spinbox.value())

        def set_slider(spinbox_value: float) -> None:
            with QSignalBlocker(slider):
                slider.setValue(round(spinbox_value * scale))

        slider.valueChanged.connect(set_spinbox)
        spinbox.valueChanged.connect(set_slider)

        container = QWidget()
        layout = QGridLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(slider, 0, 0)
        layout.addWidget(spinbox, 0, 1)
        return spinbox, container

    @staticmethod
    def _path_picker(line_edit: QLineEdit, button: QPushButton) -> QWidget:
        container = QWidget()
        layout = QGridLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(line_edit, 0, 0)
        layout.addWidget(button, 0, 1)
        return container


def main() -> int:
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
