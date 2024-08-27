import sys
import numpy as np
from scipy.integrate import quad_vec
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget,\
    QSlider, QLabel, QLineEdit, QDesktopWidget, QPushButton, QHBoxLayout, QColorDialog
import pyqtgraph as pg


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Фонон")

        screen_geometry = QDesktopWidget().screenGeometry()
        self.setGeometry(100, 100, int(screen_geometry.width() * 0.9), int(screen_geometry.height() * 0.9))

        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)

        self.layout = QHBoxLayout(self.central_widget)
        self.interface_layout = QVBoxLayout()

        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setBackground('black')
        self.layout.addWidget(self.plot_widget)
        self.layout.setSpacing(0)
        self.interface_layout.setContentsMargins(1, 1, 1, 1)
        self.interface_layout.setSpacing(0)

        # Интерфейс
        self.k_b_label = QLabel("Введите К_b:")
        self.k_b_label.setFixedSize(int(screen_geometry.width() * 0.15), 15)
        self.k_b_input = QLineEdit()
        self.k_b_input.setFixedSize(int(screen_geometry.width() * 0.15), 25)
        self.k_b_input.setText("-1")
        self.k_b_input.textChanged.connect(self.update_plot)
        self.interface_layout.addWidget(self.k_b_label)
        self.interface_layout.addWidget(self.k_b_input)

        self.k_e_label = QLabel("Введите К_e:")
        self.k_e_label.setFixedSize(int(screen_geometry.width() * 0.15), 15)
        self.k_e_input = QLineEdit()
        self.k_e_input.setFixedSize(int(screen_geometry.width() * 0.15), 25)
        self.k_e_input.setText("1")
        self.k_e_input.textChanged.connect(self.update_plot)
        self.interface_layout.addWidget(self.k_e_label)
        self.interface_layout.addWidget(self.k_e_input)

        self.k0_slider = QSlider()
        self.k0_slider.setOrientation(1)
        self.k0_slider.setMinimum(-100)
        self.k0_slider.setMaximum(100)
        self.k0_slider.setValue(0)
        self.k0_slider.setFixedSize(int(screen_geometry.width() * 0.15), 20)
        self.k0_slider.valueChanged.connect(self.update_plot)
        self.k0_label = QLabel("K0")
        self.k0_label.setFixedSize(int(screen_geometry.width() * 0.15), 15)
        self.interface_layout.addWidget(self.k0_label)
        self.interface_layout.addWidget(self.k0_slider)
        self.k0_value_label = QLabel(f"Текущее значение K0: {self.k0_slider.value() / 100.0}")
        self.k0_value_label.setFixedSize(int(screen_geometry.width() * 0.15), 15)
        self.interface_layout.addWidget(self.k0_value_label)

        self.delta_k_label = QLabel("Введите Δk:")
        self.delta_k_label.setFixedSize(int(screen_geometry.width() * 0.15), 15)
        self.delta_k_input = QLineEdit()
        self.delta_k_input.setFixedSize(int(screen_geometry.width() * 0.15), 25)
        self.delta_k_input.textChanged.connect(self.update_plot)
        self.interface_layout.addWidget(self.delta_k_label)
        self.interface_layout.addWidget(self.delta_k_input)
        self.delta_k_value_label = QLabel(f"Текущее значение Δk: "
                                          f"{float(self.delta_k_input.text()) if self.delta_k_input.text() else 1.0}")
        self.delta_k_value_label.setFixedSize(int(screen_geometry.width() * 0.15), 15)
        self.interface_layout.addWidget(self.delta_k_value_label)

        self.b_label = QLabel("Введите b:")
        self.b_label.setFixedSize(int(screen_geometry.width() * 0.15), 15)
        self.b_input = QLineEdit()
        self.b_input.setFixedSize(int(screen_geometry.width() * 0.15), 25)
        self.b_input.setText("1")
        self.b_input.textChanged.connect(self.update_plot)
        self.interface_layout.addWidget(self.b_label)
        self.interface_layout.addWidget(self.b_input)
        self.b_value_label = QLabel(f"Текущее значение b: {float(self.b_input.text())}")
        self.b_value_label.setFixedSize(int(screen_geometry.width() * 0.15), 15)
        self.interface_layout.addWidget(self.b_value_label)

        self.z0_label = QLabel("Введите z0:")
        self.z0_label.setFixedSize(int(screen_geometry.width() * 0.15), 15)
        self.z0_input = QLineEdit()
        self.z0_input.setFixedSize(int(screen_geometry.width() * 0.15), 25)
        self.z0_input.textChanged.connect(self.update_plot)
        self.interface_layout.addWidget(self.z0_label)
        self.interface_layout.addWidget(self.z0_input)
        self.z0_value_label = QLabel(f"Текущее значениие z0: "
                                     f"{float(self.z0_input.text()) if self.z0_input.text() else 0.0}")
        self.z0_value_label.setFixedSize(int(screen_geometry.width() * 0.15), 15)
        self.interface_layout.addWidget(self.z0_value_label)

        self.z1_label = QLabel("Введите z1:")
        self.z1_label.setFixedSize(int(screen_geometry.width() * 0.15), 15)
        self.z1_input = QLineEdit()
        self.z1_input.setFixedSize(int(screen_geometry.width() * 0.15), 25)
        self.z1_input.textChanged.connect(self.update_plot)
        self.interface_layout.addWidget(self.z1_label)
        self.interface_layout.addWidget(self.z1_input)
        self.z1_value_label = QLabel(f"Текущее значение z1: "
                                     f"{float(self.z1_input.text()) if self.z1_input.text() else 0.0}")
        self.z1_value_label.setFixedSize(int(screen_geometry.width() * 0.15), 15)
        self.interface_layout.addWidget(self.z1_value_label)

        self.time_slider = QSlider()
        self.time_slider.setFixedSize(int(screen_geometry.width() * 0.15), 15)
        self.time_slider.setOrientation(1)
        self.time_slider.setMinimum(0)
        self.time_slider.setMaximum(100)
        self.time_slider.setValue(0)
        self.time_slider.valueChanged.connect(self.update_plot)
        self.time_label = QLabel("Время")
        self.time_label.setFixedSize(int(screen_geometry.width() * 0.15), 15)
        self.interface_layout.addWidget(self.time_label)
        self.interface_layout.addWidget(self.time_slider)
        self.time_value_label = QLabel(f"Текущее время: {self.time_slider.value() / 100.0}")
        self.time_value_label.setFixedSize(int(screen_geometry.width() * 0.15), 15)
        self.interface_layout.addWidget(self.time_value_label)

        # Кнопка, чтобы менять тип фонона
        self.equation_state = 0  # 0 - Акустический, 1 - Оптический
        self.equation_button = QPushButton("Поменять тип фонона")
        self.equation_button.setFixedSize(int(screen_geometry.width() * 0.15), 40)
        self.equation_button.clicked.connect(self.switch_equation)
        self.interface_layout.addWidget(self.equation_button)

        # Create a toggle button for the wave packet
        self.wave_packet_button = QPushButton("Волновой пакет")
        self.wave_packet_button.setFixedSize(int(screen_geometry.width() * 0.15), 40)
        self.wave_packet_button.setCheckable(True)
        self.wave_packet_button.clicked.connect(self.update_plot)
        self.interface_layout.addWidget(self.wave_packet_button)

        self.show_real_button = QPushButton("Показать реальную часть")
        self.show_real_button.setFixedSize(int(screen_geometry.width() * 0.15), 40)
        self.show_real_button.setCheckable(True)
        self.show_real_button.clicked.connect(self.update_plot)
        self.interface_layout.addWidget(self.show_real_button)

        self.show_imaginary_button = QPushButton("Показать мнимую часть")
        self.show_imaginary_button.setFixedSize(int(screen_geometry.width() * 0.15), 40)
        self.show_imaginary_button.setCheckable(True)
        self.show_imaginary_button.clicked.connect(self.update_plot)
        self.interface_layout.addWidget(self.show_imaginary_button)

        self.show_amplitude_button = QPushButton("Показать амплитуду")
        self.show_amplitude_button.setFixedSize(int(screen_geometry.width() * 0.15), 40)
        self.show_amplitude_button.setCheckable(True)
        self.show_amplitude_button.clicked.connect(self.update_plot)
        self.interface_layout.addWidget(self.show_amplitude_button)

        self.atom_displacements_button = QPushButton("Смещение атомов(Оптический)")
        self.atom_displacements_button.setFixedSize(int(screen_geometry.width() * 0.15), 40)
        self.atom_displacements_button.setCheckable(True)
        self.atom_displacements_button.clicked.connect(self.update_plot)
        self.interface_layout.addWidget(self.atom_displacements_button)

        self.color_button = QPushButton("Поменять цвет фона")
        self.color_button.setFixedSize(int(screen_geometry.width() * 0.15), 40)
        self.color_button.clicked.connect(self.change_background_color)
        self.interface_layout.addWidget(self.color_button)

        self.layout.addLayout(self.interface_layout)
        self.update_plot()

    def update_plot(self):
        try:
            k0 = self.k0_slider.value() / 100.0
            delta_k = float(self.delta_k_input.text()) if self.delta_k_input.text() else 1.0
            b = float(self.b_input.text()) if self.b_input.text() else 1.0  # Default value for b is 1
            z0 = float(self.z0_input.text()) if self.z0_input.text() else -1.0
            z1 = float(self.z1_input.text()) if self.z1_input.text() else 1.0
            t = self.time_slider.value() * 2 * np.pi / 100.0
            k_b = float(self.k_b_input.text()) if self.k_b_input.text() else -1.0
            k_e = float(self.k_e_input.text()) if self.k_e_input.text() else 1.0
            k = np.arange(k_b, k_e, 0.01)

            # Проверяем тип фонона
            if self.equation_state == 0:
                result = self.acoustic_phonon(k, t, k0, delta_k)
                equation_label = "Акустический фонон"
            else:
                result = self.optical_phonon(k, t, k0, delta_k, b)
                equation_label = "Оптический фонон"

            # Проверяем кнопку волнового пакета
            if self.wave_packet_button.isChecked():
                # Если кнопка нажата, то рассматриваем волновой пакет
                step = 0.025
                z = np.arange(-1, 1, step)
                # Выбрал границы реального пространства <= 30, ибо интегралы долго считаются для больших значений
                if z1 > z0 and (z1 - z0) <= 30:
                    z = np.arange(z0, z1, step)

                def func(_k, _z, _t):
                    if self.equation_state:
                        return np.exp(-1 * (_k - k0) ** 2 / delta_k ** 2)\
                               * np.cos(np.sin(np.pi * _k / 2) * _t)\
                               * np.exp(-1j * _k * _z * _t)
                    else:
                        return np.exp(-1 * (_k - k0) ** 2 / delta_k ** 2)\
                               * np.cos((1 - b * np.cos(np.pi * _k)) * _t)\
                               * np.exp(-1j * _k * _z * _t)

                res, _ = quad_vec(lambda k: func(k, z, t), z0, z1)
                self.plot_widget.clear()
                self.plot_widget.showGrid(True, True)
                self.plot_widget.addLegend(brush=pg.mkBrush(220, 220, 220, 100))
                if self.show_real_button.isChecked():
                    self.plot_widget.plot(z, res.real, name='Real', pen=pg.mkPen('r', width=3)
                                          if self.show_real_button.isChecked() else None)

                if self.show_imaginary_button.isChecked():
                    self.plot_widget.plot(z, res.imag, name='Imaginary', pen=pg.mkPen('g', width=3)
                                          if self.show_imaginary_button.isChecked() else None)

                if self.show_amplitude_button.isChecked():
                    self.plot_widget.plot(z, np.abs(res) ** 2, name='Ampl', pen=pg.mkPen('b', width=3)
                                          if self.show_amplitude_button.isChecked() else None)
                self.plot_widget.setLabel('bottom', "z")
                self.plot_widget.setLabel('left', "A(z, t)")
            elif self.atom_displacements_button.isChecked():
                equation_label = "Оптический фонон"
                if delta_k < 0.05:
                    delta_k = 0.05
                lattice_spacing = 0.1
                i_values = np.arange(-80, 80, 1)
                z_values = lattice_spacing * i_values
                x_values = (-1) ** np.abs(i_values) * np.exp(-lattice_spacing ** 2 * delta_k * 2 * i_values ** 2)

                self.plot_widget.clear()
                self.plot_widget.showGrid(True, True)

                # Добавляю огибающие
                pos_envelope_x = [z_values[i] for i in range(len(z_values)) if x_values[i] >= 0]
                pos_envelope_y = [max(0, x_values[i]) for i in range(len(z_values)) if x_values[i] >= 0]
                neg_envelope_x = [z_values[i] for i in range(len(z_values)) if x_values[i] < 0]
                neg_envelope_y = [min(0, x_values[i]) for i in range(len(z_values)) if x_values[i] < 0]

                pos_envelope_curve = pg.PlotCurveItem(x=pos_envelope_x, y=pos_envelope_y,
                                                      pen=pg.mkPen('g', width=2))
                neg_envelope_curve = pg.PlotCurveItem(x=neg_envelope_x, y=neg_envelope_y,
                                                      pen=pg.mkPen('r', width=2))

                self.plot_widget.addItem(pos_envelope_curve)
                self.plot_widget.addItem(neg_envelope_curve)

                # Добавляю вертикали для каждого из значений
                for i in range(len(z_values)):
                    line = pg.PlotCurveItem(x=[z_values[i], z_values[i]], y=[0, x_values[i]],
                                            pen=pg.mkPen('b', width=1))
                    self.plot_widget.addItem(line)
                dots = pg.ScatterPlotItem(x=z_values, y=x_values, pen=None, symbol='o',
                                          brush=pg.mkBrush('orange'), size=7)
                self.plot_widget.addItem(dots)

                self.plot_widget.setLabel('bottom', "z")
                self.plot_widget.setLabel('left', "x(z)")
            else:
                self.plot_widget.clear()
                self.plot_widget.showGrid(True, True)
                self.plot_widget.plot(k, result, pen=pg.mkPen('b', width=3))
                self.plot_widget.setLabel('bottom', "k")
                self.plot_widget.setLabel('left', "a(k, t)")
            self.plot_widget.setTitle(equation_label)
            # Меняю значения при изменении
            self.k0_value_label.setText(f"Текущее значение K0: {k0:.2f}")
            self.delta_k_value_label.setText(f"Текущее значение Δk: {delta_k:.2f}")
            self.b_value_label.setText(f"Текущее значение b: {b:.2f}")
            self.z0_value_label.setText(f"Текущее значение z0: {z0:.2f}")
            self.z1_value_label.setText(f"Текущее значение z1: {z1:.2f}")
            self.time_value_label.setText(f"Текущее значение Время: {self.time_slider.value() / 100.0:.2f}")

        except ValueError:
            pass

    def switch_equation(self):
        # Меняю тип фонона
        self.equation_state = 1 - self.equation_state
        self.update_plot()

        if self.equation_state == 0:
            self.equation_button.setText("Поменять тип фонона (Акустический)")
        else:
            self.equation_button.setText("Поменять тип фонона (Оптический)")

    @staticmethod
    def acoustic_phonon(k, t, k0, delta_k):
        return np.exp(-1 * (k - k0) ** 2 / (delta_k ** 2 + 1e-50)) * np.cos(np.sin(np.pi / 2 * k) * t)

    @staticmethod
    def optical_phonon(k, t, k0, delta_k, b):
        return np.exp(-1 * (k - k0) ** 2 / (delta_k ** 2 + 1e-50)) * np.cos((1 - b * np.cos(np.pi * k)) * t)

    def change_background_color(self):
        color_dialog = QColorDialog(self)
        color = color_dialog.getColor()

        if color.isValid():
            # Меняю цвет фона
            color_name = color.name()
            self.plot_widget.setBackground(color_name)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
