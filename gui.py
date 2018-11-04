#!/usr/bin/env python3

import sys
from PyQt5.QtWidgets import QStatusBar, QMenuBar, QWidget, QDesktopWidget, QVBoxLayout, qApp, QFileDialog, QApplication, \
    QSizePolicy, QAction, QMessageBox, QMainWindow
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import os
import numpy as np


class MainApp(QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainApp, self).__init__(*args, **kwargs)
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        self.main_widget = QWidget(self)

        #create action to plot measurement data
        plot_measurement_data_action = QAction("Plot &measurement data", self)
        plot_measurement_data_action.setStatusTip("Select a file containing measurements to plot them")
        plot_measurement_data_action.triggered.connect(self.plot_measurement_data)

        # create action to load input file for inversion and do inversion
        load_dat_file_action = QAction("&Load .dat input", self)
        load_dat_file_action.setStatusTip("Select a file containing measurements for inversion")
        load_dat_file_action.triggered.connect(self.do_inversion)
        load_dat_file_action.showStatusText(self)

        # create action to plot density model
        plot_density_action = QAction("Plot &density model", self)
        plot_density_action.setStatusTip("Select a file containing density inversion to plot the model")
        plot_density_action.triggered.connect(self.plot_inversion_results)

        #create quit action, closes application
        quit_action = QAction("&Quit", self)
        quit_action.setStatusTip("Leave application")
        quit_action.triggered.connect(qApp.quit)

        #create menubar and attach actions
        self.menubar = QMenuBar(self)
        file_menu = self.menubar.addMenu("&File")
        file_menu.addAction(plot_measurement_data_action)
        file_menu.addAction(load_dat_file_action)
        file_menu.addAction(plot_density_action)
        file_menu.addAction(quit_action)
        layout.addWidget(self.menubar)

        # Add matplotlib canvas
        self.p = PlotCanvas(self, width=8, height=6)
        plot_toolbar = self.p.create_toolbar(self)
        layout.addWidget(plot_toolbar)
        layout.addWidget(self.p)

        # Add statusbar that displays current information
        self.setStatusBar(QStatusBar(self))
        self.statusBar().showMessage("Select a file...")


        self.main_widget.setLayout(layout)
        self.setCentralWidget(self.main_widget)
        self.setWindowTitle("Borehole Gravimetry Inversion")
        self.setGeometry(50, 50, 1000, 1000)
        self.center_window()

    def center_window(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def do_inversion(self):
        fname = self.get_dat_input_filepath()
        if fname is None:
            return
        if "density" in fname.split("/")[-1]:
            QMessageBox.question(self, "Wrong file selected!", "Select a file without density in filename", QMessageBox.Ok)
            self.statusBar().showMessage("Invalid file")
            return
        self.call_inversion_function(fname)
        result_fname = fname.replace(".dat", "_density.dat")
        depth, dens = self.load_density_from_file(result_fname)
        self.p.plot(depth, dens, title="Density model", ylabel='Density (g/cm³)', linestyle='r-')
        self.setWindowTitle("Density model for {}".format(fname))

    def plot_inversion_results(self):
        fname = self.get_dat_result_filepath()
        if fname is None:
            return
        if "density" not in fname:
            QMessageBox.question(self, "Wrong file selected!", "Select a file with density in filename", QMessageBox.Ok)
            self.statusBar().showMessage("Invalid file")
            return
        depth, dens = self.load_density_from_file(fname)
        self.p.plot(depth, dens, title='Density model', ylabel='Density (g/cm³)', linestyle='r-')
        self.setWindowTitle("Density model for {}".format(fname))

    def plot_measurement_data(self):
        fname = self.get_dat_input_filepath()
        if fname is None:
            return
        if "density" in fname.split("/")[-1]:
            QMessageBox.question(self, "Wrong file selected!", "Select a file without density in filename", QMessageBox.Ok, QMessageBox.Ok)
            self.statusBar().showMessage("Invalid file")
            return
        depth, grav = self.load_density_from_file(fname)
        self.p.plot(depth, grav, title="Gravity measurement", ylabel="Gravity (mGal)", linestyle='.')
        self.setWindowTitle("Measurement data for {}".format(fname))

    def get_dat_input_filepath(self):
        return self.get_dat_filepath("Open .dat input file")

    def get_dat_result_filepath(self):
        return self.get_dat_filepath("Open .dat result file")

    def get_dat_filepath(self, title):
        dlg = QFileDialog()
        # dont allow creating new files
        dlg.setFileMode(QFileDialog.ExistingFile)
        # only show dat files
        dlg.setNameFilter("*.dat")
        dlg.setWindowTitle(title)
        # start in directory of this script
        dlg.setDirectory(__file__.rpartition("/")[0])
        if dlg.exec_():
            filenames = dlg.selectedFiles()
            return filenames[0]

    def call_inversion_function(self, filepath):
        """
        Call inversion function on file given in filepath
        :param filepath:
        :type filepath:
        :return:
        :rtype:
        """
        progname = "Programm.exe" if os.name == 'nt' else "Programm"
        prog_path = os.path.join(".", "cmake-build-debug", progname)
        os.system("{} {}".format(prog_path, filepath))

    def load_density_from_file(self, filepath):
        with open(filepath, 'r') as f:
            depth, density = np.loadtxt(filepath, unpack=True)
        return depth, density


class PlotCanvas(FigureCanvas):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.fig.set_tight_layout(True)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def create_toolbar(self, parent):
        return NavigationToolbar(self, parent)

    def plot(self, x_values, y_values, title, ylabel, *args, **kwargs):
        ax = self.fig.gca()
        # clear previous ax content, so replotting works
        ax.clear()
        ax.plot(x_values, y_values, linestyle)
        ax.set_title(title)
        ax.set_xlabel("Depth (m)")
        ax.set_ylabel(ylabel)
        self.draw()


def main():
    app = QApplication(sys.argv)
    plot_application = MainApp()
    plot_application.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()