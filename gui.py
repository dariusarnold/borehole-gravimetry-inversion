#!/usr/bin/env python3

import functools
import os
import sys

import numpy as np
from PyQt5.QtWidgets import QStatusBar, QMenuBar, QWidget, QDesktopWidget, QVBoxLayout, qApp, QFileDialog, QApplication, \
    QSizePolicy, QAction, QMessageBox, QMainWindow, QInputDialog
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


def get_file_ending(filepath):
    """
    Return file ending of given file
    :param filepath:
    :return:
    """
    filename = filepath.split(os.sep)[-1]
    file_ending = filename.split(".")[-1] if "." in filename else ""
    return file_ending


class MainApp(QMainWindow):
    def __init__(self, *args, **kwargs):
        self.discretization_steps = 10000
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
        load_dat_file_action = QAction("&Invert for density model", self)
        load_dat_file_action.setStatusTip("Select a file containing measurements for inversion")
        load_dat_file_action.triggered.connect(self.do_inversion)
        load_dat_file_action.showStatusText(self)

        # create action to plot density model
        plot_inversion_results_action = QAction("Plot &inversion result", self)
        plot_inversion_results_action.setStatusTip("Select a file containing density inversion to plot the model")
        plot_inversion_results_action.triggered.connect(functools.partial(self.plot_inversion_results, None))

        # create quit action, closes application
        quit_action = QAction("&Quit", self)
        quit_action.setStatusTip("Leave application")
        quit_action.triggered.connect(qApp.quit)

        # simple settings action asks user for number of discretization steps
        get_steps_action = QAction("&Discretization steps", self)
        get_steps_action.setStatusTip("Enter number of steps to use for discretization")
        get_steps_action.triggered.connect(self.get_steps_from_user)

        #create menubar and attach actions
        self.menubar = QMenuBar(self)
        file_menu = self.menubar.addMenu("&File")
        settings_menu = self.menubar.addMenu("&Settings")
        settings_menu.addAction(get_steps_action)
        file_menu.addAction(plot_measurement_data_action)
        file_menu.addAction(load_dat_file_action)
        file_menu.addAction(plot_inversion_results_action)
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
        """
        Complete inversion sequence with file chooser to select input data
        """
        # Open file dialog to get input file
        fname = self.get_dat_input_filepath()
        if not isinstance(fname, str):
            return
        # check whether correct file was selected
        if get_file_ending(fname) != 'dat':
            QMessageBox.question(self, "Wrong file selected!", "Select a .dat file containing measurement data", QMessageBox.Ok)
            self.statusBar().showMessage("Invalid file")
            return
        self.call_inversion_function(fname)
        # transform filename since the inversion program changes file endings when saving results
        result_fname = fname.replace(".dat", ".dens")
        self.plot_inversion_results(result_fname)

    def plot_inversion_results(self, fname=None):
        """
        Plot inversion results from fname.
        :param fname: if not given, opens dialog to let user select one
        """
        if not isinstance(fname, str):
            # if no fname was given when calling the function, get one from user
            fname = self.get_dens_result_filepath()
            # None is returned when user pressed cancel or closed the filedialog
            if fname is None: return
        if get_file_ending(fname) != 'dens':
            QMessageBox.question(self, "Wrong file selected!", "Select a .dens file containing inversion results", QMessageBox.Ok)
            self.statusBar().showMessage("Invalid file")
            return
        depth, dens = self.load_density_from_file(fname)
        self.p.plot(depth, dens, 'Density model', 'Density (g/cm³)', 'r-')
        self.setWindowTitle("Density model for {}".format(fname))

    def plot_measurement_data(self):
        """
        Plot measurement data against depth. Data should be whitespace seperated in two columns, depth and data
        """
        fname = self.get_dat_input_filepath()
        if fname is None: return
        if get_file_ending(fname) != 'dat':
            QMessageBox.question(self, "Wrong file selected!", "Select a .dat file containing measurement data", QMessageBox.Ok, QMessageBox.Ok)
            self.statusBar().showMessage("Invalid file")
            return
        depth, grav = self.load_density_from_file(fname)
        self.p.plot(depth, grav, "Gravity measurement", "Gravity (mGal)", 'r+', markersize=10)
        self.setWindowTitle("Measurement data for {}".format(fname))

    def get_steps_from_user(self):
        result, _ = QInputDialog.getInt(self, "Steps settings", "Enter number of discretization steps", min=1, value=self.discretization_steps)
        self.discretization_steps = result

    def get_dat_input_filepath(self):
        return self.get_filepath("Open .dat input file", "*.dat")

    def get_dens_result_filepath(self):
        return self.get_filepath("Open .dens result file", "*.dens")

    def get_filepath(self, filedialog_title, file_extension_filter):
        """
        Get filepath from user with a file chooser dialog
        :param filedialog_title: Window title of file dialog
        :param file_extension_filter: string restricting the files shown. Use "*.txt" to only show .txt files
        :return: filepath selected from user
        """
        dlg = QFileDialog()
        # dont allow creating new files
        dlg.setFileMode(QFileDialog.ExistingFile)
        # only show files matching the filter
        dlg.setNameFilter(file_extension_filter)
        dlg.setWindowTitle(filedialog_title)
        # start in directory of this script
        dlg.setDirectory(__file__.rpartition("/")[0])
        if dlg.exec_():
            filenames = dlg.selectedFiles()
            return filenames[0]

    def call_inversion_function(self, filepath):
        """
        Call inversion function on file given in filepath
        :param filepath: path to file which contains input data for inversion
        :type filepath: str
        :return:
        :rtype:
        """
        progname = "Programm.exe" if os.name == 'nt' else "Programm"
        prog_path = os.path.join(".", "cmake-build-debug", progname)
        call_string = "{programm_path} {input_path} {discretization_steps}"
        # discretzation_steps: Number of discretization steps that are applied to discretize the resulting density distribution.
        call_string = call_string.format(programm_path=prog_path, input_path=filepath, discretization_steps=self.discretization_steps)
        print(call_string)
        os.system(call_string)

    def load_density_from_file(self, filepath):
        depth, density = np.loadtxt(filepath, unpack=True)
        return depth, density


class PlotCanvas(FigureCanvas):

    def __init__(self, parent=None, width=5, height=4, dpi=144):
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
        ax.plot(x_values, y_values, *args, **kwargs)
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