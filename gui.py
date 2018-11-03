#!/usr/bin/env python3

import sys
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QPushButton, QLabel, QTextEdit, QFileDialog, QApplication
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import os
import numpy as np


class MainApp(QWidget):
    def __init__(self, parent=None):
        super(MainApp, self).__init__(parent)
        self.init_ui()

    def init_ui(self):
        # Add button that loads input data file, does inversion and displays results
        layout = QVBoxLayout()
        self.btn = QPushButton("Load .dat input file")
        self.btn.clicked.connect(self.do_inversion)
        layout.addWidget(self.btn)

        # Add button that allows to select output density file to plot
        self.btn2 = QPushButton("Select density to plot")
        self.btn2.clicked.connect(self.plot_inversion_results)
        layout.addWidget(self.btn2)

        # Add matplotlib canvas
        self.p = PlotCanvas(self, width=8, height=6)
        layout.addWidget(self.p)

        self.setLayout(layout)
        self.setWindowTitle("Borehole Gravimetry Inversion")
        self.setGeometry(50, 50, 1000, 1000)

    def do_inversion(self):
        fname = self.get_dat_input_file()
        #fname = fname.replace("/", "\\")
        self.call_inversion_function(fname)
        result_fname = fname.replace(".dat", "_density.dat")
        depth, dens = self.load_density_from_file(result_fname)
        self.p.plot(depth, dens)

    def plot_inversion_results(self):
        fname = self.get_dat_result_file()
        depth, dens = self.load_density_from_file(fname)
        self.p.plot(depth, dens)

    def get_dat_input_file(self):
        return self.get_dat_file("Open .dat input file")

    def get_dat_result_file(self):
        return self.get_dat_file("Open .dat result file")

    def get_dat_file(self, title):
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
            depth, density = np.loadtxt(filepath, delimiter="\t", unpack=True)
        return depth, density


class PlotCanvas(FigureCanvas):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        #FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.toolbar = NavigationToolbar(self, self)

    def plot(self, depth, density):
        ax = self.figure.add_subplot(111)
        # clear previous ax content, so replotting works
        ax.clear()
        ax.plot(depth, density, 'r-')
        ax.set_title('Density model')
        ax.set_xlabel("Depth (m)")
        ax.set_ylabel("Density (g/cm³)")
        self.draw()

def main():
    app = QApplication(sys.argv)
    plot_application = MainApp()
    plot_application.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()