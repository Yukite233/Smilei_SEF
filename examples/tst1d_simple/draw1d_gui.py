import sys
import random
import matplotlib
import h5py as h5
import numpy as np
from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QVBoxLayout, QHBoxLayout, QSizePolicy, QMessageBox, QWidget, QSplitter


from hdftreemodel import *
from mycanvas import *


class ApplicationWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("application main window")

        self.file_menu = QMenu('&File', self)
        self.file_menu.addAction('&Quit', self.fileQuit,
                                 QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.help_menu = QMenu('&Help', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)

        self.help_menu.addAction('&About', self.about)



        self.filename = "Fields_global.h5"
        self.f=h5.File(self.filename)
        dset = self.f["/1d_global/potential"]

        print "Attention: the dataset must be 4 dimensions"

        #> the main layout: left part (for the hdf5 tree view) and right part (for the plot area)
        self.main_widget = QWidget(self)

        #> if not putting splitter into a layout, the widgets in splitter do not fill the main windows
        #> (may exceed the app windows, so that the figures are partially shown ).
        layout = QVBoxLayout(self.main_widget)
        hSplitter = QSplitter(self.main_widget)
        layout.addWidget(hSplitter)


        #> the left part: the hdf5 tree view
        h5tree = QWidget()
        treeview = QTreeView(h5tree)
        self.model = HDFTreeModel([])
        self.model.openFile(self.filename, 'r+')
        treeview.setModel(self.model)

        treeview.doubleClicked.connect(self.redraw)

        hSplitter.addWidget(treeview)


        #> the right part: the plot area
        plotArea = QWidget(self.main_widget)
        sizePolicy = QSizePolicy();
        sizePolicy.setHorizontalPolicy(QSizePolicy.Expanding);
        sizePolicy.setVerticalPolicy(QSizePolicy.Expanding);
        plotArea.setSizePolicy(sizePolicy);

        hSplitter.addWidget(plotArea)



        plotVboxlayout = QVBoxLayout(plotArea)
        self.sc = MyStaticMplCanvas(self.main_widget, dset, width=5, height=4, dpi=100)
        self.dc = MyDynamicMplCanvas(self.main_widget, dset, width=5, height=4, dpi=100)
        plotVboxlayout.addWidget(self.sc)
        plotVboxlayout.addWidget(self.dc)


        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        self.statusBar().showMessage("All hail matplotlib!", 2000)

    def redraw(self, index):
    	t = 0
    	item = self.model.getItem(index)
    	if (item is not None) and item.isDataset():
    		dset = item.h5node
    		self.sc.compute_initial_figure(dset, t)
    		self.dc.compute_initial_figure(dset, t)
    		self.sc.draw()

    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()

    def about(self):
        QMessageBox.about(self, "About",
"""embedding_in_qt5.py example
Copyright 2015 BoxControL

This program is a simple example of a Qt5 application embedding matplotlib
canvases. It is base on example from matplolib documentation, and initially was
developed from Florent Rougon and Darren Dale.

http://matplotlib.org/examples/user_interfaces/embedding_in_qt4.html

It may be used and modified with no restriction; raw copies as well as
modified versions may be distributed without limitation."""
)

if __name__ == '__main__':
    app = QApplication(sys.argv)

    aw = ApplicationWindow()
    aw.setWindowTitle("PyQt5 Matplot Example")
    aw.show()
    #sys.exit(qApp.exec_())
    app.exec_()
