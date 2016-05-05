import sys
import random
import matplotlib
import h5py as h5
import numpy as np
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5 import QtCore
from PyQt5.QtWidgets import QSizePolicy


class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, dset=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes1 = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes1.hold(False)

        self.dset = dset
        self.time = 0
        val = self.dset[...]
        self.compute_initial_figure(self.dset, self.time)

        #
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self, dset, time):
        pass

class MyStaticMplCanvas(MyMplCanvas):
    """Simple canvas with a sine plot."""
    def compute_initial_figure(self, dset, time):

        #> ==================================================================
        self.dset = dset
        self.time = time

        val = self.dset[...]
        val_1d = np.transpose(val[time, 0, 0, :])

        dx = 1.0
        nx = val_1d.size
        x = np.linspace(0,100.0,nx)

        self.axes1.plot(x, val_1d)
        self.axes1.axis([x.min(),x.max(),val_1d.min(),val_1d.max()])
        self.axes1.set_xlabel('x(mm)')
        self.axes1.set_ylabel('Charge density')





class MyDynamicMplCanvas(MyMplCanvas):
    """A canvas that updates itself every second with a new plot."""
    def __init__(self, *args, **kwargs):
        MyMplCanvas.__init__(self, *args, **kwargs)
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        timer.start(500)

    def compute_initial_figure(self, dset, time):
        #> ==================================================================
    	self.dset = dset
    	self.time = time


        val = self.dset[...]
        val_1d = np.transpose(val[self.time, 0, 0, :])

        dx = 1.0
        nx = val_1d.size
        x = np.linspace(0,100.0,nx)

	self.ntime = val[:, 0, 0, 0].size



        self.axes1.plot(x, val_1d)
        self.axes1.axis([x.min(),x.max(),val_1d.min(),val_1d.max()])
        self.axes1.set_xlabel('x(mm)')
        self.axes1.set_ylabel('Charge density')

    def update_figure(self):

	self.time = self.time + 1
	if self.time == self.ntime:
		self.time = 0

        val = self.dset[...]
        val_1d = np.transpose(val[self.time, 0, 0, :])

        dx = 1.0
        nx = val_1d.size
        x = np.linspace(0,100.0,nx)

        self.axes1.plot(x, val_1d)
        self.axes1.axis([x.min(),x.max(),val_1d.min(),val_1d.max()])
        self.axes1.set_xlabel('x(mm)')
        self.axes1.set_ylabel('Charge density')
        self.draw()
