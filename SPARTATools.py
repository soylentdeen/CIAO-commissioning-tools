import scipy
import numpy
import matplotlib.pyplot as pyplot
import pyfits

class framesViewer( object ):
    def __init__(self):
        pyplot.ion()
        self.fig = pyplot.figure(0)
        self.fig.clear()
        self.ax = self.fig.add_axes([0.1, 0.1, 0.8, 0.8])
        self.canvas = self.ax.imshow(numpy.zeros([2,2]))
        self.colorbar = self.fig.colorbar(self.canvas)
    
    def loadFile(self, datafile):
        self.data = pyfits.getdata(datafile)
        self.nframes = self.data.size

    def showFrame(self, i):
        frame = self.data[i]
        nx = frame[1]
        ny = frame[2]
        self.ax.cla()
        rand = numpy.random.rand(nx, ny)*frame[3].reshape([nx, ny])
        self.canvas = self.ax.imshow(rand)
        self.colorbar.update_bruteforce(self.canvas)
        self.fig.canvas.draw()
