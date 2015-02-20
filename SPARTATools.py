import scipy
import numpy
#import matplotlib.pyplot as pyplot
import pyfits
from scipy.linalg import *

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

def calculateCommandMatrix(HOIM, TTIM, nFiltModes=20):
    A = scipy.matrix(HOIM.data)
    dims = A.shape
    U,S,V = svd(A)
    D = 1.0/(S[0:-nFiltModes])
    S[-nFiltModes:] = 0.0
    newS = numpy.zeros([dims[0], dims[1]])
    I = [i for i in range(dims[1])]
    for i in range(len(D)):
        newS[i][i] = D[i]

    HOCM = scipy.matrix(V.T.dot(newS.T.dot(U.T)), dtype=numpy.float32)
    
    singular_values = newS.diagonal()
    svs = singular_values[singular_values.nonzero()[0]]

    TTCM = scipy.matrix(TTIM.data).getI()
    controlMatrix = numpy.resize(HOCM, (62,136))
    controlMatrix[-2] = TTCM[0]
    controlMatrix[-1] = TTCM[1]
    return controlMatrix
