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

def pseudoinv(mat, nFilt, calm=False):
   DtD = numpy.dot(mat.T, mat)
   if calm:
       avdiag = numpy.average(numpy.diag(DtD))
       for i in range(40,60):
           DtD[i,i] += avdiag*0.05


class modalBasis ( object ):
   def __init__(self, HOIM, TTIM, delta, fullNoll, Z2S, nFilt):
       self.HOIM = HOIM
       self.TTIM = TTIM
       self.delta = delta
       self.fullNoll = fullNoll
       self.Z2S = Z2S
       self.nFilt = nFilt
       self.nHOAct = 60   # Number of High-Order Actuators
       self.nTTAct = 2    # Number of Tip/Tilt Actuators
       self.nSubap = 68   # Number of Sub-Apertures
       self.createModalBasis()

   def createModalBasis(self):
       # Generalized inverse of the High-Order Interaction Matrix
       HO_inv = pseudoinv( self.HOIM, self.nFilt)

       # Computation of set of DM voltages that reproduce a pure Tilt/Tip
       # WARNING: These tip/tilts may contain some piston
       TT2HO = numpy.dot(HO_inv, self.TTIM)

       # Computation of PISTON mode as the least visible with highest
       # variance
       lam, eigenvect = diagonalisation(numpy.dot(self.HOIM.T, self.HOIM))
       self.pistonMode = eigenvect[:,-1]
       self.pistonProj = numpy.dot(self.delta, self.pistonMode)

       self.TT2HO = self.filterOutPiston(TT2HO) 

       # Compute the reverse projection HO2TT
       self.HO2TT = numpy.dot(TT2HO, pseudoinv(TT2HO,0))

       # Compute modes as KL in the measurement space
       turb, eigss = self.KLbasisSlopesSpace()
       self.nUseful = self.nHOAct - self.nFilt
       modes = numpy.dot(HO_inv, eigss)[:,2:self.nUseful]

       # Make all these modes orthogonal to piston
       modes = self.filterOutPiston(modes) 

       # Normalize the basis to the same value as TT2HO
       var = numpy.sum(modes*numpy.dot(self.delta, modes), axis=0)
       varTT2HO = numpy.average( numpy.sum(self.TT2HO*numpy.dot(self.delta,
                                 self.TT2HO), axis=0) )
       self.modes = numpy.dot(modes, numpy.diag(numpy.sqrt(varTT2HO/var)) )

       # Compute Modes to Voltages
       self.M2V = numpy.zeros((self.nHOACt+self.nTTAct, self.nUseful))
       self.M2V[:self.nHOAct,selfnTTAct: ] = self.modes
       self.M2V[self.nHOAct, 0] = 1.
       self.M2V[self.nHOAct+1, 1] = 1.

       # Make a single interaction matrix
       self.IM = numpy.zeros((self.nSubap, self.nHOAct+self.nTTAct))
       self.IM[:,0:self.nHOAct] = self.HOIM
       self.IM[:,self.nHOAct:self.nHOAct+self.nTTAct] = self.TTIM

       # Compute the command Matri9x
       self.M2S = numpy.dot(self.IM, self.M2V)  # Modal Interaction Matrix
       self.S2M = pseudoinv(self.M2S, 0)        # Modal Control Matrix

       #            AWF
       # Compute the projection of a voltage vector onto the command space
       self.M2VHO = numpy.zeros((self.nHOAct, self.nUseful))
       self.M2VHO[:,0:2] = self.TT2HO
       self.M2VHO[:,self.nTTAct:self.nTTAct] = self.M2V[:self.nHOAct,self.nTTAct:]
       projAWF = numpy.identity(self.nHOAct] - numpy.dot(self.M2VHO,
                      pseudoinv(self.M2VHO, 0))
       lam, mod = diagonalisation(projAWF)
       self.AWFbasis = numpy.dot(mod[:,0:self.nFilt], 
                      numpy.diag(1.0/numpy.sqrt(lam[0:self.nFilt])))

       self.createSMAbasis()

   def createSMAbasis(self):
       m = self.filterOutPiston(numpy.identity(self.nHOAct), self.pistonMode,
                                self.pistonProj)
       lam, mo = diagonalisation(numpy.dot(m.T, numpy.dot(self.delta, m)))
       mo = numpy.dot(m, mo)

       self.SMAbasis = numpy.zeros(delta.shape)
       self.SMAbasis[:,0] = self.pistonMode
       self.SMAbasis[:,1:] = mo[:,:,-1]


