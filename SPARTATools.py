import scipy
import numpy
#import matplotlib.pyplot as pyplot
import pyfits
import os
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
   lam, u = diagonalisation(DtD)
   for i in range(lam.shape[0]-nFilt):
       u[:,i] /= numpy.sqrt(lam[i])
   if (nFilt>0):
       u[:,-nFilt:] = 0.0
   return numpy.dot(numpy.dot(u, u.T), mat.T)

def diagonalisation(A):
   """
   Input arg:
   <A> : square, symmetric matrix

   Output:
   tuple with sorted eigenvalues and correpsonding eigenvectors

   This function is equivalent to numpy.linalg.eigh() (=diagonalisation
   of a symmetric matrix), plust sorting the eigenvalues
   """

   (eigenvalue, M) = numpy.linalg.eigh(A)

   sortindex = numpy.argsort(eigenvalue)[::-1]
   eigenvalue = eigenvalue[sortindex]
   M = M[:,sortindex]
   return (eigenvalue, M)

def computeNewBestFlat(outfile):
   datadir = "/diska/data/SPARTA/"
   datafile = datadir+"Amarillo_7/Amarillo_7.fits"

   data = pyfits.getdata(datafile)
   avg = numpy.average(data.field(5), axis=0)
   hdu = pyfits.PrimaryHDU(avg)
   hdu.writeto(outfile, clobber=True)

def computeIntensities(outfile):
   datadir = "/diska/data/SPARTA/"
   datafile= datadir+'hbonnet_cm/hbonnet_cm.fits'

   data = pyfits.getdata(datafile)
   intensities = data.field(3)
   avg = numpy.average(intensities, axis=0)

   hdu = pyfits.PrimaryHDU(avg)
   hdu.writeto(outfile, clobber=True)

def computeAITRefSlopes(outfile):
   data = numpy.ones(136, dtype=numpy.float32)*3.5
   hdu = pyfits.PrimaryHDU(data)
   hdu.writeto(outfile, clobber=True)

def computeDisturbanceFrame(actNum, nFrames, filename, range=0.2, max=0.3):
   nAct = 60
   frame = numpy.zeros((nFrames, nAct), dtype=numpy.float32)
   disturbance = numpy.random.randn(nframes)*range
   disturbance[disturbance > max] = max
   disturbance[disturbance < -max] = -max
   frame[:,actNum] = disturbance
   hdu = pyfits.PrimaryHDU(frame)
   hdu.writeto(filename, clobber=True)

class modalBasis ( object ):
   def __init__(self, HOIM, TTIM, nFilt):
       self.datapath = os.path.dirname(__file__)+'/Data/'
       self.HOIM = HOIM
       self.TTIM = TTIM
       self.delta = pyfits.getdata(self.datapath+"delta_MACAO.fits")
       self.fullNoll = pyfits.getdata(self.datapath+"ZernikeCovar_Noll.fits")
       self.Z2S = pyfits.getdata(self.datapath+"Z2S_136s_119z.fits").T
       self.nFilt = nFilt
       self.nHOAct = 60   # Number of High-Order Actuators
       self.nTTAct = 2    # Number of Tip/Tilt Actuators
       self.nSubap = 136   # Number of Sub-Apertures * 2
       self.createModalBasis()

   def createModalBasis(self):
       # Generalized inverse of the High-Order Interaction Matrix
       self.HO_inv = pseudoinv( self.HOIM, self.nFilt)

       # Computation of set of DM voltages that reproduce a pure Tilt/Tip
       # WARNING: These tip/tilts may contain some piston
       TT2HO = numpy.dot(self.HO_inv, self.TTIM)

       # Computation of PISTON mode as the least visible with highest
       # variance
       lam, eigenvect = diagonalisation(numpy.dot(self.HOIM.T, self.HOIM))
       self.pistonMode = eigenvect[:,-1]
       self.pistonProj = numpy.dot(self.delta, self.pistonMode)

       self.TT2HO = self.filterOutPiston(TT2HO) 

       # Compute the reverse projection HO2TT
       #self.HO2TT = numpy.dot(TT2HO, pseudoinv(TT2HO,0))
       self.HO2TT = numpy.dot(pseudoinv(self.TTIM,0), self.HOIM)

       # Compute modes as KL in the measurement space
       turb, eigss = self.KLbasisSlopesSpace()
       self.nUseful = self.nHOAct - self.nFilt
       modes = numpy.dot(self.HO_inv, eigss)[:,2:self.nUseful]

       # Make all these modes orthogonal to piston
       modes = self.filterOutPiston(modes) 

       # Normalize the basis to the same value as TT2HO
       var = numpy.sum(modes*numpy.dot(self.delta, modes), axis=0)
       varTT2HO = numpy.average( numpy.sum(self.TT2HO*numpy.dot(self.delta,
                                 self.TT2HO), axis=0) )
       self.modes = numpy.dot(modes, numpy.diag(numpy.sqrt(varTT2HO/var)) )

       # Compute Modes to Voltages
       self.M2V = numpy.zeros((self.nHOAct+self.nTTAct, self.nUseful))
       self.M2V[:self.nHOAct,self.nTTAct: ] = self.modes
       self.M2V[self.nHOAct, 0] = 1.
       self.M2V[self.nHOAct+1, 1] = 1.

       # Make a single interaction matrix

       self.IM = numpy.zeros((self.nSubap, self.nHOAct+self.nTTAct))
       self.IM[:,0:self.nHOAct] = self.HOIM
       self.IM[:,self.nHOAct:self.nHOAct+self.nTTAct] = self.TTIM

       # Compute the command Matrix
       self.M2S = numpy.dot(self.IM, self.M2V)  # Modal Interaction Matrix
       self.S2M = pseudoinv(self.M2S, 0)        # Modal Control Matrix

       #            AWF
       # Compute the projection of a voltage vector onto the command space
       self.M2VHO = numpy.zeros((self.nHOAct, self.nUseful))
       self.M2VHO[:,0:2] = self.TT2HO
       self.M2VHO[:,self.nTTAct:] = self.M2V[:self.nHOAct,self.nTTAct:]
       projAWF = numpy.identity(self.nHOAct) - numpy.dot(self.M2VHO,
                      pseudoinv(self.M2VHO, 0))
       lam, mod = diagonalisation(projAWF)
       self.AWFbasis = numpy.dot(mod[:,0:self.nFilt], 
                      numpy.diag(1.0/numpy.sqrt(lam[0:self.nFilt])))

       print asdf
       self.createSMAbasis()
       self.computeSystemControlMatrix()

   def createSMAbasis(self):
       m = self.filterOutPiston(numpy.identity(self.nHOAct))
       lam, mo = diagonalisation(numpy.dot(m.T, numpy.dot(self.delta, m)))
       mo = numpy.dot(m, mo)

       SMAbasis = numpy.zeros(self.delta.shape)
       SMAbasis[:,0] = self.pistonMode
       SMAbasis[:,1:] = mo[:,:-1]

       self.SMAbasis = SMAbasis

   def filterOutPiston(self, modes):
       """
       Input args :
       <modes>      2D array of floats, is an array containing N modes, it is an array 60xN
       <pistonMode> 1D array of floats, is the piston mode defined on actuator voltages (60 components)
       <delta>      2D array of floats, is the geometric cov matrix of the DM (60x60), fixed.
       The function will suppress piston from all the modes of matrix <modes>.
       """

       # Compute the scalar product between pistonMode and (delta.pistonMode)
       # this is a projection of the piston mode on itself
       pnorm = numpy.dot(self.pistonMode, self.pistonProj)

       # compute the scalar products between each of the modes and
       # delta.pistonMode, this is a projection of each mode on piston
       proj = numpy.dot(modes.T, self.pistonProj)
       # normalize by the value of the projection of the piston on itself
       proj /= pnorm
       # now we know the amount of each piston in each mode, we need
       # to subtract it...
       (m,n) = modes.shape
       fmodes = numpy.zeros((m, n)) # create the ouput, filtered modes
       for i in range(n):
           fmodes[:,i] = modes[:,i] - self.pistonMode*proj[i]
       return fmodes

   def KLbasisSlopesSpace(self):
       """
       Input args:
       <fullNoll> : matrix (120, 120), given as an input (covariance
                    matrix of zernike modes)
       <Z2S> : matrix zernike to slopes
       <miaHO> : HODM interaction matrix
       <miaHO1> : filtered pseudo inverse of miaHO

       Outputs:
       The procedure returns
       - the eigenvalues of the KL modes of the
       measurements projected onto the commanded subspace of miaHO1
       - the modes
       """
       n = self.fullNoll.shape[0]
       ns, nz = self.Z2S.shape

       # see what's the common number of zernike of both
       nz = min(n, nz)
       # cut matrices to the right size
       fullNoll = self.fullNoll[0:nz, 0:nz]
       Z2S = self.Z2S[:, 0:nz]
       # compute covariance matrix of the slopes Css = Z2S * Czz * Z2S.T
       covSlopes = numpy.dot( Z2S, numpy.dot(fullNoll, Z2S.T))

       proj = numpy.dot(self.HOIM, self.HO_inv)
       covSlopes = numpy.dot(proj, numpy.dot(covSlopes, proj.T))

       turb, eigss = diagonalisation(covSlopes)
       return turb, eigss


   def computeSystemControlMatrix(self):
       """
       stores the command matrix of the system in self.CM
       """
       # Matrix-multiply modes * S2M
       # in a "Normal System"
       CM = numpy.dot(self.M2V, self.S2M)

       # Transform the tip-tilt correction matrix into a
       # DM-tilt correction matrix
       tiltHO = numpy.dot(self.TT2HO, CM[60:62,:])

       # Add this DM-tilt correction matrix to the DM one
       CM[0:60,:] += tiltHO

       self.CM = CM
