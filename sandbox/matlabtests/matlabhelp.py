import pyfits
import numpy
import SPARTATools

datapath = '/home/deen/Data/GRAVITY/SPARTA_Data/'
delta = pyfits.getdata(datapath+'delta_junk.fits')
fullNoll = pyfits.getdata(datapath+'ZernikeCovar_Noll.fits')
HOIM = pyfits.getdata(datapath+'HOIM.fits')
TTIM = pyfits.getdata(datapath+'TTIM.fits')

nfilt = 50
nHOAct = 60
nTTAct = 2
nSubap = 136

HO_inv = SPARTATools.pseudoinv(HOIM, nfilt)

TT2HO = numpy.dot(HO_inv, TTIM)

lam, ev = SPARTATools.diagonalisation(numpy.dot(HOIM.T, HOIM))
pistonMode = ev[:,-1]
pistonProj = numpy.dot(delta, pistonMode)
