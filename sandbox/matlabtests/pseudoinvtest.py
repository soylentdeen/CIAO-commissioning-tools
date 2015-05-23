import scipy
from scipy.linalg import *
import numpy
import pyfits
import SPARTATools
import matplotlib.pyplot as pyplot


def pseudoInverse(A, numFilteredModes=50):
    dims = A.shape
    U,S,V = svd(A)
    D = 1.0/(S[0:-numFilteredModes])
    #S[-numFilteredModes+1:-1] = 0.0
    S[-numFilteredModes:] = 0.0
    newS = numpy.zeros([dims[0], dims[1]])
    I = [i for i in range(dims[1])]
    for i in range(len(D)):
        newS[i][i] = D[i]
    
    S = newS.copy()
    
    retval = scipy.matrix(V.T.dot(S.T.dot(U.T)), dtype=numpy.float32)
    
    return retval

datapath = '/home/deen/Data/GRAVITY/SPARTA_Data/'
HOIM = pyfits.getdata(datapath+'HOIM.fits')

nfilt = 50

A = SPARTATools.pseudoinv(HOIM, nfilt)
B = pseudoInverse(HOIM, nfilt)
