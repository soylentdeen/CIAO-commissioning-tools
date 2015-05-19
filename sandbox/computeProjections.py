import scipy
import numpy
import matplotlib.pyplot as pyplot
import pyfits
import VLTTools

ciao = VLTTools.VLTConnection(simulate=False)

HOIM = '/diska/home/ciaomgr/data/HORecnCalibrat.RESULT_IM.fits'
TTIM = '/diska/home/ciaomgr/data/TTRecnCalibrat.RESULT.IM.fits'

ciao.set_InteractionMatrices(HOIM, TTIM)
ciao.calc_CommandMatrix(nFiltModes=10)
ciao.modalBasis.calcS2Z()
print "This is where we will Compute the Modal Basis from the IMs"
