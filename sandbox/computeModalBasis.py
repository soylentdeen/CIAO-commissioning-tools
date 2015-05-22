import scipy
import numpy
import matplotlib.pyplot as pyplot
import pyfits
import VLTTools

ciao = VLTTools.VLTConnection(simulate=False)

ciao.get_InteractionMatrices()
ciao.calc_CommandMatrix(nFiltModes=20)
print "This is where we will Compute the Modal Basis from the IMs"
