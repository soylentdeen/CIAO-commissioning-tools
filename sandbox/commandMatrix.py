import scipy
import numpy
import matplotlib.pyplot as pyplot
import pyfits
import VLTTools

ciao = VLTTools.VLTConnection(simulate=False)
#ciao.measure_InteractionMatrices()
ciao.get_InteractionMatrices()
ciao.calc_CommandMatrix()
ciao.set_CommandMatrix()
ciao.save_CommandMatrixPlot()
#print "This is a test"
