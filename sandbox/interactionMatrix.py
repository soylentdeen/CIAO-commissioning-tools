import scipy
import numpy
import pyfits
import VLTTools
import SPARTATools

ciao = VLTTools.VLTConnection(simulate=False)
ciao.measure_InteractionMatrices()
#ciao.get_InteractionMatrices()
#ciao.calc_CommandMatrix()
#ciao.set_CommandMatrix()
