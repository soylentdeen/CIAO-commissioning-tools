import scipy
import pyfits
import numpy
import VLTTools
import sys

ciao = VLTTools.VLTConnection(simulate=False)

ciao.calc_CommandMatrix()
ciao.applyZernike([1.0, 0.0, 0.0, 0.0])
