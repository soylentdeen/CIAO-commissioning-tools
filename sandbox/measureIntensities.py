import scipy
import pyfits
import numpy
import VLTTools

ciao = VLTTools.VLTConnection(simulate=False)

ciao.averageIntensities()
