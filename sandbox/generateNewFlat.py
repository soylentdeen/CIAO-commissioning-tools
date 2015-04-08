import scipy
import pyfits
import numpy
import VLTTools
import sys

ciao = VLTTools.VLTConnection(simulate=False)

ciao.averageActuatorPositions(sys.argv[1])
