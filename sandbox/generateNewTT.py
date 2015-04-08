import scipy
import pyfits
import numpy
import VLTTools
import sys

ciao = VLTTools.VLTConnection(simulate=False)

ciao.averageTTPositions(sys.argv[1])
