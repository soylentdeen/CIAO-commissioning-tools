import scipy
import numpy
import VLTTools

ciao = VLTTools.VLTConnection(simulate=False)

ciao.replaceSlopesWithCurrent()
