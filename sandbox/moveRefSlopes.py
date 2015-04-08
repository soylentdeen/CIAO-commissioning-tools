import scipy
import pyfits
import VLTTools
import SPARTATools

ciao = VLTTools.VLTConnection(simulate=False)

df = '/diska/home/ciaomgr/data/HOIM_test_data/RP2_0_0_.fits'

#ciao.updateReferenceSlopes(df)
ciao.updateRefSlopes(3.0, 2.0)
