import scipy
import numpy
import pyfits
import VLTTools
import SPARTATools
import os
import glob
import time

ciao = VLTTools.VLTConnection(simulate=False)

i = 0
datadir = "Derotator_test_data/"

logfile = open(os.path.expanduser('~')+'/data/'+datadir+'logfile.dat', 'w')


#ciao.setupFieldLens()

for angle in range(0, 400, 40):

    ciao.moveDerotator(angle)
    print 'Moved to ', angle

    time.sleep(0.5)
    ciao.setup_HOIM()
    ciao.measure_HOIM(config=True)
    ciao.saveMap(mapname="HORecnCalibrat.RESULT_IM", filename=datadir+"IM_"+str(angle)+"_.fits")


logfile.close()

