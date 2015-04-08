import scipy
import numpy
import pyfits
import VLTTools
import SPARTATools
import os
import glob
import time


i = 0
#datadir = "Derotator_test_data/"
datadir = "/diska/data/SPARTA/2015-04-08/Derotator_9/"
ciao = VLTTools.VLTConnection(simulate=False, datapath=datadir)

#logfile = open(os.path.expanduser('~')+'/data/'+datadir+'logfile.dat', 'w')


#ciao.setupFieldLens()
i = 0
for angle in numpy.arange(38)*9.73:

    ciao.moveDerotator(angle)
    print 'Moved to ', angle

    ciao.measureCircularBuffer("derot_circbuff_"+str(i)+"_")
    ciao.setup_HOIM()
    ciao.measure_HOIM(config=True)
    ciao.saveMap(mapname="HORecnCalibrat.RESULT_IM", filename="IM_"+str(i)+"_.fits")
    i += 1


#logfile.close()
#
