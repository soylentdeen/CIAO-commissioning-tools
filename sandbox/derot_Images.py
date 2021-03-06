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
datadir = "/diska/data/SPARTA/2015-04-10/Derotator_7/"
ciao = VLTTools.VLTConnection(simulate=False, datapath=datadir)

#logfile = open(os.path.expanduser('~')+'/data/'+datadir+'logfile.dat', 'w')

short = True

if short:
    angles = numpy.arange(10)*40.0
else:
    #angles = numpy.arange(38)*9.73
    angles = numpy.arange(25)*15.0

#ciao.setupFieldLens()
i = 0
for angle in angles:

    ciao.moveDerotator(angle)
    print 'Moved to ', angle

    ciao.measurePixelFrames("derot_circbuff_"+str(i)+"_")
    #ciao.setup_HOIM(cycles=1)
    #ciao.measure_HOIM(config=True)
    #ciao.saveMap(mapname="HORecnCalibrat.RESULT_IM", filename="IM_"+str(i)+"_.fits")
    i += 1


#logfile.close()
#
