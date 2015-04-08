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
datadir = "/diska/data/SPARTA/2015-04-08/Parabola_1/"
ciao = VLTTools.VLTConnection(simulate=False, datapath=datadir)

#logfile = open(os.path.expanduser('~')+'/data/'+datadir+'logfile.dat', 'w')


#ciao.setupFieldLens()

tipRef = 505820
tiltRef = 1772160

i = 0
for angle in [0, 45]:
    ciao.moveDerotator(angle)
    print 'Moved to ', angle
    for TT in zip([0, 1000, -1000, 0, 0], [0, 0, 0, 1000, -1000]):
        Tip = tipRef+TT[0]
        Tilt = tiltRef+TT[1]
        ciao.movePM(Tip, Tilt)
        ciao.measureCircularBuffer("PM_circbuff_"+str(i)+"_"+str(TT[0])+"_"+str(TT[1])+"_"+str(angle)+"_")
        ciao.setup_HOIM()
        ciao.measure_HOIM(config=True)
        ciao.saveMap(mapname="HORecnCalibrat.RESULT_IM", filename="IM_"+str(i)+"_"+str(TT[0])+"_"+str(TT[1])+"_"+str(angle)+"_.fits")
        i+= 1


#logfile.close()
#
