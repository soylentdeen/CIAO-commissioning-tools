import scipy
import numpy
import pyfits
import VLTTools
import SPARTATools
import os
import glob
import time

datdir = "/diska/data/SPARTA/2015-05-19/PupilConjugation_2/"
ciao = VLTTools.VLTConnection(simulate=False, datapath=datdir)

"""
Which variables do we want to vary?

AMPLITUDE:  0.01 - 0.15
NOISE_THRESHOLD: 0.01 - 0.1
SKIP_TIME:  0.01 - 0.05
WAVE_PERIOD: 0.1 - 0.5
CYCLES: 1 - 5

"""
i = 0
datadir = "HOIM_test_data/"

logfile = open(os.path.expanduser('~')+'/data/'+datadir+'logfile.dat', 'w')

#refpixelsFiles = glob.glob(os.path.expanduser('~')+'/data/'+datadir+'RP2*.fits')

rpx = [0.0, -0.35, -0.35, -0.5, 0.0, 0.35, 0.35, 0.5]
rpy = [-0.5, -0.35, 0.35, 0.0, 0.5, -0.35, 0.35, 0.0]

#for f in refpixelsFiles:
for x, y in zip(rpx, rpy):
    print("%s   %s" % (x, y))
    ciao.updateRefSlopes(x+3.5, y+3.5)
    print("Updated the reference Slopes")
    time.sleep(1.0)
    ciao.measureNewTTRefPositions("TWHydra")
    
    print("Updated the Tip and Tilt Reference Positions")
    time.sleep(5.0)
    ciao.setup_HOIM()
    ciao.measure_HOIM(config=True)
    #raw_input("Record Interaction Matrix, Press Enter when Done")
    ciao.saveMap(mapname="HORecnCalibrat.RESULT_IM", filename="IM_"+str(x)+"_"+str(y)+"_.fits")


"""
                    print("%i %.2f %.2f %.2f %.2f %.2f\n" % 
                                  (i, amplitude, noise, skip, period, cycles))
                    ciao.setup_HOIM(amplitude=amplitude, noise=noise, skip=skip,
                                   period=period, cycles=cycles)
                    ciao.measure_HOIM(config=True)
                    #ciao.get_HOIM()
                    ciao.saveMap(mapname="HORecnCalibrat.RESULT_IM", filename=datadir+"HOIM_"+str(i)+".fits")
                    i += 1
                    logfile.write("%i %.2f %.2f %.2f %.2f %.2f\n" % 
                                  (i, amplitude, noise, skip, period, cycles))
                    

"""
logfile.close()

