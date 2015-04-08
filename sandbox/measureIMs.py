import scipy
import numpy
import pyfits
import VLTTools
import SPARTATools
import os
import glob
import time

ciao = VLTTools.VLTConnection(simulate=False)

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

refpixelsFiles = glob.glob(os.path.expanduser('~')+'/data/'+datadir+'RP2*.fits')

for f in refpixelsFiles:
    junk = f.split('_')
    x =junk[-3]
    y = junk[-2]
    print("%s   %s" % (x, y))
    ciao.updateReferenceSlopes(f)
    print("Updated the reference Slopes")
    time.sleep(1.0)
    ciao.measureNewTTRefPositions("TWHydra")
    
    print("Updated the Tip and Tilt Reference Positions")
    time.sleep(5.0)
    ciao.setup_HOIM()
    ciao.measure_HOIM(config=True)
    #raw_input("Record Interaction Matrix, Press Enter when Done")
    ciao.saveMap(mapname="HORecnCalibrat.RESULT_IM", filename=datadir+"IM_"+x+"_"+y+"_.fits")


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

