import scipy
import numpy
import pyfits
import VLTTools
import SPARTATools
import os

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

for amplitude in numpy.linspace(0.05, 0.25, num=5):
    for noise in numpy.linspace(0.01, 0.1, num=5):
        for skip in numpy.linspace(0.01, 0.05, num=5):
            for period in numpy.linspace(0.1, 0.5, num=5):
                for cycles in numpy.arange(1)+1:
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
                    


logfile.close()
