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
datadir = "/diska/data/SPARTA/2015-05-20/Derotator_2/"
ciao = VLTTools.VLTConnection(simulate=False, datapath=datadir)

#logfile = open(os.path.expanduser('~')+'/data/'+datadir+'logfile.dat', 'w')

short = False

ciao.moveDerotator(0.0)
ciao.set_Tip(-0.017)
ciao.set_Tilt(0.03)

ciao.measureNewTTRefPositions("TWHydra")

TT = [ciao.get_Tip(), ciao.get_Tilt()]

ciao.setup_HOIM(cycles=5)
ciao.measure_HOIM(config=True)
ciao.get_InteractionMatrices()
ciao.calc_CommandMatrix(nFiltModes=20)

ciao.set_HO_gain(-0.01)

ciao.measureNewHORefPositions("BetaPic")


print TT

if short:
    angles = numpy.arange(10)*40.0
else:
    #angles = numpy.arange(38)*9.73
    #angles = numpy.arange(25)*15.0
    angles = numpy.arange(26)*14.4

#ciao.setupFieldLens()
i = 0
for angle in angles:

    ciao.moveDerotator(angle)
    print 'Moved to ', angle

    ciao.setup_TTIM(cycles=3)
    ciao.measure_TTIM(config=True)
    ciao.get_InteractionMatrices()
    ciao.calc_commandMatrix(nFiltModes=20)

    ciao.measureNewTTRefPositions("TWHydra")

    ciao.measureCircularBuffer("derot_circbuff_"+str(i)+"_")
    ciao.setup_HOIM(cycles=3)
    ciao.measure_HOIM(config=True)
    ciao.saveMap(mapname="HORecnCalibrat.RESULT_IM", filename="IM_"+str(i)+"_.fits")
    i += 1


#logfile.close()
#
