import scipy
import numpy
import matplotlib.pyplot as pyplot
import pyfits
import VLTTools
import time

datadir = "/diska/data/SPARTA/2015-06-03/DMModulation_2/"

fname = "Modulation_20150603.fits"

ciao = VLTTools.VLTConnection(simulate=False, datapath=datadir)
#ciao.disturb(disturbance=fname, baseName="Test")

short = True
#scratch = False
scratch = True

if scratch:
    ciao.moveDerotator(0.0)
    ciao.set_Tip(-0.017)
    ciao.set_Tilt(0.03)

    ciao.setup_TTIM(cycles=3)
    ciao.measure_TTIM(config=True)
    ciao.setup_HOIM(cycles=1)
    ciao.measure_HOIM(config=True)
    ciao.get_InteractionMatrices()

    ciao.calc_CommandMatrix(nFiltModes=20)

    ciao.set_TT_gain(-0.1)
    ciao.measureNewTTRefPositions()

if short:
    angles = numpy.arange(10)*40.0
else:
    angles = numpy.arange(26)*14.4
#print "This is a test"

i = 0
for angle in angles:
    ciao.moveDerotator(angle)
    print "Moved to :",  angle
    if True:
        ciao.setup_TTIM(cycles=3)
        ciao.measure_TTIM(config=True)
        ciao.get_InteractionMatrices()
        ciao.calc_CommandMatrix(nFiltModes=20)
        ciao.measureNewTTRefPositions()

    time.sleep(3)
    ciao.setup_HOIM(cycles=3)
    ciao.measure_HOIM(config=True)
    ciao.saveMap(mapname="HORecnCalibrat.RESULT_IM", filename="IM_"+str(i)+"_.fits")
    ciao.disturb(disturbance=fname, baseName="Position_"+str(i))
    i += 1
