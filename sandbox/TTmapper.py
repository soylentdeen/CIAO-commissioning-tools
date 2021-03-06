import scipy
import numpy
import pyfits
import VLTTools
import SPARTATools
import os
import glob
import time
import errno

datadir = "/diska/data/SPARTA/2015-05-19/TTMap_3/"
ciao = VLTTools.VLTConnection(simulate=False, datapath=datadir)

#ciao.measureNewTTRefPositions("Alderan")
#ciao.measureNewHORefPositions("BetaPic")

i = 0

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


pupilShiftFile = open(os.path.expanduser('~')+'/data/TTMapping/TTMCommandSequence_2.txt', 'r')

Tip = []
Tilt = []

for line in pupilShiftFile:
    l = line.split()
    Tip.append(float(l[0]))
    Tilt.append(float(l[1]))



ciao.measureNewTTRefPositions("TWHydra")
TT = [ciao.get_Tip()[0], ciao.get_Tilt()[0]]

TTLoop = numpy.array([2.0, 5.0])
#TTLoop = numpy.array([5.0])
DrotLoop = numpy.arange(5)*72.0
#DrotLoop = numpy.array([0.0])

for drot in DrotLoop:
    ciao.moveDerotator(drot)
    for factor in TTLoop:
        #directory = "Derot_"+str(drot)+"/TTfactor_"+str(factor)
        #make_sure_path_exists(datadir+directory)
        for tip, tilt in zip(Tip, Tilt):
            ciao.set_Tip(factor*tip+TT[0])
            ciao.set_Tilt(factor*tilt+TT[1])
            print 'Moved to Tip, Tilt : ', factor*tip, factor*tilt

            time.sleep(3.5)
            ciao.measurePixelFrames("TTMapPixels_"+str(i)+"_"+str(factor*tip)+"_"+str(factor*tilt))
            i += 1

ciao.set_Tip(TT[0])
ciao.set_Tilt(TT[1])
