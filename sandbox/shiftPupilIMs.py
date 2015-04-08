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
datadir = "FieldLens_test_data/"

logfile = open(os.path.expanduser('~')+'/data/'+datadir+'logfile.dat', 'w')

pupilShiftFile = open(os.path.expanduser('~')+'/data/'+datadir+'FieldLensPath.txt', 'r')

#ciao.setupFieldLens()

ciao.fieldLens.x = -50
ciao.fieldLens.y = -50

for line in pupilShiftFile:
    l = line.split()
    x = float(l[0])
    y = float(l[1])

    ciao.moveFieldLens(x, y)
    print 'Moved to ', x, y

    time.sleep(0.5)
    ciao.setup_HOIM()
    ciao.measure_HOIM(config=True)
    ciao.saveMap(mapname="HORecnCalibrat.RESULT_IM", filename=datadir+"FLT2_IM_"+str(x)+"_"+str(y)+"_"+str(i)+"_.fits")
    i += 1


logfile.close()

