import scipy
import numpy
import pyfits
import VLTTools
import SPARTATools
import os
import glob
import time

datadir = "/diska/data/SPARTA/2015-04-09/FieldLens_2/"
ciao = VLTTools.VLTConnection(simulate=False, datapath=datadir)

i = 0

pupilShiftFile = open(os.path.expanduser('~')+'/data/FieldLens_test_data/FieldLensPath.txt', 'r')

xpos = []
ypos = []

for line in pupilShiftFile:
    l = line.split()
    xpos.append(float(l[0]))
    ypos.append(float(l[1]))


ciao.fieldLens.x = -316.8
ciao.fieldLens.y = -55.8
for j in ['0', '1', '2', '3']:
    for x, y in zip(xpos, ypos):
        ciao.moveFieldLens(x, y)
        print 'Moved to ', x, y

        time.sleep(0.5)
        ciao.setup_HOIM()
        ciao.measure_HOIM(config=True)
        ciao.saveMap(mapname="HORecnCalibrat.RESULT_IM", filename="FLT"+j+"_IM_"+str(x)+"_"+str(y)+"_"+str(i)+"_.fits")
        i += 1

