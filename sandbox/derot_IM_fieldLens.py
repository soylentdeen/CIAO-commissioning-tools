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
datadir = "/diska/data/SPARTA/2015-04-09/Derotator_20/"
ciao = VLTTools.VLTConnection(simulate=False, datapath=datadir)

offsetFile = '/diska/data/SPARTA/2015-04-09/Derotator_19/FieldLensModel.txt'
offsetX = []
offsetY = []

for line in open(offsetFile, 'r'):
    l = line.split()
    offsetX.append(float(l[0]))
    offsetY.append(float(l[1]))

ciao.fieldLens.initialize(0.0, 0.0)

#x = -260.0
#y = +30.0

i = 0
for angle, dx, dy in zip(numpy.arange(10)*40.0, offsetX, offsetY):

    ciao.moveDerotator(angle)
    print 'Moved Derotator to ', angle
    ciao.moveFieldLens(dx, dy)
    print "Moved Field Lens to :", dx, dy

    #"""
    ciao.measureCircularBuffer("derot_circbuff_"+str(i)+"_")
    ciao.setup_HOIM(cycles=1)
    ciao.measure_HOIM(config=True)
    ciao.saveMap(mapname="HORecnCalibrat.RESULT_IM", filename="IM_"+str(i)+"_.fits")
    i += 1
    #"""

ciao.moveFieldLens(0.0, 0.0)

