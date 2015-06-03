import os
import glob
import numpy
import pyfits
import warnings
import select
import logging
import subprocess
import time
import SPARTATools
import scipy
import matplotlib.pyplot as pyplot

class Derotator( object ):
    def __init__(self, parent):
        self.angle = 0.0
        self.encoder = 0
        self.parent = parent

    def initialize(self):
        command = ""
        self.parent.sendCommand(command)

    def moveToAngle(self, angle):
        command = "msgSend -n wci1ao ciiControl SETUP \" -function INS.DROT.ENC "+str(angle*1000)+" \""
        self.parent.sendCommand(command)
        self.angle = angle
        self.encoder = angle*1000 % 360000

    def moveAngleRel(self, deltaAngle):
        command = "msgSend -n wci1ao ciiControl SETUP \" -function INS.DROT.ENCREL "+str(deltaAngle*1000)+" \""
        self.parent.sendCommand(command)
        self.angle += deltaAngle
        self.angle %= 360.0
        self.encoder += deltaAngle*1000
        self.encoder %= 360000


class ParabolicMirror( object ):
    def __init__(self, parent):
        self.Tip = 0.0
        self.Tilt = 0.0
        self.parent = parent
    
    def initialize(self, Tip, Tilt):
        self.Tip = Tip
        self.Tilt = Tilt

    def moveToTip(self, Tip):
        command = "msgSend -n wci1ao ciiControl SETUP \" -function INS.PMTIP.ENC "+str(Tip)+" \""
        self.parent.sendCommand(command)
        self.Tip = Tip

    def moveToTilt(self, Tilt):
        command = "msgSend -n wci1ao ciiControl SETUP \" -function INS.PMTIL.ENC "+str(Tilt)+" \""
        self.parent.sendCommand(command)
        self.Tilt = Tilt

class FieldLens( object ):
    def __init__(self, parent):
        self.x = 0.0
        self.y = 0.0
        self.parent = parent
    
    def initialize(self, x, y):
        self.x = x
        self.y = y

    def moveToX(self, x):
        dx = x - self.x
        self.x = x
        command = "msgSend -n wci1ao ciiControl SETUP \" -function INS.FLDL.DX "+str(dx)+" INS.FLDL.DY 0 \""
        self.parent.sendCommand(command)
        """
        newY = self.y
        while newX != x:
            difference = x - newX
            sign = numpy.abs(difference)/difference
            stepsize = sign*min( 10.0, numpy.abs(difference))
            command = "msgSend -n wci1ao ciiControl SETUP \" -function INS.FLDL.DX "+str(stepsize)+" INS.FLDL.DY 0 \""
            #print command
            self.parent.sendCommand(command)
            newX += stepsize
        self.x = newX
        """

    def moveToY(self, y):
        dy = y - self.y
        self.y = y
        command = "msgSend -n wci1ao ciiControl SETUP \" -function INS.FLDL.DX 0 INS.FLDL.DY "+str(dy)+" \""
        self.parent.sendCommand(command)
        """
        newY = self.y
        while newY != y:
            difference = y - newY
            sign = numpy.abs(difference)/difference
            stepsize = sign*min(10.0, numpy.abs(difference))
            command = "msgSend -n wci1ao ciiControl SETUP \" -function INS.FLDL.DX 0 INS.FLDL.DY "+str(stepsize)+" \""
            #print command
            self.parent.sendCommand(command)
            newY += stepsize
        self.y = newY
        """

class VLTConnection( object ):
    """
    VLTConnection:  This object allows python to log into a computer
    running the VLT SPARTA Light software and do the following:
        - Send a new flat pattern to the DM
        - Retrieve data from the RTC (slopes, intensities, etc...)
        - what else?

    """
    def __init__(self, simulate=True, datapath=None):
        if datapath:
            self.datapath = datapath
        else:
            self.datapath = os.path.expanduser('~')+'/data/'
        self.CDMS = CDMS()
        self.sim = simulate
        self.modalBasis = None
        logging.basicConfig(level=logging.DEBUG)
        self.fieldLens = FieldLens(self)
        self.derotator = Derotator(self)
        self.PM = ParabolicMirror(self)

    def simulate(self):
        self.sim = True

    def goLive(self):
        self.sim = False

    def sendCommand(self, command, response=False):
        if not(self.sim):
            #logging.debug("Executing '%s'" % command)
            results = subprocess.check_output(command, shell=True)
            if response:
                return results
        else:
            print("In Simulation mode.  I would have executed the following command")
            print("'%s'" % command)
            if response:
                return "SIMULATION"

    def parse(self, text, type):
        lines = text.split('\n')
        retval = numpy.array(lines[1].split(','), dtype=type)
        return retval

    def applyPAF(self, paf):
        name = paf.name
        for key in paf.parameters.keys():
            par = paf.parameters[key]
            if isinstance(par, numpy.int):
                flag = '-i'
            elif isinstance(par, numpy.float):
                flag = '-d'
            else:
                flag = '-s'
            command = "cdmsSetProp "+name+" "+key+" "+flag+" "+str(par)
            self.sendCommand(command)


    def saveMap(self, mapname, filename):
        localfile = self.datapath+filename
        command = "cdmsSave -f "+localfile+" "+mapname
        self.sendCommand(command)

    def updateMap(self, mapname):
        localfile = self.datapath+self.CDMS.maps[mapname].outfile
        command = "cdmsSave -f "+localfile+" "+mapname
        self.sendCommand(command)
        self.CDMS.maps[mapname].load(localfile)

    def transmitMap(self, mapname, update=None):
        localfile = self.datapath+self.CDMS.maps[mapname].outfile
        self.CDMS.maps[mapname].write(path=self.datapath)
        command = "cdmsLoad -f "+localfile+" "+mapname+" --rename"
        self.sendCommand(command)
        if update:
            command = "msgSend \"\" spaccsServer EXEC \" -command "+update+".update ALL\""
            self.sendCommand(command)
        
    def zeroMap(self, mapName):
        self.CDMS.maps[mapName].scale(0.0)
        self.transmitMap(mapName)

    def mirrorCDMS(self):
        for mapname in self.CDMS.maps.keys():
            self.updateMap(mapname)

    def applyZernike(self, coeffs):
        self.updateMap("HOCtr.ACT_POS_REF_MAP")
        offsets = self.modalBasis.getZernikeOffsets(coeffs)
        self.CDMS.maps["HOCtr.ACT_POS_REF_MAP"].delta(offsets)
        self.transmitMap("HOCtr.ACT_POS_REF_MAP", update="HOCtr")

    def set_Tip(self, tip):
        self.CDMS.maps['TTCtr.ACT_POS_REF_MAP'].data[0][0] = tip
        self.transmitMap('TTCtr.ACT_POS_REF_MAP', update='TTCtr')
        #self.sendCommand("msgSend \"\" CDMSGateway SETMAP \"-object TTCtr.ACT_POS_REF_MAP -function 0,0= "+str("%.3g" % tip)+"\"")
        #self.sendCommand("msgSend \"\" spaccsServer EXEC \"-command TTCtr.update ALL\"")

    def set_Tilt(self, tilt):
        self.CDMS.maps['TTCtr.ACT_POS_REF_MAP'].data[0][1] = tilt
        self.transmitMap('TTCtr.ACT_POS_REF_MAP', update='TTCtr')
        #self.sendCommand("msgSend \"\" CDMSGateway SETMAP \"-object TTCtr.ACT_POS_REF_MAP -function 0,1= "+str("%.3g" % tilt)+"\"")
        #self.sendCommand("msgSend \"\" spaccsServer EXEC \"-command TTCtr.update ALL\"")

    def get_TipTilt(self):
        tip = self.sendCommand("msgSend \"\" CDMSGateway GETMAP \"-object TTCtr.ACT_POS_REF_MAP -function 0,0 1,0\"", response=True)
        return self.parse(tip, numpy.float32)

    def get_Tip(self):
        tip = self.sendCommand("msgSend \"\" CDMSGateway GETMAP \"-object TTCtr.ACT_POS_REF_MAP -function 0,0\"", response=True)
        return self.parse(tip, numpy.float32)

    def get_Tilt(self):
        tilt = self.sendCommand("msgSend \"\" CDMSGateway GETMAP \"-object TTCtr.ACT_POS_REF_MAP -function 1,0\"", response=True)
        return self.parse(tilt, numpy.float32)

    def set_TT_gain(self, gain):
        self.sendCommand("msgSend \"\" CDMSGateway SETMAP \"-object TTCtr.TERM_B -function 0,0="+str("%.2g" % gain)+"\"")
        self.sendCommand("msgSend \"\" spaccsServer EXEC \"-command TTCtr.update ALL\"")

    def set_HO_gain(self, gain):
        self.sendCommand("msgSend \"\" CDMSGateway SETMAP \"-object HOCtr.TERM_B -function 0,0="+str("%.2g" % gain)+"\"")
        self.sendCommand("msgSend \"\" spaccsServer EXEC \"-command HOCtr.update ALL\"")


    """
    def calc_CommandMatrix(self, nFiltModes=20):
        self.CDMS.maps['Recn.REC1.CM'].replace(
                  SPARTATools.calculateCommandMatrix(
                         self.CDMS.maps['HORecnCalibrat.RESULT_IM'],
                         self.CDMS.maps['TTRecnCalibrat.RESULT.IM'],
                         nFiltModes))
    #"""

    #"""

    def getTTFocus(self):
        TT = self.sendCommand("msgSend \"\" spaccsServer EXEC \" -command LoopMonitor.measureTipTilt 10 \"", response=True)
        TT = self.parse(TT, numpy.float32)

        focus = self.sendCommand("msgSend \"\" spaccsServer EXEC \" -command LoopMonitor.measureFocus 10 \"", response=True)
        focus = self.parse(focus, numpy.float32)

        return numpy.append(TT, focus)
        
    def getIntensities(self):
        self.sendCommand("msgSend \"\" spaccsServer EXEC \" -command AcqOptimiser.measureFrames 10 \"", response=False)


        return intensities

    def calc_CommandMatrix(self, nFiltModes=20):
        self.modalBasis = SPARTATools.modalBasis(self.CDMS.maps['HORecnCalibrat.RESULT_IM'].data, self.CDMS.maps['TTRecnCalibrat.RESULT.IM'].data, nFiltModes)
        self.modalBasis.computeSystemControlMatrix()
        self.CDMS.maps['Recn.REC1.CM'].replace(self.modalBasis.CM)
        self.CDMS.maps['HOCtr.TT_TO_HO'].replace(self.modalBasis.TT2HO)
        self.CDMS.maps['HOCtr.HO_TO_TT'].replace(self.modalBasis.HO2TT)
        self.CDMS.maps['HOCtr.SMA_BASIS'].replace(self.modalBasis.SMAbasis)
        #self.CDMS.maps['HOCtr.AWF_IM_KERNEL'].replace(self.modalBasis.AWFbasis)
        self.CDMS.maps['HOCtr.PRA_PISTON_MODE'].replace(self.modalBasis.pistonMode)
        self.CDMS.maps['HOCtr.PRA_PISTON_PROJECTION'].replace(self.modalBasis.pistonProj)
        self.CDMS.maps['LoopMonitor.SLOPES2FOCUS'].replace(self.modalBasis.S2Z[2,:])
        self.CDMS.maps['LoopMonitor.SLOPES2TT'].replace(self.modalBasis.S2Z[0:2,:])
        self.CDMS.maps['LoopMonitor.SLOPES2FOCUS'].write()
        self.CDMS.maps['LoopMonitor.SLOPES2TT'].write()
        self.transmitMap('HOCtr.TT_TO_HO', update='HOCtr')
        self.transmitMap('HOCtr.HO_TO_TT', update='HOCtr')
        self.transmitMap('Recn.REC1.CM', update='Recn')
        self.transmitMap('HOCtr.SMA_BASIS', update='HOCtr')
        #self.transmitMap('HOCtr.AWF_IM_KERNEL', update='HOCtr')
        self.transmitMap('HOCtr.PRA_PISTON_MODE', update='HOCtr')
        self.transmitMap('HOCtr.PRA_PISTON_PROJECTION', update='HOCtr')
        self.transmitMap('LoopMonitor.SLOPES2FOCUS', update='LoopMonitor')
        self.transmitMap('LoopMonitor.SLOPES2TT', update='LoopMonitor')
        self.CDMS.maps['Recn.REC1.CM'].replace(self.modalBasis.CM)
    #"""
        
    def dumpCommandMatrix(self, nFiltModes=20):
        self.modalBasis = SPARTATools.modalBasis(self.CDMS.maps['HORecnCalibrat.RESULT_IM'].data, self.CDMS.maps['TTRecnCalibrat.RESULT.IM'].data, nFiltModes)
        self.modalBasis.computeSystemControlMatrix()
        hdu = pyfits.PrimaryHDU(numpy.array(self.modalBasis.M2V, dtype=numpy.float32))
        hdu.writeto("M2V.fits", clobber=True)
        hdu = pyfits.PrimaryHDU(numpy.array(self.modalBasis.TT2HO, dtype=numpy.float32))
        hdu.writeto("TT2HO.fits", clobber=True)
        hdu = pyfits.PrimaryHDU(numpy.array(self.modalBasis.S2M, dtype=numpy.float32))
        hdu.writeto("S2M.fits", clobber=True)
        hdu = pyfits.PrimaryHDU(numpy.array(self.modalBasis.S2Z, dtype=numpy.float32))
        hdu.writeto("S2Z.fits", clobber=True)
        zeroCM = numpy.zeros(self.modalBasis.CM.shape, dtype=numpy.float32)
        hdu = pyfits.PrimaryHDU(self.modalBasis.CM)
        hdu.writeto("CM.fits", clobber=True)
        zeroCM = numpy.zeros(62, dtype=numpy.int)
        hdu = pyfits.PrimaryHDU(zeroCM)
        hdu.writeto("RTC_HODM_UNUSED_ACT_MAP.fits", clobber=True)
        #zeroCM = numpy.zeros(self.modalBasis.CM.shape, dtype=numpy.float32)
        #hdu = pyfits.PrimaryHDU(zeroCM)
        hdu.writeto("HODM_SLUG_ACT_GAINS.fits", clobber=True)

        CM = numpy.dot(self.modalBasis.M2V, self.modalBasis.S2M)
        tiltHO = numpy.dot(self.modalBasis.TT2HO, CM[60:62,:])
        aux_proj = numpy.zeros((62, 136), dtype=numpy.float32)
        aux_proj[0:60,:] += tiltHO
        hdu = pyfits.PrimaryHDU(aux_proj)
        hdu.writeto("CLMatrixOptimiser.AUX_PROJ.fits", clobber=True)
        self.CDMS.maps["RTC.MODAL_GAINS"].write("RTC.MODAL_GAINS.fits")


    def updateRefSlopes(self, x, y):
        slopes = numpy.zeros(136)
        slopes[0::2] += x
        slopes[1::2] += y
        self.CDMS.maps["Acq.DET1.REFSLP"].replace(slopes)
        self.transmitMap("Acq.DET1.REFSLP", update='Acq')

    def replaceSlopesWithCurrent(self, rec5rdingName="BetaPic"):
        self.measureCircularBuffer(recordingName=recordingName)
        outfile="gradients.fits"
        time.sleep(2.0)
        SPARTATools.computeGradients(outfile, recordingName)
        self.updateReferenceSlopes(outfile)

    def updateReferenceSlopes(self, filename):
        self.CDMS.maps["Acq.DET1.REFSLP"].load(filename)
        self.transmitMap("Acq.DET1.REFSLP", update='Acq')

    def measureCircularBuffer(self, recordingName="Arcturus", nframes=100):
        command = "msgSend \"\" spaccsServer SETUP \"-function LoopRecorder.FILE_BASENAME "+recordingName+"\""
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer SETUP \"-function LoopRecorder.FILE_DIRNAME "+self.datapath+"\""
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer SETUP \"-function LoopRecorder.REQUESTED_FRAMES "+str(nframes)+"\""
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer EXEC \" -command LoopRecorder.run\""
        self.sendCommand(command)
        print "Took Circular Buffer"

    def measurePixelFrames(self, recordingName="Arcturus", nframes=100):
        command = "msgSend \"\" spaccsServer SETUP \"-function PixelRecorder.FILE_BASENAME "+recordingName+"\""
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer SETUP \"-function PixelRecorder.FILE_DIRNAME "+self.datapath+"\""
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer SETUP \"-function PixelRecorder.REQUESTED_FRAMES "+str(nframes)+"\""
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer EXEC \" -command PixelRecorder.run\""
        self.sendCommand(command)
        print "Took Pixel Frame"

    def measureNewHORefPositions(self, recordingName='TWHydra', nframes=25):
        command = "msgSend \"\" spaccsServer EXEC \" -command TTCtr.openLoop\""
        self.sendCommand(command)
        time.sleep(2.0)
        command = "msgSend \"\" spaccsServer EXEC \" -command HOCtr.closeLoop\""
        self.sendCommand(command)
        time.sleep(5.0)

        """
        command = "msgSend \"\" CommandGateway EXEC \" LoopMonitor.measureRefDMPos "+str(nframes)+"\""
        self.sendCommand(command)
        #"""
        print "Recording Frames"
        self.measureCircularBuffer(recordingName=recordingName)
        time.sleep(2.0)
        #"""
        print "Opening Loop"
        command = "msgSend \"\" spaccsServer EXEC \" -command HOCtr.openLoop\""
        self.sendCommand(command)
        time.sleep(2.0)
        self.averageHOPositions(recordingName)

    def measureNewTTRefPositions(self, recordingName='TWHydra', nframes=25):
        command = "msgSend \"\" spaccsServer EXEC \" -command HOCtr.openLoop\""
        self.sendCommand(command)
        time.sleep(2.0)
        command = "msgSend \"\" spaccsServer EXEC \" -command TTCtr.closeLoop\""
        self.sendCommand(command)
        time.sleep(5.0)
        
        """
        command = "msgSend \"\" CommandGateway EXEC \" LoopMonitor.measureRefTTMPos "+str(nframes)+"\""
        self.sendCommand(command)
        #"""
        print "Recording Frames"
        self.measureCircularBuffer(recordingName=recordingName)
        time.sleep(2.0)
        #"""
        print "Opening Loop"
        command = "msgSend \"\" spaccsServer EXEC \" -command TTCtr.openLoop\""
        self.sendCommand(command)
        time.sleep(2.0)
        self.averageTTPositions(recordingName)

    def averageActuatorPositions(self, recordingName):
        outfile= self.datapath+"new_flat.fits"
        SPARTATools.computeNewBestFlat(outfile, self.datapath, recordingName)
        command = "cdmsLoad -f "+outfile+" HOCtr.ACT_POS_REF_MAP --rename"
        self.sendCommand(command)

    def averageTTPositions(self, recordingName):
        outfile=self.datapath+"new_TT_flat.fits"
        SPARTATools.computeNewTTFlat(outfile, self.datapath, recordingName)
        command = "cdmsLoad -f "+outfile+" TTCtr.ACT_POS_REF_MAP --rename"
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer EXEC \" -command TTCtr.update ALL\""
        self.sendCommand(command)

    def averageHOPositions(self, recordingName):
        outfile=self.datapath+"new_HO_flat.fits"
        SPARTATools.computeNewHOFlat(outfile, self.datapath, recordingName)
        command = "cdmsLoad -f "+outfile+" HOCtr.ACT_POS_REF_MAP --rename"
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer EXEC \" -command HOCtr.update ALL\""
        self.sendCommand(command)
    
    def averageIntensities(self):
        outfile=self.datapath+"averageIntensities.fits"
        SPARTATools.computeIntensities(outfile)

    def set_CommandMatrix(self):
        #self.transmitMap('Recn.REC1.CM')
        self.transmitMap('Recn.REC1.CM', update='Recn')

    def save_CommandMatrixPlot(self):
        self.updateMap("Recn.REC1.CM")
        fig = pyplot.figure(0)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.imshow(self.CDMS.maps['Recn.REC1.CM'].data)
        fig.savefig("CM.png")

    def disturb(self, disturbance="disturbance.fits", baseName = "TWHydra"):
        command = "msgSend \"\" spaccsServer SETUP \"-function LoopRecorder.FILE_DIRNAME "+self.datapath+"\""
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer SETUP \"-function PixelRecorder.FILE_DIRNAME "+self.datapath+"\""
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer SETUP \"-function PixelRecorder.FILE_BASENAME "+baseName+"_pixels\""
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer SETUP \"-function LoopRecorder.FILE_BASENAME "+baseName+"_frames\""
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer SETUP \"-function LoopRecorder.REQUESTED_FRAMES 0\""
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer SETUP \"-function PixelRecorder.REQUESTED_FRAMES 0\""
        self.sendCommand(command)

        self.CDMS.paf["HOCtrDisturb.CFG.DYNAMIC"].update("FILENAME", self.datapath+"DM"+disturbance)
        self.CDMS.paf["TTCtrDisturb.CFG.DYNAMIC"].update("FILENAME", self.datapath+"TTM"+disturbance)
        startTime = numpy.int(time.time()+20)
        self.CDMS.paf["HOCtrDisturb.CFG.DYNAMIC"].update("START_AT", startTime)
        self.CDMS.paf["TTCtrDisturb.CFG.DYNAMIC"].update("START_AT", startTime)
        self.applyPAF(self.CDMS.paf["HOCtrDisturb.CFG.DYNAMIC"])
        self.applyPAF(self.CDMS.paf["TTCtrDisturb.CFG.DYNAMIC"])

        command="msgSend \"\" spaccsServer EXEC \"-command LoopRecorder.run\""
        self.sendCommand(command)
        command="msgSend \"\" spaccsServer EXEC \"-command PixelRecorder.run\""
        self.sendCommand(command)
        
        command="msgSend \"\" spaccsServer EXEC \" -command HOCtrDisturb.run \""
        #command = "spaciaortdfwDisturbPubl -d HODM -f " + self.datapath+disturbance
        self.sendCommand(command)
        command="msgSend \"\" spaccsServer EXEC \" -command TTCtrDisturb.run\""
        self.sendCommand(command)

        command = "dbRead \"<alias>SPARTA:HOCtrDisturb.percent_complete\""
        complete = 0
        time.sleep(5.0)
        while complete < 100:
            wait = self.sendCommand(command, response=True)
            complete = numpy.float(wait.split()[-1])
            print complete
            time.sleep(5.0)

        command="msgSend \"\" spaccsServer EXEC \"-command LoopRecorder.idle\""
        self.sendCommand(command)
        command="msgSend \"\" spaccsServer EXEC \"-command PixelRecorder.idle\""
        self.sendCommand(command)
    
    def disturbHO(self, disturbType="SINE", rng=0.5, max=0.95, period=40.0, actNum=5):
        fname = self.datapath+"Disturbances/disturbanceFrame.fits"
        SPARTATools.computeHODisturbanceFrame(20000,fname, rng=rng, max=max, disturbType=disturbType, period=period, actNum=actNum)
        command = "spaciaortdfwDisturbPubl -d HODM -f "+fname
        self.sendCommand(command)

    def disturbTT(self, tip=0.05, tilt=0.05, waveshape='SQUARE'):
        fname = self.datapath+"disturbanceTTFrame.fits"
        SPARTATools.computeTTDisturbanceFrame(1000, fname, tip, tilt, waveshape=waveshape)
        command = "spaciaortdfwDisturbPubl -d ITTM -f "+fname+" -m 1000"
        self.sendCommand(command)

    def measure_HOIM(self, config=None):
        if config:
            self.applyPAF(self.CDMS.paf["HORecnCalibrat.CFG.DYNAMIC"])

        #command = "msgSend \"\" spaccsServer EXEC \" -command HOCtrUpload.run\""
        #self.sendCommand(command)
        command = "msgSend \"\" spaccsServer EXEC \" -command HORecnCalibrat.update ALL\""
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer EXEC \" -command HORecnCalibrat.run\""
        self.sendCommand(command)
        #command = "msgSend \"\" spaccsServer EXEC \" -command HORecnCalibrat.waitIdle\""
        command = "dbRead \"<alias>SPARTA:HORecnCalibrat.percent_complete\""
        complete = 0
        time.sleep(5.0)
        while complete < 100:
            wait = self.sendCommand(command, response=True)
            complete = numpy.float(wait.split()[-1])
            #print complete
            time.sleep(5.0)

    def moveFieldLens(self, x, y):
        #Move in X
        self.fieldLens.moveToX(x)
        #Move in Y
        self.fieldLens.moveToY(y)

    def movePM(self, Tip, Tilt):
        #Move in X
        self.PM.moveToTip(Tip)
        #Move in Y
        self.PM.moveToTilt(Tilt)

    def moveDerotator(self, angle):
        #Move to angle
        self.derotator.moveToAngle(angle)

    def measure_TTIM(self, config=None):
        if config:
            self.applyPAF(self.CDMS.paf["TTRecnCalibrat.CFG.DYNAMIC"])

        #command = "msgSend \"\" spaccsServer EXEC \" -command HOCtrUpload.run\""
        #self.sendCommand(command)
        command = "msgSend \"\" spaccsServer EXEC \" -command TTRecnCalibrat.update ALL\""
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer EXEC \" -command TTRecnCalibrat.run\""
        self.sendCommand(command)
        #command = "msgSend \"\" spaccsServer EXEC \" -command HORecnCalibrat.waitIdle\""
        command = "dbRead \"<alias>SPARTA:TTRecnCalibrat.percent_complete\""
        complete = 0
        time.sleep(5.0)
        while complete < 100:
            wait = self.sendCommand(command, response=True)
            complete = numpy.float(wait.split()[-1])
            #print complete
            time.sleep(5.0)

    def setup_HOIM(self, amplitude=1.0, noise=0.05, skip=0.05,
                   period=0.2, mode_cycles=1, cycles=3):
        self.CDMS.paf["HORecnCalibrat.CFG.DYNAMIC"].update("ACTUATION_MATRIX", "HORecnCalibrat.USER_60")
        self.CDMS.paf["HORecnCalibrat.CFG.DYNAMIC"].update("ACTUATION_MATRIX_INV", "HORecnCalibrat.USER_INV_60")
        self.CDMS.paf["HORecnCalibrat.CFG.DYNAMIC"].update("TIME_UNIT",
                  "SECONDS")
        self.CDMS.paf["HORecnCalibrat.CFG.DYNAMIC"].update("WAVE_PERIOD",
                  period)
        self.CDMS.paf["HORecnCalibrat.CFG.DYNAMIC"].update("AMPLITUDE",
                  amplitude)
        self.CDMS.paf["HORecnCalibrat.CFG.DYNAMIC"].update("NOISE_THRESHOLD",
                  noise)
        self.CDMS.paf["HORecnCalibrat.CFG.DYNAMIC"].update("SKIP_TIME", skip)
        self.CDMS.paf["HORecnCalibrat.CFG.DYNAMIC"].update("CYCLES", cycles)
        self.CDMS.paf["HORecnCalibrat.CFG.DYNAMIC"].update("MODE_CYCLES", mode_cycles)
    
    def setup_TTIM(self, amplitude=0.05, skip=0.2,
                   period=2.0, mode_cycles=1, cycles=3):
        self.CDMS.paf["TTRecnCalibrat.CFG.DYNAMIC"].update("TIME_UNIT",
                  "SECONDS")
        self.CDMS.paf["TTRecnCalibrat.CFG.DYNAMIC"].update("WAVE_PERIOD",
                  period)
        self.CDMS.paf["TTRecnCalibrat.CFG.DYNAMIC"].update("AMPLITUDE",
                  amplitude)
        self.CDMS.paf["TTRecnCalibrat.CFG.DYNAMIC"].update("SKIP_TIME", skip)
        self.CDMS.paf["TTRecnCalibrat.CFG.DYNAMIC"].update("CYCLES", cycles)
        self.CDMS.paf["TTRecnCalibrat.CFG.DYNAMIC"].update("MODE_CYCLES", mode_cycles)
        self.CDMS.paf["TTRecnCalibrat.CFG.DYNAMIC"].update("INTERACTION_MATRIX", 'TTRecnCalibrat.RESULT.IM')

    def get_HOIM(self):
        self.updateMap('HORecnCalibrat.RESULT_IM')

    def get_InteractionMatrices(self):
        self.updateMap('HORecnCalibrat.RESULT_IM')
        self.updateMap('TTRecnCalibrat.RESULT.IM')

    def set_InteractionMatrices(self, HOIM_file, TTIM_file):
        self.CDMS.maps["HORecnCalibrat.RESULT_IM"].load(HOIM_file)
        self.CDMS.maps["TTRecnCalibrat.RESULT.IM"].load(TTIM_file)
        #self.transmitMap("HORecnCalibrat.RESULT_IM")
        #self.transmitMap("TTRecnCalibrat.RESULT.IM")


    def changePixelTapPoint(self, tp):
        try:
            if tp == "RAW":
                pass
            elif tp == "CALIB":
                pass
            elif tp == "BACKGROUND":
                pass
            else:
                print("Error!  Unrecognized tap point!")
                escape
            command="cdmsSetProp Acq.CFG.DYNAMIC DET1.PIXEL_TAP -s \""+tp+"\""
            self.sendCommand(command)
            command = "msgSend \"\" spaccsServer EXEC \" -command Acq.update ALL\""
            self.sendCommand(command)
        except:
            print("Error!  Invalid tap point!")

    def measureBackground(self, nframes):
        self.changePixelTapPoint("RAW")
        command="msgSend \"\" CommandGateway EXEC \"AcqOptimiser.measureBackground "+str(nframes)+"\""
        self.sendCommand(command)
        self.updateMap('Acq.DET1.BACKGROUND')
        self.changePixelTapPoint("CALIB")

class CDMS_Map( object ):
    def __init__(self, name, ax1, ax2, dtype, filltype, bscale):
        if dtype == "float32":
            self.dtype = numpy.float32
        elif dtype == "float16":
            self.dtype = numpy.float16
        elif dtype == "int32":
            self.dtype = numpy.int32
        elif dtype == "int16":
            self.dtype = numpy.int16
        else:
            print "Error!"
        if filltype == 0.0:
            self.data = numpy.zeros((ax1, ax2), dtype=self.dtype)
        elif filltype >= 1.0:
            self.data = numpy.ones((ax1, ax2), dtype=self.dtype)*filltype
        elif filltype == -1.0:
            self.data = numpy.arange(ax1, dtype=self.dtype)
        else:
            print "Error! I can't understand the fill type!"
        self.data_template = self.data.copy()
        self.bscale = bscale
        self.outfile = name+'.fits'

    def replace(self, newmap):
        self.data = self.dtype(newmap).copy()

    def load(self, file):
        self.data = pyfits.getdata(file)
        
    def revert(self):
        self.data = self.data_template.copy()

    def delta(self, offsets):
        self.data += offsets

    def scale(self, factor):
        self.data *= factor

    def write(self, path=''):
        self.hdu = pyfits.PrimaryHDU(self.data)
        if self.bscale == 'minmax':
            self.hdu.scale(option='minmax')
        elif self.bscale == 'True':
            self.hdu.scale()
        warnings.resetwarnings()
        warnings.filterwarnings('ignore', category=UserWarning, append=True)
        self.hdu.writeto(path+self.outfile, clobber=True)
        warnings.resetwarnings()
        warnings.filterwarnings('always', category=UserWarning, append=True)


class PAF_File( object ):
    def __init__(self, filename, name):
        self.name = name
        self.file = filename
        self.parameters = {}
        file = open(filename, 'r')
        for line in file:
            l = line.split()
            if len(l) > 0:
                if (l[0].find(name) == 0) & (l[0][0] != '#'):
                    parameter = l[0][len(name)+1:]
                    if (l[1][:-1].find('.') != -1):
                        try:
                            val = numpy.float(l[1][:-1])
                        except:
                            val = l[1][:-1]
                    else:
                        try:
                            val = numpy.int(l[1][:-1])
                        except:
                            val = l[1][:-1]
                    self.parameters[parameter] = val
    
    def update(self, parameter, value):
        self.parameters[parameter] = value

class CDMS( object ):
    def __init__(self):
        self.maps = {}
        self.populateMapDefs()
        self.paf = {}
        self.populatePAF()

    def populateMapDefs(self):
        definitionFile = os.path.dirname(__file__)+'/CDMS_Map_Definitions.dat'
        df = open(definitionFile, 'r')
        for line in df:
            l = line.split(',')
            name = l[0]
            ax1 = int(l[1])
            ax2 = int(l[2])
            dtype = l[3].strip()
            filltype = float(l[4])
            bscale = bool(l[5])
            self.maps[name] = CDMS_Map(name, ax1, ax2, dtype, filltype, bscale)

    def populatePAF(self):
        pafdirectory = os.path.dirname(__file__)+'/PAF/'
        paffiles = glob.glob(pafdirectory+'*.paf')
        for paf in paffiles:
            name = paf[len(pafdirectory):-4]
            self.paf[name] = PAF_File(paf, name)
