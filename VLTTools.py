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

class VLTConnection( object ):
    """
    VLTConnection:  This object allows python to log into a computer
    running the VLT SPARTA Light software and do the following:
        - Send a new flat pattern to the DM
        - Retrieve data from the RTC (slopes, intensities, etc...)
        - what else?

    """
    def __init__(self, simulate=True):
        self.datapath = os.path.expanduser('~')+'/data/'
        self.CDMS = CDMS()
        self.sim = simulate
        self.modalBasis = None
        logging.basicConfig(level=logging.DEBUG)

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
        retval = numpy.array(lines[1].split(), dtype=type)
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
        
    def mirrorCDMS(self):
        for mapname in self.CDMS.maps.keys():
            self.updateMap(mapname)

    def set_Tip(self, tip):
        self.CDMS.maps['TTCtr.ACT_POS_REF_MAP'].data[0][0] = tip
        self.sendCommand("msgSend \"\" CDMSGateway SETMAP \"-object TTCtr.ACT_POS_REF_MAP -function 0,0="+str("%.2g" % tip)+"\"")
        self.sendCommand("msgSend \"\" spaccsServer EXEC \"-command TTCtr.update ALL\"")

    def set_Tilt(self, tilt):
        self.CDMS.maps['TTCtr.ACT_POS_REF_MAP'].data[0][1] = tilt
        self.sendCommand("msgSend \"\" CDMSGateway SETMAP \"-object TTCtr.ACT_POS_REF_MAP -function 0,1="+str("%.2g" % tilt)+"\"")
        self.sendCommand("msgSend \"\" spaccsServer EXEC \"-command TTCtr.update ALL\"")

    def get_TipTilt(self):
        tip = self.sendCommand("msgSend \"\" CDMSGateway GETMAP \"-object TTCtr.ACT_POS_REF_MAP -function 0,0 1,0\"", response=True)
        return self.parse(tip, numpy.float32)

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
        self.transmitMap('HOCtr.TT_TO_HO', update='HOCtr')
        self.transmitMap('HOCtr.HO_TO_TT', update='HOCtr')
        self.transmitMap('Recn.REC1.CM', update='Recn')
        self.transmitMap('HOCtr.SMA_BASIS', update='HOCtr')
        #self.transmitMap('HOCtr.AWF_IM_KERNEL', update='HOCtr')
        self.transmitMap('HOCtr.PRA_PISTON_MODE', update='HOCtr')
        self.transmitMap('HOCtr.PRA_PISTON_PROJECTION', update='HOCtr')
        #self.CDMS.maps['Recn.REC1.CM'].replace(self.modalBasis.CM)
    #"""
        
    def updateReferenceSlopes(self, filename):
        self.CDMS.maps["Acq.DET1.REFSLP"].load(filename)
        self.transmitMap("Acq.DET1.REFSLP", update='Acq')

    def measureNewTTRefPositions(self, recordingName):
        command = "msgSend \"\" spaccsServer EXEC \" -command TTCtr.closeLoop\""
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer EXEC \" -command LoopRecorder.run\""
        self.sendCommand(command)
        time.sleep(5.0)
        command = "msgSend \"\" spaccsServer EXEC \" -command TTCtr.openLoop\""
        self.sendCommand(command)
        self.averageTTPositions(recordingName)

    def averageActuatorPositions(self, recordingName):
        outfile= self.datapath+"new_flat_16.fits"
        SPARTATools.computeNewBestFlat(outfile. recordingName)
        command = "cdmsLoad -f "+outfile+" HOCtr.ACT_POS_REF_MAP --rename"
        self.sendCommand(command)

    def averageTTPositions(self, recordingName):
        outfile=self.datapath+"new_TT_flat.fits"
        SPARTATools.computeNewTTFlat(outfile, recordingName)
        command = "cdmsLoad -f "+outfile+" TTCtr.ACT_POS_REF_MAP --rename"
        self.sendCommand(command)
        command = "msgSend \"\" spaccsServer EXEC \" -command TTCtr.update ALL\""
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
    
    def disturb(self):
        fname = self.datapath+"Disturbances/disturbanceFrame.fits"
        SPARTATools.computeDisturbanceFrame(5, 2000, fname, range=0.2, max=0.3)
        command = "spaciaortdftDisturbPubl -d HODM -f "+fname
        self.sendCommand(command)

    def measure_HOIM(self, config=None):
        if config:
            self.applyPAF(self.CDMS.paf["HORecnCalibrat.CFG.DYNAMIC"])

        command = "msgSend \"\" spaccsServer EXEC \" -command HOCtrUpload.run\""
        self.sendCommand(command)
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
            print complete
            time.sleep(5.0)

    def measure_TTIM(self, config=None):
        if config:
            self.applyPAF(self.CDMS.paf["TTRecnCalibrat.CFG.DYNAMIC"])
        command = "msgSend \"\" spaccsServer EXEC \" -command TTRecnCalibrat.run\""
        self.sendCommand(command)

    def setup_HOIM(self, amplitude=1.0, noise=0.05, skip=0.01,
                   period=0.2, cycles=10):
        self.CDMS.paf{"HORecnCalibrat.CFG.DYNAMIC"].update("ACTUATION_MATRIX", "HORecnCalibrat.USER_60")
        self.CDMS.paf{"HORecnCalibrat.CFG.DYNAMIC"].update("ACTUATION_MATRIX_INV", "HORecnCalibrat.USER_INV_60")
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
    
    def get_HOIM(self):
        self.updateMap('HORecnCalibrat.RESULT_IM')

    def get_InteractionMatrices(self):
        self.updateMap('HORecnCalibrat.RESULT_IM')
        self.updateMap('TTRecnCalibrat.RESULT.IM')

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
        except:
            print("Error!  Invalid tap point!")

    def measureBackground(self, nframes):
        command="msgSend \"\" CommandGateway EXEC \"AcqOptimiser.measureBackground "+str(nframes)+"\""
        self.sendCommand(command)
        self.updateMap('Acq.DET1.BACKGROUND')

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
