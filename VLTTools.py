import os
import numpy
import pyfits
import warnings
import select
import logging
import subprocess
import SPARTATools

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
        logging.basicConfig(level=logging.DEBUG)

    def simulate(self):
        self.sim = True

    def goLive(self):
        self.sim = False

    def sendCommand(self, command, response=False):
        if not(self.sim):
            logging.debug("Executing '%s'" % command)
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


    def calc_CommandMatrix(self, nFiltModes=20):
        self.CDMS.maps['Recn.REC1.CM'].replace(
                  SPARTATools.calculateCommandMatrix(
                         self.CDMS.maps['HORecnCalibrat.RESULT_IM'],
                         self.CDMS.maps['TTRecnCalibrat.RESULT.IM'],
                         nFiltModes))
        
    def set_CommandMatrix(self):
        print("This is where the Command matrix is uploaded to SL and updated")
        #self.transmitMap('Recn.REC1.CM')
        self.transmitMap('Recn.REC1.CM', update='Recn')

    def measure_InteractionMatrices(self):
        print("This is where the Interaction Matrices will be measured")

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

    def updateAcq(self):
        stdin, stdout, stderr = self.ssh.exec_command("msgSend \"\" spaccsServer EXEC \"-command Acq.update ALL\"")
        while not stdout.channel.exit_status_ready():
            if stdout.channel.recv_ready():
                rl, wl, xl = select.select([stdout.channel], [], [], 0.0)
                if len(rl) > 0:
                    print stdout.channel.recv(1024)
    

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


class CDMS( object ):
    def __init__(self):
        self.maps = {}
        self.populateMapDefs()

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

