import os
import numpy
import pyfits
import warnings
import select
import logging
import subprocess

class VLTConnection( object ):
    """
    VLTConnection:  This object allows python to log into a computer
    running the VLT SPARTA Light software and do the following:
        - Send a new flat pattern to the DM
        - Retrieve data from the RTC (slopes, intensities, etc...)
        - what else?

    """
    def __init__(self, simulate=True):
        self.localpath = './data/'
        self.remotepath = './local/test/'
        self.CDMS = CDMS()
        self.sim = simulate
        self.logging = logging.basicConfig(level=logging.DEBUG)

    def simulate(self):
        self.sim = True

    def goLive(self):
        self.sim = False

    def sendCommand(self, command, response=False):
        if not(self.sim):
            self.log.debug("Executing '%s'" % command)
            exitCode = os.system(command)

    def parse(self, text):
        print text
        return 0.0

    def set_Tip(self, tip):
        self.sendCommand("msgSend \"\" CDMSGateway SETMAP \"-object TTCtr.ACT_POS_REF_MAP -function 0,0="+str("%.2g" % tip)+"\"")
        self.sendCommand("msgSend \"\" spaccsServer EXEC \"-command TTCtr.update ALL\"")

    def set_Tilt(self, tilt):
        self.sendCommand("msgSend \"\" CDMSGateway SETMAP \"-object TTCtr.ACT_POS_REF_MAP -function 0,1="+str("%.2g" % tilt)+"\"")
        self.sendCommand("msgSend \"\" spaccsServer EXEC \"-command TTCtr.update ALL\"")

    def get_Tip(self):
        tip = self.sendCommand("msgSend \"\" CDMSGateway GETMAP \"-object TTCtr.ACT_POS_REF_MAP -function 0,0\"", response=True)
        return self.parse(tip)

    def get_Tilt(self):
        tilt = self.sendCommand("msgSend \"\" CDMSGateway GETMAP \"-object TTCtr.ACT_POS_REF_MAP -function 0,1\"", response=True)
        return self.parse(tilt)

    def set_TT_gain(self, gain):
        self.sendCommand("msgSend \"\" CDMSGateway SETMAP \"-object TTCtr.TERM_B -function 0,0="+str("%.2g" % gain)+"\"")
        self.sendCommand("msgSend \"\" spaccsServer EXEC \"-command TTCtr.update ALL\"")

    def set_HO_gain(self, gain):
        self.sendCommand("msgSend \"\" CDMSGateway SETMAP \"-object HOCtr.TERM_B -function 0,0="+str("%.2g" % gain)+"\"")
        self.sendCommand("msgSend \"\" spaccsServer EXEC \"-command HOCtr.update ALL\"")

    def get_HO_ACT_POS_REF_MAP(self):
        name = self.CDMS.maps["HOCtr.ACT_POS_REF_MAP"].outfile
        command = "cdmsSave -f "+self.remotepath+name+" HOCtr.ACT_POS_REF_MAP"
        self.sendCommand(command)
        if not(self.sim):
            self.ftp.get(self.remotepath+name, self.localpath+name)
        return pyfits.getdata(self.localpath+name)

    def get_TT_ACT_POS_REF_MAP(self):
        name = self.CDMS.maps["TTCtr.ACT_POS_REF_MAP"].outfile
        command = "cdmsSave -f "+self.remotepath+name+" TTCtr.ACT_POS_REF_MAP"
        self.sendCommand(command)
        if not(self.sim):
            self.ftp.get(self.remotepath+name, self.localpath+name)
        return pyfits.getdata(self.localpath+name)

    def set_gain(self, gain):
        termA = numpy.array([-1], dtype='float32')
        termB = gain*(numpy.array([-1.0, 0.0], dtype='float32'))
        self.CDMS.maps["HOCtr.TERM_A"].replace(termA)
        self.CDMS.maps["HOCtr.TERM_B"].replace(termB)
        self.CDMS.maps["HOCtr.TERM_A"].write(path=self.localpath)
        self.CDMS.maps["HOCtr.TERM_B"].write(path=self.localpath)
        nameA = self.CDMS.maps["HOCtr.TERM_A"].outfile
        nameB = self.CDMS.maps["HOCtr.TERM_B"].outfile
        self.ftp.put(self.localpath+nameA, self.remotepath+nameA)
        self.ftp.put(self.localpath+nameB, self.remotepath+nameB)
        stdin, stdout, stderr = self.ssh.exec_command("cdmsLoad -f "+self.remotepath+nameA+" HOCtr.TERM_A --rename")
        while not stdout.channel.exit_status_ready():
            if stdout.channel.recv_ready():
                rl, wl, xl = select.select([stdout.channel], [], [], 0.0)
                if len(rl) > 0:
                    print stdout.channel.recv(1024)
        stdin, stdout, stderr = self.ssh.exec_command("cdmsLoad -f "+self.remotepath+nameB+" HOCtr.TERM_B --rename")
        while not stdout.channel.exit_status_ready():
            if stdout.channel.recv_ready():
                rl, wl, xl = select.select([stdout.channel], [], [], 0.0)
                if len(rl) > 0:
                    print stdout.channel.recv(1024)
        termB = gain*(numpy.array([-1.0, 0.0], dtype='float32'))
        self.CDMS.maps["TTCtr.TERM_A"].replace(termA)
        self.CDMS.maps["TTCtr.TERM_B"].replace(termB)
        self.CDMS.maps["TTCtr.TERM_A"].write(path=self.localpath)
        self.CDMS.maps["TTCtr.TERM_B"].write(path=self.localpath)
        nameA = self.CDMS.maps["TTCtr.TERM_A"].outfile
        nameB = self.CDMS.maps["TTCtr.TERM_B"].outfile
        self.ftp.put(self.localpath+nameA, self.remotepath+nameA)
        self.ftp.put(self.localpath+nameB, self.remotepath+nameB)
        stdin, stdout, stderr = self.ssh.exec_command("cdmsLoad -f "+self.remotepath+nameA+" TTCtr.TERM_A --rename")
        while not stdout.channel.exit_status_ready():
            if stdout.channel.recv_ready():
                rl, wl, xl = select.select([stdout.channel], [], [], 0.0)
                if len(rl) > 0:
                    print stdout.channel.recv(1024)
        stdin, stdout, stderr = self.ssh.exec_command("cdmsLoad -f "+self.remotepath+nameB+" TTCtr.TERM_B --rename")
        while not stdout.channel.exit_status_ready():
            if stdout.channel.recv_ready():
                rl, wl, xl = select.select([stdout.channel], [], [], 0.0)
                if len(rl) > 0:
                    print stdout.channel.recv(1024)

    def make_TT_unscr(self):
        self.CDMS.maps["TTCtr.SEC_ACT_UNSCR_MAP"].write(path=self.localpath)
        name = self.CDMS.maps["TTCtr.SEC_ACT_UNSCR_MAP"].outfile
        self.ftp.put(self.localpath+name, self.remotepath+name)
        stdin, stdout, stderr = self.ssh.exec_command("cdmsLoad -f "+self.remotepath+name+" TTCtr.SEC_ACT_UNSCR_MAP --rename")
        while not stdout.channel.exit_status_ready():
            if stdout.channel.recv_ready():
                rl, wl, xl = select.select([stdout.channel], [], [], 0.0)
                if len(rl) > 0:
                    print stdout.channel.recv(1024)

    def set_CommandMatrix(self, pattern):
        self.CDMS.maps["Recn.REC1.CM"].replace(pattern)
        self.CDMS.maps["Recn.REC1.CM"].write(path=self.localpath)
        name = self.CDMS.maps["Recn.REC1.CM"].outfile
        self.ftp.put(self.localpath+name, self.remotepath+name)
        stdin, stdout, stderr = self.ssh.exec_command("cdmsLoad -f "+self.remotepath+name+" Recn.REC1.CM --rename")
        while not stdout.channel.exit_status_ready():
            if stdout.channel.recv_ready():
                rl, wl, xl = select.select([stdout.channel], [], [], 0.0)
                if len(rl) > 0:
                    print stdout.channel.recv(1024)

    def get_InteractionMatrices(self):
        HOname = self.CDMS.maps["HORecnCalibrat.RESULT_IM"].outfile
        TTname = self.CDMS.maps["TTRecnCalibrat.RESULT.IM"].outfile
        command = "cdmsSave -f "+self.remotepath+HOname+" HORecnCalibrat.RESULT_IM"
        self.sendCommand(command)
        command = "cdmsSave -f "+self.remotepath+TTname+" TTRecnCalibrat.RESULT.IM"
        self.sendCommand(command)
        self.ftp.get(self.remotepath+HOname, self.localpath+HOname)
        self.ftp.get(self.remotepath+TTname, self.localpath+TTname)
        return HOname, TTname
        

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
            stdin, stdout, stderr = self.ssh.exec_command("cdmsSetProp Acq.CFG.DYNAMIC DET1.PIXEL_TAP -s \""+tp+"\"")
            while not stdout.channel.exit_status_ready():
                if stdout.channel.recv_ready():
                    rl, wl, xl = select.select([stdout.channel], [], [], 0.0)
                    if len(rl) > 0:
                        print stdout.channel.recv(1024)
        except:
            print("Error!  Invalid tap point!")

    def measureBackground(self, nframes):
        stdin, stdout, stderr = self.ssh.exec_command("msgSend \"\" CommandGateway EXEC \"AcqOptimiser.measureBackground "+str(nframes)+"\"")
        while not stdout.channel.exit_status_ready():
            if stdout.channel.recv_ready():
                rl, wl, xl = select.select([stdout.channel], [], [], 0.0)
                if len(rl) > 0:
                    print stdout.channel.recv(1024)

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
        definitionFile = '/home/deen/Code/Python/BlurryApple/Tools/CDMS_Map_Definitions.dat'
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

