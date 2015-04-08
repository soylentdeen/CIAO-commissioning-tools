#!/usr/bin/env python

import logging
import subprocess
import os


logging.basicConfig(level= logging.DEBUG)

def execCmd(cmd):
    logging.debug("Executing '%s'" % cmd)
    exitCode= os.system(cmd)
    assert exitCode == os.EX_OK


class Request(object):

    def __init__(self, config):
        self.dc9= config["DC9"]
        self.dit= config["DIT"]
        self.fileBaseName= config["FILE_BASE_NAME"]
        self.ndsamples= config["NDSAMPLES"]
        self.nsamppix= config["NSAMPPIX"]
        self.ndit= config["NDIT"]
        self.readoutModeName= config["READOUT_MODE"]


class DarkFrameTaker(object):

    PARAM_FILTER_POS_NAME= "INS.FILT1.NAME"

    def __init__(self, config):
        self._envName= config["environment"]
        self._icsProcName= config["icsProcess"]
        self._ngcProcName= config["ngcProcess"]
        self._fwClosedPosName= config["filterWheelClosedPositionName"]


    def moveFilterWheelToClosedPosition(self):
        logging.info("Moving filter wheel to closed position ...")
        execCmd("msgSend %s %s SETUP \"-function %s %s\"" % (
             self._envName, self._icsProcName, self.PARAM_FILTER_POS_NAME,
             self._fwClosedPosName))
        logging.info("Filter wheel has been moved to closed position.")


    def getFilterWheelPositionName(self):
        proc= subprocess.Popen(
            ["msgSend", self._envName, self._icsProcName, "STATUS",
             "-function " + self.PARAM_FILTER_POS_NAME],
            stdout= subprocess.PIPE,
            stderr= subprocess.STDOUT)
        proc.wait()
        assert os.EX_OK == proc.returncode

        outputLines= []
        for each in proc.stdout:
            outputLines.append(each)
            logging.debug(each)
            if each.startswith(self.PARAM_FILTER_POS_NAME):
                spacePos= each.find(' ')
                return each[spacePos+1:-1].strip()
        raise Exception("No filter wheel position: %s" '\n'.join(outputLine))


    def takeDarkFrame(self, req):
        logging.info("Taking dark frame ...")
        execCmd("msgSend %s %s SETUP \"-function DET.READ.CURNAME %s\"" % (
            self._envName, self._ngcProcName, req.readoutModeName))
       	execCmd("msgSend %s %s SETUP \"-function DET.SEQ1.DIT %f\"" % (
            self._envName, self._ngcProcName, req.dit))
        execCmd("msgSend %s %s SETUP \"-function DET.NDIT %d\"" % (
            self._envName, self._ngcProcName, req.ndit))
        execCmd("msgSend %s %s SETUP \"-function DET.NDSAMPLES %d\"" % (
            self._envName, self._ngcProcName, req.ndsamples))
        execCmd("msgSend %s %s SETUP \"-function DET.NSAMPPIX %d\"" % ( 
            self._envName, self._ngcProcName, req.nsamppix))
        execCmd("msgSend %s %s SETUP \"-function DET.CLDC1.DC9 %f\"" % (
            self._envName, self._ngcProcName, req.dc9))

        logging.info("Listing available readout modes ...")
        execCmd("msgSend %s %s STATUS \"-function DET.READ.AVAIL\"" % (
            self._envName, self._ngcProcName))
        logging.info("Available readout modes have been listed.")

        execCmd("msgSend %s %s SETUP \"-function DET.FRAM.FILENAME %s\"" % (
            self._envName, self._ngcProcName, req.fileBaseName))
        execCmd("msgSend %s %s SETUP \"-function DET.FRAM1.GEN T\"" % (self._envName, self._ngcProcName))
        execCmd("msgSend %s %s SETUP \"-function DET.FRAM1.STORE T\"" % (self._envName, self._ngcProcName))
        execCmd("msgSend %s %s SETUP \"-function DET.FRAM3.GEN T\"" % (self._envName, self._ngcProcName))
        execCmd("msgSend %s %s SETUP \"-function DET.FRAM3.STORE T\"" % (self._envName, self._ngcProcName))
        # TODO execCmd("msgSend %s %s FRAME \"-gen T\"" % (
        #    self._envName, self._ngcProcName))

        execCmd("msgSend %s %s START \"\"" % (self._envName, self._ngcProcName))
        execCmd("msgSend %s %s WAIT \"\"" % (self._envName, self._ngcProcName))

        # execCmd("msgSend %s %s START "-expo0w
        logging.info("Dark frame has been taken.")

config= {
    "environment": "wci1ao",
    "icsProcess": "ciiControl",
    "filterWheelClosedPositionName": "CLOSED1",
    "ngcProcess": "ngcircon_NGCIR1",
 
}


taker= DarkFrameTaker(config)
taker.moveFilterWheelToClosedPosition()
fwPosName= taker.getFilterWheelPositionName()
print "FW pos name:", fwPosName
assert config["filterWheelClosedPositionName"] == fwPosName

print "OK1"; import sys; sys.exit(1)

"""
for eachDC9 in xrange(-6, 6): 
    taker.takeDarkFrame(Request({"DC9": eachDC9,
                                 "NDIT": 10,
                                 "DIT": 0.005,
                                 "FILE_BASE_NAME": "/diska/insroot/CIAO/SYSTEM/DETDATA/dark_fow8-2_COMM%d" % eachDC9,
                                 "NDSAMPLES": 8,
                                 "NSAMPPIX": 4,
                                 
				 "READOUT_MODE": "Fowler", 
				 #"READOUT_MODE": "FowlerNsamp_Windows", 
                                 #"READOUT_MODE": "Double",
                                }))
"""
