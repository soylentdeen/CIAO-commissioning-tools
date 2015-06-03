import scipy
import VLTTools

datadir = '/diska/data/SPARTA/2015-06-01/temp/'

ciao = VLTTools.VLTConnection(simulate=False, datapath=datadir)

ciao.measureBackground(9)
