import scipy
import VLTTools

datadir = '/diska/data/SPARTA/2015-06-01/temp/'

ciao = VLTTools.VLTConnection(simulate=False, datapath=datadir)

ciao.set_Tip(-0.017)
ciao.set_Tilt(0.03)

ciao.setup_TTIM(cycles=3)
ciao.measure_TTIM(config=True)
ciao.setup_HOIM(cycles=1)
ciao.measure_HOIM(config=True)

ciao.get_InteractionMatrices()
ciao.calc_CommandMatrix(nFiltModes=20)

ciao.set_TT_gain(-0.1)
ciao.measureNewTTRefPositions()
ciao.set_HO_gain(-0.1)
ciao.measureNewHORefPositions()
