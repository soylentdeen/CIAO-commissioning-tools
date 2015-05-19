import scipy
import matplotlib.pyplot as pyplot
import VLTTools
import numpy
import time


datadir = "/diska/data/SPARTA/2015-04-29/TTFocus/"
ciao = VLTTools.VLTConnection(simulate=False, datapath=datadir)

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

Tip = []
Tilt = []
Focus = []

length = 120
tip = numpy.zeros(length)
tilt = numpy.zeros(length)
focus = numpy.zeros(length)

tipPlot, = ax.plot(tip, c='b')
tiltPlot, = ax.plot(tilt, c='g')
focusPlot, = ax.plot(focus, c='r')
fig.show()
pix2m = (3.14159/(180.0*3600.0)) * 8.0*0.5/2.0

n = 13.49


while True:
    TTF = ciao.getTTFocus() * pix2m
    tip = numpy.append(tip,TTF[0]*4*n*1000.0)
    tilt = numpy.append(tilt, TTF[2]*4*n*1000.0)
    focus = numpy.append(focus, TTF[4]*16*(3.0)**(0.5)*n**2.0*1000.0)

    tipPlot.set_ydata(tip[-length:])
    tiltPlot.set_ydata(tilt[-length:])
    focusPlot.set_ydata(focus[-length:])
    mn = numpy.min([numpy.min(tipPlot.get_data()[1]), numpy.min(tiltPlot.get_data()[1]), numpy.min(focusPlot.get_data()[1])])
    mx = numpy.max([numpy.max(tipPlot.get_data()[1]), numpy.max(tiltPlot.get_data()[1]), numpy.max(focusPlot.get_data()[1])])
    ax.set_ybound(lower=mn*1.2, upper=mx*1.2)

    fig.canvas.draw()

    TTF *= 1e9
    # print("focus: %.4f, tip: %.4f, tilt: %.4f" % (focus[-1], tip[-1], tilt[-1])) 
    print("focus: %.0f, tip: %.0f, tilt: %.0f" % (TTF[4], TTF[0], TTF[2]))

    time.sleep(1)
