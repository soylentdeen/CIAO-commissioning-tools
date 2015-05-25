import scipy
import scipy.fftpack as fftpack
import matplotlib.pyplot as pyplot
import numpy
import pyfits
import NGCTools
import scipy.optimize

detector = NGCTools.detector()
LA = detector.lensletArray

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

makeCircle = lambda p, t: p[0]*numpy.sin(6.28*t/p[1] + p[2])

def fitCircle(TT, time):
   dT = numpy.average(time[1:101] - time[:100])
   yhat = fftpack.rfft(numpy.array(TT[:,0], dtype=numpy.float32))
   idx = (yhat**2).argmax()
   freqs = fftpack.rfftfreq(len(TT[:,0]), d=dT/numpy.pi)

   params = (numpy.max(TT[:,0])-numpy.min(TT[:,0]))/2.0, freqs[idx]/numpy.pi, 0.0

   errfunc = lambda p, t, meas: numpy.r_[
             makeCircle(p, t) - meas[:,0],
             makeCircle(p+[0.0, 0.0, numpy.pi/2.0], t) - meas[:,1]
           ]
   p, success = scipy.optimize.leastsq(errfunc, params, args = (time, TT))
   print p
   return p

scfac = 1.0

i = 0
j = 0
xc = []
yc = []

for i in range(9):
    for j in range(9):
        if LA.apertureMap[i][j]:
            xc.append(i*8+4)
            yc.append(j*8+4)

#datadir = '/diska/data/SPARTA/2015-05-22/'
datadir = '../Data/'

#numbers = [16, 17, 18, 19]
numbers = [25, 26, 27, 28]
#numbers = [3, 4, 5, 6]
amp = [0.1, 0.05, 0.025, 0.0125]
color = ['b', 'g', 'r', 'm']

for n, a, c in zip(numbers, amp, color):
    LoopData = pyfits.getdata(datadir+'LoopData_'+str(n)+'/LoopData_'+str(n)+'.fits')
    
    Centroids = LoopData.field(4)
    TTCommands = LoopData.field(6)
    seconds = LoopData.field(1)
    microsec = LoopData.field(2)
    time = numpy.float32(seconds-numpy.min(seconds) + microsec*1e-6)

    coeffs = fitCircle(TTCommands, time)

    for centx, centy, x, y in zip(Centroids.T[0::2], Centroids.T[1::2], xc, yc):
        #ax.plot(TTCommands[:,0], centx)
        #ax.plot(TTCommands[:,1], centy)
        meas = numpy.array([centx-numpy.mean(centx), 
                            centy-numpy.mean(centy)])
        expected = numpy.array([a*numpy.sin(time/coeffs[1]+coeffs[2]),
                                a*numpy.cos(time/coeffs[1]+coeffs[2])])
        #ax.plot((time% coeffs[1])+x, (meas-expected)+y, color=c)
        ax.plot(scfac*(centx-numpy.mean(centx))+x, scfac*(centy-numpy.mean(centy))+y, color = c)
        #ax.plot(a*10*scfac*numpy.sin(numpy.arange(21)*6.28/20)+x, a*10*scfac*numpy.cos(numpy.arange(21)*6.28/20)+y, 'k')

ax.set_title("Linear Centroiding Estimation Algorithm Testing")
ax.set_xlabel("X difference CoG - LCE (Pixels)")
ax.set_ylabel("Y difference CoG - LCE (Pixels)")
#ax.set_xbound(-5e-8, 1e-7)
#ax.set_ybound(-5e-8, 1e-7)
ax.legend()
fig.show()
fig.savefig("LCE_Test.png")
