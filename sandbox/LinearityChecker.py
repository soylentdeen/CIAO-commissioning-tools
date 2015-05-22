import scipy
import matplotlib.pyplot as pyplot
import numpy
import pyfits
from NGCTools import lensletArray
import scipy.optimize

LA = lensletArray()

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

def makeCircle(params, t):
    r = params[0]
    k = params[1]
    lag = params[2]
    x = r*numpy.sin(t/k+lag)
    y = r*numpy.cos(t/k+lag)
    print x[30]
    return numpy.array([x, y]).T
    

def fitCircle(TT, time):
   params = 0.5, 30, 0.1
   errfunc = lambda p: numpy.ravel(makeCircle(p, time) - TT)
   p, success = scipy.optimize.leastsq(errfunc, params)
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

datadir = '/diska/data/SPARTA/2015-05-22/'

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
    time = seconds-numpy.min(seconds) + microsec*1e-6

    coeffs = fitCircle(TTCommands, time)

    for centx, centy, x, y in zip(Centroids.T[0::2], Centroids.T[1::2], xc, yc):
        #ax.plot(TTCommands[:,0], centx)
        #ax.plot(TTCommands[:,1], centy)
        meas = numpy.array([centx-numpy.mean(centx), 
                            centy-numpy.mean(centy)])
        expected = numpy.array([a*numpy.sin(time/coeffs[1]+coeffs[2]),
                                a*numpy.cos(time/coeffs[1]+coeffs[2])])
        print asdf
        ax.plot((time% coeffs[1])+x, (meas-expected)+y, color=c)
        #ax.plot(scfac*(centx-numpy.mean(centx))+x, scfac*(centy-numpy.mean(centy))+y, color = c)
        #ax.plot(a*10*scfac*numpy.sin(numpy.arange(21)*6.28/20)+x, a*10*scfac*numpy.cos(numpy.arange(21)*6.28/20)+y, 'k')

ax.set_title("Linear Centroiding Estimation Algorithm Testing")
ax.set_xlabel("X difference CoG - LCE (Pixels)")
ax.set_ylabel("Y difference CoG - LCE (Pixels)")
#ax.set_xbound(-5e-8, 1e-7)
#ax.set_ybound(-5e-8, 1e-7)
ax.legend()
fig.show()
fig.savefig("LCE_Test.png")
