import pyfits
import matplotlib.pyplot as pyplot
import sys
import numpy
import termios
import tty
import time

class _Getch:
   def __enter__(self):
       self.fd = sys.stdin.fileno()
       self.old_settings = termios.tcgetattr(self.fd)
       try:
          tty.setraw(sys.stdin.fileno())
          ch = sys.stdin.read(1)
       except:
          ch = ''
       #finally:
       #   termios.tcsetattr(self.fd, termios.TCSADRAIN, self.old_settings)
       return ch
   
   def __exit__(self, type, value, traceback):
       termios.tcsetattr(self.fd, termios.TCSADRAIN, self.old_settings)
       
       

#CMFile = "/diska/home/ciaomgr/data/Recn.REC1.CM.fits"
CMFile = "/diska/home/ciaomgr/data/controlMatrix.fits"
CM = pyfits.getdata(CMFile)
TT2HOfile = "/diska/home/ciaomgr/data/HOCtr.TT_TO_HO.fits"
TT2HO = pyfits.getdata(TT2HOfile)
loopFile = "/diska/data/SPARTA/30March2015/"+sys.argv[1]+'/'+sys.argv[1]+'.fits'
loopData = pyfits.getdata(loopFile)

t = loopData.field(1)
gradients = loopData.field(4)
hodm = loopData.field(5)
ttm = loopData.field(6)
actuators = numpy.append(hodm, ttm, axis=1)
st = numpy.std(hodm, axis=0)
print "Max Value: ", numpy.max(numpy.abs(hodm))
fig = pyplot.figure(1)
fig.clear()
fig.show()
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.4])
ax2 = fig.add_axes([0.1, 0.5, 0.8, 0.4])
ActPos1, = ax1.plot(actuators[0,:], c='r')
ActPos2, = ax1.plot(numpy.dot(TT2HO, ttm[0,:]), c='g')
#ActPos3, = ax1.plot(actuators[2,:])
#ActPos4, = ax1.plot(actuators[3,:])
CalcPos, = ax1.plot(actuators[0,:]+numpy.dot(CM, gradients[0,:]), c='b')
Gradx, = ax2.plot(gradients[0,0::2])
Grady, = ax2.plot(gradients[0,1::2])


i = 205
next = 'j'
print("h=+10, j=+1, k=-1, l=-10, q = quit")
while (next != 'q'):
    ActPos1.set_ydata(actuators[i,:]-actuators[i+1,:])
    ActPos2.set_ydata(numpy.dot(TT2HO, ttm[i,:]-ttm[i+1,:]))
    #ActPos3.set_ydata(actuators[i+2,:])
    #ActPos4.set_ydata(actuators[i+3,:])
    CalcPos.set_ydata(numpy.dot(CM, gradients[i,:])*0.3)
    Gradx.set_ydata(gradients[i,0::2])
    Grady.set_ydata(gradients[i,1::2])
    ax2.set_title("Frame Number " + str(i))
    ax1.set_ybound(lower=-0.5, upper=0.5)
    ax2.set_ybound(lower=-1.0, upper=1.0)
    #next = raw_input("j - forward, k - backward, q - quit")
    fig.canvas.draw()
    with _Getch() as rc:
        next = rc
    if next == 'j':
        i += 1
    elif next == 'k':
        i -= 1
    elif next == 'h':
        i += 10
    elif next == 'l':
        i -= 10
#ax1 = fig.add_axes([0.1, 0.1, 0.4, 0.8])
#ax2 = fig.add_axes([0.5, 0.1, 0.4, 0.8])

#ax1.matshow(hodm, aspect='auto')
#ax2.plot(st)
#ax2.plot(hodm[:,5])

#fig.show()
