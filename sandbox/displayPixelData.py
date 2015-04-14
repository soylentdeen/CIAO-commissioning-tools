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
       
fig = pyplot.figure(1)
fig.clear()
fig.show()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
imagePlot = ax.matshow(numpy.zeros((72,72)), vmin=0, vmax=300)
i = 0
next = 'j'
print("j=+1, k=-1, q = quit")
while (next != 'q'):
    pixelFile = "/home/deen/Data/GRAVITY/2015-04-10/Derotator/derot_circbuff_"+str(i)+"_/derot_circbuff_"+str(i)+"_.fits"
    pixelData = pyfits.getdata(pixelFile)

    t = pixelData.field(1)
    pixels = pixelData.field(3)
    image = numpy.average(pixels, axis=0)
    imagePlot.set_clim(vmin= 0.0, vmax = numpy.max(image)/1.2)
    imagePlot.set_data(image.reshape(72,72))
    fig.canvas.draw()
    print i
    with _Getch() as rc:
        next = rc
    if next == 'j':
        i += 1
    elif next == 'k':
        i -= 1


#ax1.matshow(hodm, aspect='auto')
#ax2.plot(st)
#ax2.plot(hodm[:,5])

#fig.show()
