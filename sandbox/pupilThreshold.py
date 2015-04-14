import pyfits
import matplotlib.pyplot as pyplot
import sys
import numpy
import termios
import tty
import time

def pasteSubApertures(subAps):
    image = numpy.zeros((72,72))
    apertureMap = [[False,False,True,True,True,True,True,False,False],
               [False, True, True, True, True, True, True, True, False],
               [True, True, True, True, True, True, True, True, True],
               [True, True, True, True, True, True, True, True, True],
               [True, True, True, True, False, True, True, True, True],
               [True, True, True, True, True, True, True, True, True],
               [True, True, True, True, True, True, True, True, True],
               [False, True, True, True, True, True, True, True, False],
               [False, False, True, True, True, True, True, False, False]]
    k = 0
    for i in range(9):
        for j in range(9):
            if apertureMap[i][j]:
                stamp = subAps[k]
                image[i*8+1:(i+1)*8-1, j*8+1:(j+1)*8-1] = stamp
                #image[i*8+1:(i+1)*8-1, j*8+1:(j+1)*8-1] = numpy.floor(stamp)
                k += 1

    return image

def cutoutSubApertures(pixels):
    image = pixels.reshape(72,72)/numpy.median(pixels)
    apertureMap = [[False,False,True,True,True,True,True,False,False],
               [False, True, True, True, True, True, True, True, False],
               [True, True, True, True, True, True, True, True, True],
               [True, True, True, True, True, True, True, True, True],
               [True, True, True, True, False, True, True, True, True],
               [True, True, True, True, True, True, True, True, True],
               [True, True, True, True, True, True, True, True, True],
               [False, True, True, True, True, True, True, True, False],
               [False, False, True, True, True, True, True, False, False]]
    retval = numpy.zeros((68, 6, 6))
    k = 0
    for i in range(9):
        for j in range(9):
            if apertureMap[i][j]:
                stamp = image[i*8+1:(i+1)*8-1,j*8+1:(j+1)*8-1]
                #stamp /= numpy.median(stamp)
                retval[k,:,:] += stamp
                k += 1
            
    return retval


pupilImages = numpy.zeros((10, 68, 6, 6))
threshold = 0.1
for i in range(10):
    pixelFile = "/home/deen/Data/GRAVITY/2015-04-10/DerotatorPupilTest/derot_circbuff_"+str(i)+"_/derot_circbuff_"+str(i)+"_.fits"
    pixelData = pyfits.getdata(pixelFile)

    pixels = pixelData.field(3)
    image = numpy.average(pixels, axis=0)
    pupilImages[i,:,:,:] = cutoutSubApertures(image)

errors = numpy.std(pupilImages, axis=0)
pupils = numpy.average(pupilImages, axis=0)
minflux = numpy.min(pupilImages, axis=0)
maxflux = numpy.max(pupilImages, axis=0)

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
pupilPlot = ax.matshow(pasteSubApertures(maxflux-minflux))
fig.show()
fig.savefig("DerotatorEffectOnPupil.png")
