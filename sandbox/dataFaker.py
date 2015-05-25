import scipy
import numpy
import NGCTools
import matplotlib.pyplot as pyplot

detector = NGCTools.detector()

zern = [0.0, 0.0, 0.0, 0.0, 0.0]
zerntip = [1.0, 0.0, 0.0, 0.0, 0.0]
zerntilt = [0.0, 1.0, 0.0, 0.0, 0.0]
pokes = numpy.zeros(60)
pupil = [0.0, 0.0]
angle = 0.0

single = False

if single:
    detector.generateFrame(zern, pupil, pokes, angle)
else:
    detector.generateFrame([0.0, 0.0, 0.0, 0.0, 0.0], pupil, pokes, angle)
    for i in range(21):
        tip = 10.0*numpy.sin(i*6.28/20.0)
        tilt = 10.0*numpy.cos(i*6.28/20.0)
        print tip, tilt
        detector.generateFrame([tip, tilt, 0.0, 0.0, 0.0], pupil, pokes, angle)
    #print 'flat'
    #detector.generateFrame(zern, pupil, pokes, angle)
    #print 'tip'
    #detector.generateFrame(zerntip, pupil, pokes, angle)
    #print 'tilt'
    #detector.generateFrame(zerntilt, pupil, pokes, angle)

fig = pyplot.figure(0)
fig.clear()
if single:
    ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax1.matshow(detector.z[0])
else:
    ax1 = fig.add_axes([0.1, 0.1, 0.4, 0.4])
    ax2 = fig.add_axes([0.1, 0.5, 0.4, 0.4])
    ax3 = fig.add_axes([0.5, 0.1, 0.4, 0.4])
    ax4 = fig.add_axes([0.5, 0.5, 0.4, 0.4])

    ax1.matshow(detector.z[0])
    ax2.matshow(detector.z[1])
    ax3.matshow(detector.z[2])
    ax4.matshow(detector.z[3])

detector.saveRawFrames("TipTiltMotion.fits")

fig.show()
fig.savefig("TTMotion.png")
