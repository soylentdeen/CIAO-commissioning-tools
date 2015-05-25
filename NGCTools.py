import numpy
import scipy
import pyfits
import os
from scipy.misc import factorial as fac
import scipy.interpolate as interp
import scipy.fftpack as fftpack
import matplotlib.pyplot as pyplot

def twoDgaussian(x, y, center, stdev, A):
    retval = A * (numpy.exp(-(x-center[0])**2.0/stdev[0])*
                  numpy.exp(-(y-center[1])**2.0/stdev[1]))
    print center, A, numpy.max(retval), numpy.max(x), numpy.min(x)
    return retval

class zernikeMode(object):
    """
       Class representing a Zernike mode
    """
    def __init__(self, noll, mag):
        """
        input:  noll - Noll index
                mag - magnitude of Zernike mode - Units?
        """
        self.mag = mag
        if (noll == 2):
            self.n = 1
            self.m = 1
        elif (noll == 3):
            self.n = 1
            self.m = -1
        elif (noll == 4):
            self.n = 2
            self.m = 0
        elif (noll == 5):
            self.n = 2
            self.m = -2
        elif (noll == 6):
            self.n = 2
            self.m = 2
        else:
            self.n = 0
            self.m = 0

    def zernike_rad(self, rho):
        n=abs(self.n)
        m=abs(self.m)
        if (numpy.mod(n-m, 2) == 1):
            return rho*0.0
        
        wf = rho*0.0
        for k in range((n-m)/2+1):
            wf += rho**(n-2.0*k) * (-1.0)**k * fac(n-k) / ( fac(k) * fac( (n+m)/2.0 - k ) * fac( (n-m)/2.0 - k ) )
        
        return wf

    def setMag(self, mag):
        self.mag = mag

    def zernike(self, rho, phi, norm=False):
        nc = self.mag
        if (norm):
            nc = (2*(self.n+1)/(1+(self.m==0)))**0.5
        if (self.m > 0): return nc*self.zernike_rad(rho) * numpy.cos(self.m*phi)
        if (self.m < 0): return nc*self.zernike_rad(rho) * numpy.sin(-self.m*phi)
        return nc*self.zernike_rad(rho)

class pupil( object ):
    def __init__(self, x, y, innerRadius, outerRadius):
        self.x = x
        self.y = y
        self.inrad2 = innerRadius**2.0
        self.outrad2 = outerRadius**2.0

    def calculateFlux(self, x, y):
        pythag = (self.x-x)**2.0 + (self.y - y)**2.0
        if ( (pythag < self.outrad2) & (pythag > self.inrad2)):
            return 1.0
        else:
            return 0.0

    def calculateApertureIllumination(self, corner, size):
        area = 0.0
        illuminated = 0.0
        for x in numpy.linspace(corner[0], corner[0]+size, 50):
            for y in numpy.linspace(corner[1], corner[1]+size, 50):
                area += 1.0
                pythag = (self.x-x)**2.0 + (self.y - y)**2.0
                if ( (pythag < self.outrad2) & (pythag > self.inrad2)):
                    illuminated += 1.0
        return illuminated/area

    def getDecenter(self):
        return self.x, self.y

    def setDecenter(self, x, y):
        self.x = x
        self.y = y

class deformableMirror( object ):
    def __init__(self, parent):
        self.nActuators = 60
        self.influenceFunctions = pyfits.getdata(parent.datadir+'IF_cube_HR.fits')
        nx = len(self.influenceFunctions[0][0])
        ny = len(self.influenceFunctions[0])
        xcoords = numpy.linspace(-888.0, 888.0, nx)
        ycoords = numpy.linspace(-888.0, 888.0, ny)
        self.interpFunctions = []
        self.actuatorPositions = numpy.zeros(len(self.influenceFunctions))
        for inflfunc in self.influenceFunctions:
            self.interpFunctions.append(interp.interp2d(xcoords, ycoords,
                            inflfunc, kind='cubic'))
        
    def setMirror(self, actPos):
        self.actuatorPositions = actPos

    def calcPosition(self, x, y):
        retval = 0.0
        for IF, amp in zip(self.interpFunctions, self.actuatorPositions):
            retval += amp*IF(x, y)[0]
        return retval

class derotator( object ):
    def __init__(self):
        self.angle = 0.0

    def setAngle(self, angle):
        self.angle = angle
        self.sine = numpy.sin(self.angle)
        self.cosine = numpy.cos(self.angle)

    def getMirrorPosition(self, x, y):
        newx = x*self.cosine-y*self.sine
        newy = x*self.sine+y*self.cosine
        return newx, newy


class detector( object ):
    """
        The Detector class allows us to simulate the detector in the cryostat.
        Some of the numbers are most likely wrong, but with tweaking we should
        be able to simulate the effects of simple Zernike aberations of the 
        wave front on the spots on the detector.
    """
    def __init__(self, beamSize=1776.0):
        # Calculates the relative size of the central obscuration (M2)
        self.datadir = os.path.dirname(__file__)+'/Data/'
        self.scramblingMap = pyfits.getdata(
                            self.datadir+"scramblemap.fits")
        self.unscramblingMap = pyfits.getdata(
                            self.datadir+"unscramblemap.fits")
        self.windowMap = pyfits.getdata(
                            self.datadir+"windowmap.fits")
        self.SLsubapmap = self.datadir+"LoopDisplaySrv.SUBAP_MAP.fits"
        self.pixPerSubAp = 8
        self.centObscScale = 1.116/8.00
        self.beamSize = beamSize
        self.lenslet = lensletArray(self)
        self.wavefront = waveFront(beamSize=beamSize)
        self.pupil = pupil(0.0, 0.0,
                     innerRadius=beamSize/2.0*self.centObscScale,
                     outerRadius=beamSize/2.0)
        self.DM = deformableMirror(self)
        self.derotator = derotator()
        self.nx = 72
        self.ny = 72
        self.readoutMode = "8x8"
        self.spacing = 24.0 #microns  Is this value correct?
        self.xpix = (numpy.arange(self.nx)-self.nx/2.0)*self.spacing
        self.ypix = (numpy.arange(self.ny)-self.ny/2.0)*self.spacing
        self.stdev = (8.0*self.spacing, 8.0*self.spacing)
        self.flux = 1000.0
        self.z = []
        self.frames = []
        self.centroids = []
        #self.weightingMap = self.makeWeightingMap()

    def scrambleFrame(self):
        #"""
        scrambledFrame = numpy.zeros(6912)
        flatframe = self.z[-1].ravel()
        for y, i in zip(flatframe, self.scramblingMap):
            scrambledFrame[i]=y

        self.frames.append(scrambledFrame)
        """

        scrambledFrame = numpy.zeros(5184)
        HRframe = self.z[-1].ravel()
        for y, i in zip(HRframe, self.unscramblingMap):
            scrambledFrame[i] = y
        #"""

        """
        largeScrambledFrame = numpy.zeros(6912)
        j = 0
        for i in range(6912):
            if self.windowMap[i] == 1:
                largeScrambledFrame[i] = scrambledFrame[j]
                j += 1


        self.frames.append(largeScrambledFrame)
        #"""
        

    def makeRamp(self):
        z = numpy.zeros((self.ny, self.nx))
        k = 0
        for i in range(self.nx):
            for j in range(self.ny):
                z[i][j] = k
                k += 1
        self.z.append(z)
        self.scrambleFrame()
        self.centroids.append([])

    def generateFrame(self, zern, pupil, actuatorPokes, angle):
        """
        Generates an image seen by the detector of a wavefront described by
        the zernike coefficients in zern

        zern = [tip, tilt, defocus, astig1, astig2]
        """
        self.derotator.setAngle(numpy.deg2rad(angle))
        self.pupil.setDecenter(pupil[0], pupil[1])
        self.wavefront.setZern(zern)
        self.DM.setMirror(actuatorPokes)
        self.expose()
        """
        centroids, amplitudes = self.calculateCentroids(zern, actuatorPokes)
        z = numpy.zeros((self.ny, self.nx))
        for xcoord in range(self.nx):
            for ycoord in range(self.ny):
                z[ycoord][xcoord] = self.flux*sum(amplitudes*
                   numpy.exp(-(self.xpix[xcoord]+self.spacing/
                       2.0-centroids[:,0])**2.0/self.stdev[0])*
                   numpy.exp(-(self.ypix[ycoord]+self.spacing/
                       2.0-centroids[:,1])**2.0/self.stdev[1]))
        self.z.append(z)
        self.scrambleFrame()
        self.centroids.append(centroids)
        #"""

    def expose(self):
        subsamp = 100.0
        FWHM_i = 1.1 * subsamp # FWHM in Pixel space
        FWHM_k = 0.88*self.pixPerSubAp*subsamp/FWHM_i
        location = int(round((self.pixPerSubAp*subsamp - FWHM_k)/2.0))
        delta = self.lenslet.spacing/2.0
        nPixPoints = self.pixPerSubAp*subsamp
        z = numpy.zeros((self.ny, self.nx))
        ptsOnDetector = numpy.linspace(-nPixPoints/2.0*self.spacing, 
                nPixPoints*self.spacing/2.0, nPixPoints)
        gridx, gridy = numpy.meshgrid(ptsOnDetector, ptsOnDetector)
        gridx /= 100.0
        gridy /= 100.0
        totalFlux = 0.0
        #fig = pyplot.figure(0)
        #fig.clear()
        #ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        for coord in self.lenslet.coordinates:
            flux = numpy.zeros((nPixPoints, nPixPoints), dtype=numpy.complex)
            #intensity = numpy.zeros((nPixPoints, nPixPoints))
            #phase = numpy.zeros((nPixPoints, nPixPoints))
            i = location
            for x in numpy.linspace(coord[0]-delta, coord[0]+delta, 
                    num=int(round(FWHM_k))):
                j = location
                for y in numpy.linspace(coord[1]-delta, coord[1]+delta, 
                        num=int(round(FWHM_k))):
                    flux[i][j] = numpy.complex(self.pupil.calculateFlux(x, y),
                            self.calcWaveFront(x, y))
                    #intensity[i][j] = self.pupil.calculateFlux(x, y)
                    #phase[i][j] = self.calcWaveFront(x, y)
                    totalFlux += 1
                    j += 1
                i += 1
            image = fftpack.fft2(flux)
            image = fftpack.fftshift((image*image.conjugate()).real)
            xc = coord[0]+gridx
            yc = coord[1]+gridy
            inviewx = (self.xpix > numpy.min(xc)) & (self.xpix < numpy.max(xc))
            inviewy = (self.ypix > numpy.min(yc)) & (self.ypix < numpy.max(yc))
            for i in range(self.nx):
                if inviewx[i]:
                    for j in range(self.ny):
                        if inviewy[j]:
                            fp = scipy.where( (xc >= self.xpix[i]) & (xc < self.xpix[i]+self.spacing) & (yc >= self.ypix[j]) & (yc < self.ypix[j]+self.spacing))
                            z[i][j] = numpy.sum(image[fp])
                            #print i, j, z[i][j], coord
                    #raw_input()
        self.z.append(z)
        self.scrambleFrame()


    def calculateCentroids(self, zern, actuatorPokes):
        """
            Calcualates the locations of the centroids under the given 
            Zernike coefficients
        """
        self.wavefront.setZern(zern)
        self.DM.setMirror(actuatorPokes)
        dx = 10.0   # Microns
        dy = 10.0   # Microns
        centroids = []
        intensities = []
        DCx, DCy = self.pupil.getDecenter()  # Decenter of Pupil
        for c in self.lenslet.coordinates:
            # Calculates the partial derivatives
            zxp = self.calcWaveFront(c[0]+DCx+dx, c[1]+DCy)
            zxm = self.calcWaveFront(c[0]+DCx-dx, c[1]+DCy)
            zyp = self.calcWaveFront(c[0]+DCx, c[1]+DCy+dy)
            zym = self.calcWaveFront(c[0]+DCx, c[1]+DCy-dy)
            delx = (zxp - zxm)/(2)
            dely = (zyp - zym)/(2)

            # Computes the normal vector to the surface
            normalx = -delx*dy
            normaly = -dely*dx
            normalz = dx*dy

            #Calculates the shift in microns on the detector
            theta_x = scipy.arctan2(normalx, normalz)
            theta_y = scipy.arctan2(normaly, normalz)
            shift_x = scipy.tan(theta_x)*self.lenslet.fl
            shift_y = scipy.tan(theta_y)*self.lenslet.fl
            centroids.append([c[0]+shift_x, c[1]+shift_y])

            intensities.append(self.pupil.calculateApertureIllumination(
                               [c[0]-self.lenslet.spacing/2.0, 
                               c[1]-self.lenslet.spacing/2.0],
                               self.lenslet.spacing))

        return numpy.array(centroids), numpy.array(intensities)

    def makeWeightingMap(self):
        weightingMap = numpy.zeros((72,72),dtype=numpy.float32)

        for coord in self.lenslet.coords:
            x = coord[0]
            y = coord[1]
            print "Not done yet!"
                    
        
    def calcWaveFront(self, x, y):
        wave = self.wavefront.calcWaveFront(x, y)
        rotatedPosition = self.derotator.getMirrorPosition(x, y)
        mirror = self.DM.calcPosition(rotatedPosition[0], rotatedPosition[1])
        return wave + mirror

    def saveRawFrames(self, filename):
        """
        Saves the raw (unscrambled) frames to a data file
        """

        self.z = numpy.array(self.z)
        hdu=pyfits.PrimaryHDU(self.z)
        hdu.writeto(filename, clobber=True)

    def saveFrames(self, filename):
        """
        Saves the frames to a SPARTA-readable data file
        """
        self.frames = numpy.array(self.frames)
        hdu = pyfits.PrimaryHDU(self.frames)
        hdu.scale('int16', bzero=32768, bscale=1)
        hdu.writeto(filename, clobber=True)

    def saveCentroids(self, filename):
        """
        Saves the calculated centroids to a fits file
        """
        self.centroids = numpy.array(self.centroids)
        hdu = pyfits.PrimaryHDU(self.centroids)
        hdu.writeto(filename, clobber=True)

class lensletArray( object ):
    """
    This class simulates the lenslet array
    """
    def __init__(self, parent, spacing=192.0, fl=2095.0):
        """
            Spacing - spacing between adjacent lenslets (in microns)
            fl - focal length of individual lenslet (in microns)
        """
        self.parent = parent
        self.SLapertureMap = pyfits.getdata(self.parent.SLsubapmap)
        """
        self.apertureMap = [[False,False,True,True,True,True,True,False,False],
               [False, True, True, True, True, True, True, True, False],
               [True, True, True, True, True, True, True, True, True],
               [True, True, True, True, True, True, True, True, True],
               [True, True, True, True, False, True, True, True, True],
               [True, True, True, True, True, True, True, True, True],
               [True, True, True, True, True, True, True, True, True],
               [False, True, True, True, True, True, True, True, False],
               [False, False, True, True, True, True, True, False, False]]
        """
        self.spacing = spacing   # Lenslet Array Spacing in microns
        self.fl = fl

        coords = []

        """
        for i in range(9):
            for j in range(9):
                if self.apertureMap[i][j]:
                    coords.append(((i-4)*spacing, (j-4)*spacing))
        """
        for i in range(9):
            for j in range(9):
                if self.SLapertureMap[i][j]:
                    coords.append(((i-4)*spacing, (j-4)*spacing))

        self.coordinates = coords


class waveFront( object ):
    """
    This object describes the wavefront as it hits the lenslet array
    """
    def __init__(self, beamSize = 1776.0):
        self.beamSize = beamSize
        self.tip = zernikeMode(2, 0.00)
        self.tilt = zernikeMode(3, 0.00)
        self.defocus = zernikeMode(4, 0.00)
        self.astig1 = zernikeMode(5, 0.00)
        self.astig2 = zernikeMode(6, 0.00)
    
    def setZern(self, zern):
        """
            Sets the magnitudes of the Zernike components.
        """
        for mag, z in zip(zern,
                [self.tip, self.tilt, self.defocus, self.astig1, self.astig2]):
            z.setMag(mag)

    def calcWaveFront(self, x, y):
        """
        Sums through the different Zernike components at a particular location
        on the wavefront (x, y) to find the local zernike magnitude.
        """
        rho = (x**2.0 + y**2.0)**(0.5)/(self.beamSize/2.0)
        phi = numpy.arctan2(y, x)
        value = 0.0
        for zern in [self.tip, self.tilt, self.defocus,self.astig1,self.astig2]:
            value += zern.zernike(rho, phi)
        return value

    def pokeActuator(self, x, y, InflFunc):
        """
            interpolates the value of an actuator poke at a certain x and y
        """
        return InflFunc(x, y)


"""
class frameBuffer( object ):
    def __init__(self, xoffset = 4.0, yoffset = 4.0):
        self.nx = 96
        self.ny = 72
        self.apertureSize = 8.0
        self.xoffset = xoffset
        self.yoffset = yoffset
        self.xpos = numpy.arange(len(self.apertureMap[0]))*self.apertureSize+self.xoffset
        self.ypos = numpy.arange(len(self.apertureMap))*self.apertureSize+self.yoffset
        self.xpix = numpy.arange(self.nx)
        self.ypix = numpy.arange(self.ny)
        self.xx, self.yy = numpy.meshgrid(self.xpix, self.ypix)
        self.stdev = (2.0, 2.0)
        self.amplitude = 1000.0
        self.frames = []
        self.centroids = []  

    def generateZernikeFrame(self, coeffs=None):
        if not(coeffs):
            coeffs = numpy.array([0.0, 0.0, 0.0, 0.0])

        for zern in coeffs:
            calculate_shifts

    def generateRandomFrames(self, nframes=1):
        for i in range(nframes):
            z = numpy.zeros((self.ny, self.nx))
            centroids = []
            for apvec, y in zip(self.apertureMap, self.ypos):
                for ap, x in zip(apvec, self.xpos):
                    if ap:
                        centroids.append((numpy.random.rand(2)-0.5)*
                                self.apertureSize/4.0+(x, y))

            for centroid in centroids:
                z += twoDgaussian(self.xx, self.yy, centroid, self.stdev,
                        self.amplitude)

            self.addFrame(z, centroids)

    def addFrame(self, frame, cents):
        self.frames.append(frame.ravel())
        self.centroids.append(cents)
        

    def saveFile(self, filename):
        self.frames = numpy.array(self.frames)
        hdu = pyfits.PrimaryHDU(self.frames)
        hdu.scale('int16', bzero=32768, bscale=1)
        hdu.writeto(filename, clobber=True)

    def saveCentroids(self, filename):
        self.centroids = numpy.array(self.centroids)
        hdu = pyfits.PrimaryHDU(self.centroids)
        hdu.writeto(filename, clobber=True)

"""

