# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 15:58:38 2014

@author: egendron
"""


import math
import numpy as np
import gral_zer
import pyfits



def rotateCoupleOfCoeffs(a, b, m, angledeg):
    """
    Inputs:
    <a>        : scalar, floating point. Coefficient of 1st zernike mode.
    <b>        : scalar, floating point. Coefficient of 2nd zernike mode.
    <m>        : scalar, integer (azimutal order of the 2 zernike)
    <angledeg> : scalar, floating point (rotation angle)

    Output:
    An array of 2 floating point values: the 2 output coefficient of rotated wavefront.

    For a couple [a,b] of Zernike coefficients with same azimutal order
    m, returns the couple of coefficients of Zernike decomposition
    rotated by an angle <angledeg> (angle expressed in degrees)
    """
    anglerad = angledeg * math.pi / 180.0
    aa = a*np.cos(m*anglerad) - b*np.sin(m*anglerad)
    bb = a*np.sin(m*anglerad) + b*np.cos(m*anglerad)
    return (aa,bb)






def checkConsistencyOfNumberOfElementsOfZernikeArray( imax ):
    """
    Input:
    <imax> : scalar, integer. A number of Zernike mode.

    Output:
    An boolean: True=error, False=no error

    The function returns a boolean error code.
    True when there is an error, False when there is no error.

     """
    # determine radial and azimutal orders n and m
    n,m = gral_zer.nm(imax)
    if( m==0 ):
        return False
    n2 = n//2    # Warning : interger division !  5//2 = 2, and not 2.500 !!
    if( (imax&1)==(n2&1) ):
    # If parity of <imax> and <n2> is the same: not good.
        return True
    else:
        # If parity is different: ok this is correct.
        return False






def rotateZernike(z, i0, angledeg):
    """

    Inputs:
    <z>        : 1D array of floating point values (zernike coefficients)
    <i0>       : scalar, integer. It is the zernike index of the 1st coeff contained in <z>.
    <angledeg> : scalar, floating point

    Output:
    An array of floating-point values (zernike coefficients of the rotated wavefront).

    For a couple [a,b] of Zernike coefficients with same azimutal order
    m, returns the couple of coefficients of Zernike decomposition
    rotated by an angle <angledeg> (angle expressed in degrees)

    """
    # number of coefficients contained in the vector z
    nzer = z.shape[0]

    # index of the last Zernike coefficient
    imax = nzer+i0-1
    # check whether <imax> is ok or not
    errorImax = checkConsistencyOfNumberOfElementsOfZernikeArray( nzer+i0-1 )
    # produces an error if not consistent
    if( errorImax ):
        print "The number of Zernike modes (%d to %d) does not entierely fill a radial order." % (i0, imax)
        print "A coefficient is missing : it is impossible to rotate the wavefront properly."
        print "I will ignore the last Zernike coefficient."
        nzer = nzer-1


    # allocates memory for the result: an array of <nzer> coefficient, same length as <z>, floats.
    z_output = np.zeros(nzer)

    k = 0
    while(k<nzer):
        i  = k+i0      # index of the Zernike mode, as described in the literature (Z_1 is piston)
        # determine radial and azimutal orders, called <n> and <m>
        n,m = gral_zer.nm(i)
        if( m==0 ):
            # do nothing, coefficient is unchanged
            z_output[k] = z[k]
        else:
            if( (i&1)==0 ):    # equivalent to "i modulo 2"
                # if i is an even number (2,4,6,..)
                tmp = rotateCoupleOfCoeffs(z[k], z[k+1], m, angledeg)
                z_output[k]   = tmp[0]
                z_output[k+1] = tmp[1]
            else:
                # if i is an odd number (1,3,5,7,...)
                tmp = rotateCoupleOfCoeffs(z[k+1], z[k], m, angledeg)
                # warning: SWAP coefficients !!!!
                z_output[k]   = tmp[1]
                z_output[k+1] = tmp[0]
            # skip the next coefficient z[k+1], that has already been processed
            k = k+1
        k = k+1

    return z_output












def flipZernike(z, i0):
    """

    Input parameters:
    <z> : 1D array of zernike coefficients
    <i0> : Zernike index of the first element of the array

    Output parameters:
    1D array of zernike coefficients

    The function returns a list of zernike coefficient that represent a
    flipped version of the wavefront wrt to the input. The flip is
    operated around the X axis (equation: y=0, i.e. "horizontal" axis).

    Note: the Yorick version provides a new array of values in
    output. For a C version of the code, it would be MUCH better to do
    the transformation in place, just modifying the coefficients in the
    same memory area, as the transformation is straightforward (just
    modifying the sign of some of the coefficients).
    """
    # number of coefficients contained in the vector z
    nzer = z.shape[0]

    # allocates memory for the result: an array of <nzer> coefficient, same length as <z>, floats.
    z_output = np.zeros(nzer)

    for k in range(nzer):
        i= k+i0   # index of the Zernike mode, as described in the literature (Z_1 is piston)
        # determine radial and azimutal orders, called <n> and <m> (well, actually, n is not required here...)
        n,m = gral_zer.nm(i)
        if( m==0 ):
            z_output[k] = z[k]   # coeff unchanged
        else:
            if( (i&1)==1 ) :
                # if i is an odd number (1,3,5..)
                z_output[k] = -z[k]   # change sign
            else :
                # if i is an even number (2,4,6..)
                z_output[k] = z[k]   # coeff unchanged

    return z_output








def transformZernikeIntoSlopes(z, i0, flip, angledeg, scalingFactor, pixelSizeArcsec, DiamVLT):
    """

    Input parameters:
    <z>    : 1D array of floating-point zernike coefficients.
             Zernike coefficients are expected to be expressed in "nanometers rms".
    <i0>   : scalar integer, Zernike index of the first element of the array
    <flip> : boolean
    <angledeg>        : scalar, floating point. Angle of rotation in degrees.
    <scalingFactor>   : scalar floating point, nominally equal to 1.000
    <pixelSizeArcsec> : scalar floating point = pixel size of WFS in arcsec = 0.51 (Pixel size of WFS is 0.51'')
    <DiamVLT>         : scalar floating point = VLT diameter in meters = 8.0 (VLT is 8m diameter)

    Output parameters:
    1D array of slopes.

    """

    # Reading the Z2S somewhere .. (the path needs to be changed according to your config)
    # <Z2S>  : transformation matrix from zernike to slopes provided by E. Gendron as a FITS
    # file, that should be part of the system configuration
    Z2S = pyfits.getdata("/home/egendron/Projets/GRAVITY/Z2S_136s_119z.fits").T

    # number of coefficients contained in the vector z
    nzer = z.shape[0]
    # Zernike index of the last coefficient of z
    imax = nzer-1+i0

    # number of zernike coefficients contained in the transformation matrix
    nzermat = Z2S.shape[1]
    # Maximum Zernike index in the transformation matrix
    imaxmat = nzermat-1+i0

    if( imax>imaxmat ):
        # hopefully, this should never happen ..
        print "The list of Zernike given in input is too large (%d modes), the transformation matrix only contains %d modes.\n" % (imax,imaxmat)
        print "I will skip all the zernike coeffs greater than %d." % imaxmat

    # allocate memory, 1D array of floating points filled with 0.00
    z_out = np.zeros(nzermat)

    # Copies zernike coefficients in another array that will match the dimension of the Z2S
    # with coeffs at the right place
    for k in range(nzermat):
        # Filling z_out with proper values
        j  = k-2+i0   # j is an index parsing the array z
        if( (j>=0) and (j<nzer) ):
            z_out[k] = z[j] * scalingFactor
        else:
            z_out[k] = 0.00

    # flipping wavefront, if required
    if( flip==True ):
        z_out = flipZernike(z_out, 2)

    # rotating wavefront by an angle <angledeg>
    if( angledeg!=0 ):
        z_out = rotateZernike(z_out, 2, angledeg)


    # transforming zernike coeffs into slopes by a matrix product
    slopes = np.dot(Z2S, z_out)

    # scaling factor to transform "nanometers rms" to "WFS pixels"
    RASC = 180*3600/3.14159265358979324          # number of arcsec in 1 radian
    unitfactor = 1e-9 * 2.0/ DiamVLT / (pixelSizeArcsec/RASC)
    slopes = slopes * unitfactor

    return slopes







def moveDerotatorToPosition( angle ):
    """

    :-)

    Super-dummy function, supposed to turn the de-rotator ...
    Of course this function needs to be replaced by something *actually* rotating the de-rotator

    """
    print "Moving de-rotator to position %g degrees\n" % float(angle)




def acquireSPARTAData( nframes ):
    """
    This function simulates the acquisition of <nframes> acquisitions of slopes
    made by SPARTA.
    Obviously, this has to be replaced by a real SPARTA acquisition, or by reading a FITS file
    of slopes from SPARTA.

    """
    NSLOPES = 136   # number of slopes returned by sparta ...
    # creating a dummy array, plenty of random numbers, just to test procedures ...
    return np.random.rand(nframes, NSLOPES)-0.5



def getCurrentReferenceSlopes( void ):
    """
    This function simulates the query of the currently applied reference slopes
    Obviously, this has to be replaced by a real SPARTA procedure
    """
    NSLOPES = 136   # number of slopes returned by sparta ...
    # creating a dummy array, filled with 0.00s, just to test procedures ...
    return np.zeros(NSLOPES)


def calibrateRefSlopesSimple(nframes):
    """
    Input parameters:
    <nframes> : scalar integer, number of frames to grab at each acquisition.
                Recommended default value (if ever required for whatever reason) : 200


    Output:
    A 1D array of floating point numbers, correponding to the reference slopes


    """

    # calls a function/procedure from sparta to get the cuently applied reference slopes
    curRefSlopes = getCurrentReferenceSlopes( )
    # acquire slopes data from SPARTA
    slopes = acquireSPARTAData( nframes )
    # time-average slopes, and store result in an array
    newRefSlopes =  np.average(slopes, axis=0) + curRefSlopes
    return newRefSlopes



def calibrateRefSlopesUsingDerotator(nframes, nsteps):
    """
    Input parameters:
    <nframes> : scalar integer, number of frames to grab at each acquisition.
                Recommended default value (if ever required for whatever reason) : 200
    <nsteps>  : number of position angles the rotator needs to do
                Recommended efault value (if ever required for whatever reason) : 12

    Output:
    A 2D array of floating point numbers, correponding to the reference slopes
    for every degree from 0 to 360.

    How it works:
    The procedure will grab <nframes> frames of slopes, for a given
    set of de-rotator angles. It will time-average them and add the current
    reference slopes back to it. Then, it will interpolate all these
    data between the measurement angles, for another set
    of <nsamp> angles spanning the range 0-360 degrees.

    What to do with the output:
    This interpolated data set will have to be saved somewhere in the system
    configuration as reference slopes.
    Then, while on-sky observing, for any angle A of the derotator, the
    reference slopes should be taken in the calibrated+interpolated dataset
    at the closest angle to A.

    """
    # number of slopes of the Shack-Hartmann
    NSLOPE = 136

    # allocate memory for result
    calibref = np.zeros((NSLOPE, nsteps))

    # Create a list of angles for the derotator
    # from 0 (included) to 360 (excluded !!) degrees.
    angle = np.arange(nsteps)*360.0/nsteps   # if nsteps=12, then angle=[0,30,60,90,120,150,180,210,240,270,300,330]


    # ................ MEASUREMENTS .......................
    #
    # calls a function/procedure from sparta to get the cuently applied reference slopes
    curRefSlopes = getCurrentReferenceSlopes( 0 )
    # loop over de-rotator positions
    for k in range(nsteps):
        # move de-rotator to position
        moveDerotatorToPosition, angle[k]

        # acquire slopes data from SPARTA
        slopes = acquireSPARTAData( nframes )
        # time-average slopes, and store result in an array
        calibref[:,k] = np.average(slopes, axis=0) + curRefSlopes



    # ................ INTERPOLATION .......................
    #
    # Now the array of slopes will be interpolated for any angle from 0
    # to 360 degree, on <nsamp> steps.
    nsamp       = 361   # this number could possibly change after some tests during AITs.
    nfourier    = nsteps//2   # integer division
    # allocate memory for the interpolated result
    calibinterp = np.zeros((NSLOPE, nsamp))
    # list of the angles where the interpolation will be calculated.
    # these angles range from 0 (included) to 360 (included) degrees.
    theta       = np.arange(nsamp)*360.0/(nsamp-1)
    # loop over slopes
    for i in range(NSLOPE):
        # loop on Fourier modes, from 0 to nfourier (included, i.e. nfourier+1 modes)
        for j in range(nfourier+1):     # warning: nfourier+1 iterations here !!
            cosinus = np.cos(angle * np.pi/180.0 * j) * 2.0 / nsteps
            sinus   = np.sin(angle * np.pi/180.0 * j) * 2.0 / nsteps
            acc = np.sum(calibref[i,:] * cosinus)
            ass = np.sum(calibref[i,:] * sinus)
            if (j==0) or (j==nfourier):
                acc /= 2.0
                ass /= 2.0

            calibinterp[i,:] += acc * np.cos(theta * np.pi/180.0 * j) + ass * np.sin(theta * np.pi/180.0 * j)


    # .................. CENTERING ...........................
    #
    # We now need to re-center all the reference slopes, in order to
    # set them at the center of the circle they've made when the
    # rotator turns. One has to subtract the tilt in X, and the tilt in
    # Y, and replace them by the average value in X and Y.
    #
    # First, computation of tilt in x and y for each theta, by
    # averaging all the slopes in X, and all slopes in Y.
    # Slopes 1 to NSLOPE/2=68 correspond to X, 69 to 136 correspond to Y.
    tx = np.average(calibinterp[0:NSLOPE/2,:], axis=0)
    ty = np.average(calibinterp[NSLOPE/2:NSLOPE,:], axis=0)
    avg_tx = np.average(tx)   # this is the average (over theta) of the x-tilt
    avg_ty = np.average(ty)   # this is the average (over theta) of the y-tilt
    tx -= avg_tx
    ty -= avg_ty
    # Then subtraction of the tilt for any angle
    for k in range(NSLOPE/2):
        j = k+NSLOPE/2   # index of slopes in Y (yorick numbering, start at 1)
        calibinterp[k,:] -= tx
        calibinterp[j,:] -= ty


    return calibinterp

