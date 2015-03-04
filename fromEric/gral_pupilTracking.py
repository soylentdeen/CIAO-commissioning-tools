# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 08:58:47 2014

@author: egendron
"""

#sys.path.insert(0, '/home/egendron/SVN_HRA/lib/python') #Add this python session lib path

import numpy as np
import pyfits
import matplotlib.pyplot as plt
import gral_zer
#from gral_gendron_utils import *



def plsh(v, n, sparta=1):
    """
    Special function for E. Gendron, useful for debugging
    and/or for data visualization.

    Might also be useful to others ... so i leave it here.

    """
    nsub = v.shape[0] / 2
    x = np.linspace(-1,1,n)
    x,y = np.meshgrid(x,x)
    r = np.sqrt(x*x+y*y)
    # defines outer and inner radiuses that will decide of validity of subapertures
    # inner radius is imposed.
    # outer one will be computed so that it will match the number of subapertures in v
    validint = 0.1
    rorder = np.sort(r.reshape(n*n))
    # number of subapertures not valid due to central obscuration
    ncentral = n*n - np.sum(r>validint)
    # determine value of external radius so that the test (validint < r < validext)
    # leads to the correct number of subapertures
    validext = rorder[ncentral + nsub]
    # get the indexes of valid subapertures in the 9x9 map
    valid = (r<validext) & (r>validint)
    ivalid = np.where(valid)
    # feeding data <v> into <vv>
    vx = np.zeros([n,n])
    vy = np.zeros([n,n])
    if( sparta==0 ):
        vy[ivalid] = v[0:nsub]
        vx[ivalid] = v[nsub:]
    else:
        vx[ivalid] = v[0::2]
        vy[ivalid] = v[1::2]
    plt.quiver(x, y, vx, vy, pivot='mid')






# Simulated acquisition of flux data from SPARTA
# PLEASE REPLACE THIS EITHER BY A REAL ACQUISITION FROM SPARTA,
# OR BY READING A FITS FILE ALREADY ACQUIRED
def gral_acquireSPARTAFluxData(nframes):
    """
    Input arguments:
    <nframes> : number of frames to be acquired by sparta

    Simulated acquisition of flux data from SPARTA
    PLEASE REPLACE THIS EITHER BY A REAL ACQUISITION FROM SPARTA,
    OR BY READING A FITS FILE ALREADY ACQUIRED
    """
    NSLOPE = 68
    data = np.random.rand(NSLOPE, nframes) * 22.0
    return data



# Simulated acquisition of slopes + volts data from SPARTA
def gral_acquireSPARTASlopeAndVoltData(nframes):
    """
    Input arguments:
    <nframes> : number of frames to be acquired by sparta
    Simulated acquisition of slopes + volts data from SPARTA
    """
    path = '/home/egendron/Projets/GRAVITY/SVN_Packages/'
    case = "gravity"
    if( case=="canary" ):
        slope = pyfits.getdata(path+'slopescanary_TS.fits')
        volt  = pyfits.getdata(path+'voltscanary_TS.fits')
        return slope, volt
    elif(case=="gravity"):
        # PLEASE REPLACE THIS EITHER BY A REAL ACQUISITION FROM SPARTA,
        # OR BY READING A FITS FILE ALREADY ACQUIRED
        #
        # FEEL FREE TO DO HERE WHAT YOU NEED TO BE ABLE TO RETURN RELEVANT
        # CLOSED-LOOP DATA WITH SLOPES+VOLTS (preferably with some turbulence
        # or with some noise : the DM has to MOVE in some way, so that the
        # algorithm can run properly )
        #
        a = pyfits.getdata(path+'LoopData_2.fits')
        slope = a.field(4)
        volt = a.field(5)
        return slope, volt





def pupilTrack_fluxMethod(nframes):
    """
    Input arguments:
    <nframes> : number of frames of flux data to be acquired by sparta
    Output:
    pupil position in X and Y, in arbitrary units to be calibrated during AITs.

    Procedure to track the pupil position using the algorithm based on flux data.
    This is the simplest and more robust algo, however some oral comments from
    some ESO people indicate this algo could fails because of differential motions
    between the MACAO DM and the VLT primary mirror. However, those doubts are not
    officially published by ESO.

    Anyway, the algo grabs some flux data from SPARTA, and analyses the flux
    differences between edge subapertures (outer and inner ones) to deduce
    the pupil position.

    """
    # acquire flux data from SPARTA
    flux = gral_acquireSPARTAFluxData(nframes)
    flux = np.average(flux, axis=1)

    # Number of subapertures across the pupil
    nssp = 9
    # generates a map of the distance of subapertures to pupil center
    # and map of the azimut angle of subapertures
    r,theta = gral_zer.mkxy(nssp,0)
    # defines outer and inner radiuses that will decide of validity of subapertures
    validext = 1.15
    validint = 0.1
    # get the indexes of valid subapertures in the 9x9 map
    valid = (r<validext) & (r>validint)
    # get indexes of subap at the outer limit of the pupil within a 9x9 map
    out = (r<validext) & (r>validext-2.0/nssp)
    # get indexes of subap at the inner limit of the pupil within a 9x9 map
    inner = (r>validint) & (r<validint+2.0/nssp)
    # now producing a 9x9 map with outer subapertures coded with 2, and
    # inner ones coded with 3
    carte = np.zeros_like(r)
    outcode = 2
    incode = 3
    carte[np.where(out)] = outcode
    carte[np.where(inner)] = incode

    # determining the indexes of outer and inner subapertures within
    # a list of 68 valid ones
    outlist = np.where(carte[valid]==outcode)
    inlist = np.where(carte[valid]==incode)

    # determining the azimut angles of outer and inner subapertures
    outangle = theta[valid][outlist]
    inangle = theta[valid][inlist]

    # Computing the total flux for outer and inner subaps
    outflux = np.mean(flux[outlist])
    influx = np.mean(flux[inlist])

    # define what is the average minimum flux per subaperture
    # this number should be adjusted during AITs
    minflux = 10.;

    # comparing average flux to the minimum limit allowed
    if( (outflux<minflux) or (influx<minflux)):
        print "The flux is too low %f %f" % (outflux, influx)
        tx = 0.00
        ty = 0.00
    else:
        # weighting the subaperture flux by cosine of their azimut angle, so
        # that the pupils on the right are +1, -1 on the left, 0 on the
        # top and bottom .. this will allow us to get the tilt in X
        tx = np.mean(flux[outlist] * np.cos(outangle))/outflux
        tx = tx - 2*np.mean(flux[inlist] * np.cos(inangle))/influx
        # same for tilt in Y, using sine instead of cosine
        ty = np.mean(flux[outlist] * np.sin(outangle))/outflux
        ty = ty - 2*np.mean(flux[inlist] * np.sin(inangle))/influx

    return tx, ty









def pupilTrack_slopesMethod(nframes, latency, Fs, norder):
    """
    Input args:
    <nframes> : number of frames to be acquired (2000 recommended)
    <latency> : time between the beginning of integration of the WFS and the
                start of the application of the command on the DM
    <Fs> : sampling frquency in Hertz
    <norder> : scalar int, 9.

    Procedure to track the pupil position using the algorithm based on slopes data.
    This algo is equivalent to the one coded by J. Kolb for the AOF.
    However, we might need to recode it to understand it more easily, as we
    know we will not get enough support from ESO.
    This algo is less robust than flux-based one, but has the advantage of
    retrieveing the position of the DM, instead of the position of the
    edge the edge of the pupil (VLT mirror M2).


    """
    # Grab nframes of slopes+volts from SPARTA
    slope, volt = gral_acquireSPARTASlopeAndVoltData(nframes)

    # express latecy in frames, and splits the integer part, and
    # the fractionnal remainder
    latency *= Fs;
    frac_latency = latency - int(latency)

    # computing the average voltages that were applied during WFS
    # integration periods
    volt  = (1-frac_latency)*volt[1:nframes-1,:] + frac_latency*volt[0:nframes-2,:]
    # truncating slope set to get only the slopes frames where we know the
    # applied voltages. The 2 firsts elements of the sequence have seen
    # voltages that are not in the volt set.
    slope = slope[2:nframes,:]
    # update <nframe> to its new value now, since datasets were truncated
    # in the 2 previous lines of code
    nframes -= 2

    # Compute (opposite of) voltage increments between 2 frames
    volt = -volt[1:nframes,:] + volt[0:nframes-1,:]
    # Compute the slopes increment
    slope = slope[1:nframes,:] - slope[0:nframes-1,:]

    # Least-square inverting the slopes matrix
    # result <s1> is transposed wrt <slope>
    s1 = np.linalg.pinv(slope, rcond=1e-2)

    # Now we compute the retrieved interaction matrix, as seen by the system
    # during the dataset. The shape of <mia> is (nb of slopes, nb of volts)
    mia = np.dot( s1, volt )

    # Read the Zernike-2-slope matrix
    # This could also be passed to the function as an argument
    path = '/home/egendron/Projets/GRAVITY/SVN_Packages/'
    Z2S = pyfits.getdata(path+'Z2S_136s_119z.fits').T
    Jmax = (norder*(norder+1))/2
    Z2S = Z2S[:,0:Jmax-1]

    # Reads *nominal* GRAVITY interaction matrix
    # The HO matrix is stored if FITS as NAXIS1=60, NAXIS2=136.
    # In python the HO matrix is read as (136,60), i.e. with dimensions
    # hopefully conformable with <mia>
    HOmi = pyfits.getdata(path+'HORecnCalibrat.RESULT_IM.fits')

    # Then the HOmi is differentiated along X and along Y
    # The 'derivative matrices' all have a shape conformable with the <mia>
    # and with the <HOmi>, i.e. (nb of slopes, nb of volts)=(136,60)
    dHOmi_dx, dHOmi_dy = derivative_MI(HOmi, Z2S)

    # finally the <mia> is searched as being the weighted sum of 3 matrices:
    # mia = A.(HOmi + Tx.dHOmi_dx + Ty.dHOmi_dy)
    # with A, Tx, Ty floating-point scalar values.
    # The constants Tx and Ty will be the pupil shifts.

    Tx, Ty = fittingProcedureOfInterMat(mia, HOmi, dHOmi_dx, dHOmi_dy)

    return Tx, Ty




def derivative_MI(HOmi, Z2S):
    """
    Input arguments:
    <HOmi> : nominal interaction matrix of the system, shape (136,60)
    <Z2S>  : matrix slopes-to-Zernike, shape (nzer, 136)

    The procedure will first reconstruct the influence functions of the DM
    as 'seen' by the HOmi, on the Zernike modes.
    Thanks to the article of R.J. Noll 1976, we know that any zernike mode
    can be differentiated, and the derivative can be expressed on zernike
    thanks to a matrix Gx (and Gx for derivative along y).
    Thanks to linearity properties, the derivative of the interaction matrix
    with respect to a translation dx of the influence functions, is equal
    to the interaction of the derivative of the influence functions.
    So, we will take the derivative of the influence functions, and compute
    the associated interaction matrices.
    Because we're using zernike modes, and because the derivative of
    zernike can be expressed on zernike themselves, we only need the Z2S to
    do all these operations.

    The function returns the 2 derivatives (1 in x, and 1 in y) of the HOmi.
    """
    # Inverts the Z2S it to get a slope-2-zernike matrix
    S2Z = np.linalg.pinv(Z2S)
    # get number of zernike and number of actuators
    nslopes,nzer = Z2S.shape
    # Compute matrices for derivating Zernike
    Gx = gral_zer.gammX(nzer, 2)[1:,:]  # remove piston
    Gy = gral_zer.gammY(nzer, 2)[1:,:]  # remove piston
    # Taking derivatives of the HOmi

    # Translate the HOmi measurements into zernike modes : the influence
    # functions will be expressed on Zernike modes (119,60)
    influFuncZer = np.dot(S2Z, HOmi)
    # Then zernike are differentiated along X and along Y and transformed
    # back to slopes
    dHOmi_dx = np.dot(Z2S, np.dot(Gx, influFuncZer))
    dHOmi_dy = np.dot(Z2S, np.dot(Gy, influFuncZer))

    return dHOmi_dx, dHOmi_dy






def fittingProcedureOfInterMat(mia, HOmi, dHOmi_dx, dHOmi_dy):
    """
    Input arguments:
    <mia> : an interaction matrix (IM), retrieved from closed-loop data
    <HOmi> : the nominal IM of the system
    <dHOmi_dx> : the derivative of the IM for a translation 'dx' of the DM
    <dHOmi_dy> : idem, in y.

    The procedure is looking for the linear combination of the 3 matrices
    A.(HOmi + Tx.dHOmi_dx + Ty.dHOmi_dy
    which is as close as possible to <mia>.
    The variables A, Tx, Ty are floating-point scalar values.
    Tx and Ty are the pupil shifts, and are returned by the procedure.

    """
    # allocate memory for a matrix 3x3..
    projMat = np.zeros([3,3])
    # .. and for a vector with 3 coeffs (i.e. triplet)
    projData = np.zeros(3)

    # The 3x3 matrix will be filled with the 'sum-product' of each couple
    # of the 3 matrices <HOmi>, <dHOmi_dx>, <dHOmi_dy>.
    # The 'sum-product' is the sum of all the products element-by-element
    # of 2 matrices.
    # As the 3x3 matrix is a symmetric one, we only
    # compute one-half, it's quicker.
    # The triplet will contain the 'sum-product' of the mia on the 3 matrices
    listMat = (HOmi, dHOmi_dx, dHOmi_dy)
    for i in range(3):
        for j in range(i,3):
            projMat[i,j] = np.sum(listMat[i]*listMat[j])
            projMat[j,i] = projMat[i,j]
        projData[i] = np.sum(listMat[i]*mia)

    # inverting <projMat>
    antiProjMat = np.linalg.pinv(projMat)

    # multiply the coeffs of the projection of the mia.
    # The fitted coeffs are [A,A*Tx,A*Ty] : we will extract Tx and Ty.
    fittedCoeffs = np.dot(antiProjMat, projData)
    A = fittedCoeffs[0]
    Tx = fittedCoeffs[1] / A
    Ty = fittedCoeffs[2] / A

    return (Tx, Ty)


