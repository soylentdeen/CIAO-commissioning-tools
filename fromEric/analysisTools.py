# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 13:46:54 2015

@author: egendron
"""

import pyfits
import numpy as np
#import gral_createModalBasis

import matplotlib.pyplot as plt



def checkSMAandCo(path, filename):
    """
    Input args:
    <path> directory where files are
    <filename> filename of the slopes+voltage data set

    """

    # reading data
    a = pyfits.getdata(path+'/'+'FtWorth_2.fits')
    vHO = a.field(5).T
    vTT = a.field(6).T
    slopes = a.field(4).T

    # reading other data
    AWFbasis = pyfits.getdata(path+'/'+'HOCtr.AWF_IM_KERNEL.fits')
    SMAbasis = pyfits.getdata(path+'/'+'HOCtr.SMA_BASIS.fits')
    HO2TT = pyfits.getdata(path+'/'+'HOCtr.HO_TO_TT.fits')
    TT2HO = pyfits.getdata(path+'/'+'HOCtr.TT_TO_HO.fits')
    pistonMode = pyfits.getdata(path+'/'+'HOCtr.PRA_PISTON_MODE.fits')
    pistonProj = pyfits.getdata(path+'/'+'HOCtr.PRA_PISTON_PROJECTION.fits')
    S2M = pyfits.getdata(path+'/'+'S2M.fits')
    M2V = pyfits.getdata(path+'/'+'M2V.fits')

    # projecting voltages onto piston mode
    pistonProj /= np.dot(pistonProj, pistonMode)  # normalizing things (just as sparta does internally)
    piston = np.dot(pistonProj, vHO) * np.sqrt(np.sum(pistonMode**2))  # projecting onto piston mode
    print 'This is PISTON component : %f ± rms %f ' % (np.average(piston) , np.std(piston))

    # transforming TT signals into HO .. not sure this is useful
    # but on the present data set, one cannot check
    vHO += np.dot(TT2HO, vTT)

    # computing voltages increments between 2 frames
    dvHO = vHO[:,1:] - vHO[:,:-1]
    dvTT = vTT[:,1:] - vTT[:,:-1]  # useless
    slopes = slopes[:,1:]          # drop first frame

    # cutting data set to avoid the beginning, which is openloop
    Nstart = 1111
    dvHO = dvHO[:,Nstart:]
    dvTT = dvTT[:,Nstart:]
    slopes = slopes[:,Nstart:]


    # Here i'm trying to re-generate the command increments that
    # occured during the closed loop, and see how well they match the
    # command increments computed using the SPARTA voltages.
    #
    MC = gral_createModalBasis.computeSystemControlMatrix( S2M, M2V, TT2HO*0 )
    sm = np.dot(MC, slopes)

    # now i'm doing a linear regression, on all actuators, to
    # find what was the close-loop gain that was applied
    na     = 60
    gain   = np.zeros(na)
    offset = np.zeros(na)
    for i in range(na):
        gain[i], offset[i] = np.polyfit(sm[i,:], -dvHO[i,:] ,1)
    avgain = np.average(gain)
    print 'Average gain found is %f ± rms %f' % (avgain, np.std(gain))
    print 'Offset variation is rms %f' % np.std(offset)

    # Subtraction between SPARTA incr-voltage conmmands, and
    # our calculation of what should be the incr-volt commands
    errcomm = sm[0:na,:] * avgain + dvHO
    print 'Error between sparta volts and computed volts : %f' % np.average(np.std(errcomm,axis=1))

    # Using the M2V, we're gonna re-generate the projection matrix onto
    # the non-commanded modes
    na = 60
    nuseful = M2V.shape[1]
    M2VHO = np.zeros((na, nuseful))
    M2VHO[:,0:2] = TT2HO
    M2VHO[:,2: ] = M2V[:na,2: ]
    projAWF_parall = np.dot(M2VHO, gral_createModalBasis.pseudoinv( M2VHO, 0 ))
    projAWF_kernel = np.identity(na) - projAWF_parall
    # using sparta 'V2 = AWF_IM_KERNEL' to check we're retrieving the
    # same projection matrix. We use the <AWFbasis> passed to sparta, and use the
    # internal sparta algorithm (the famous V2.transpose(V2)) for generating it.
    # If the <AWFbasis> is correct, both computation should give the same
    # result.
    projAWF = np.dot(AWFbasis, AWFbasis.T)
    # computing difference between sparta proj matrix and our proj matrix
    errAWF = np.sum((projAWF-projAWF_kernel)**2) / np.sum((projAWF)**2)
    print 'Error on AWF matrix is %f' % errAWF

    # now computing the fraction of the command which
    # belongs to the commanded space
    dpHO = np.dot(projAWF_parall, dvHO)
    # and the command lying out of the commanded space (belong to kernel)
    dkHO = np.dot(projAWF_kernel, dvHO)
    # see where kernel-commands were sent
    kpow = np.sqrt(np.sum(dkHO**2, axis=0)) > 5e-7
    plt.plot( kpow )
    print 'Nb of occurences SMA sent kernel-commands: %d' % np.sum(kpow)


