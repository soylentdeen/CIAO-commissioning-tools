# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 13:49:59 2015

@author: egendron

"""


#
import pyfits

import numpy as np
import matplotlib.pyplot as plt





def pseudoinv( mat, nfilt, calm=False ):
    """
    Input args:
    <mat> : a matrix (2D of floats)
    <nfilt> : number of modes to be filtered in the inversion
    <calm> : optionnal boolean

    The procedure will compute the pseudo-inverse of <mat> using the
    formula:
    D+  =  (Dt * D)**-1  *  Dt
    which is the solution for v=D+.m that minimizes ||m-Dv||**2

    The inversion of the matrix (Dt * D) is obtained by diagonalisation,
    and inversion of the eigenvalues (except the <nfilt> last ones, set to
    0).


    When the variable <calm>==True, the inversion is done using
    D+  =  (Dt * D + diag)**-1  *  Dt
    where diag is a diagonal matrix with non-null values on diagonal
    elements in the range 40:60, with a value of 5% of the average value of
    the diagonal of (Dt*D). This is equivalent to minimize the quantity
    ||m-Dv||**2 + ||B.v||**2 and the solution is
    D+  =  (Dt * D + Bt * B)**-1  *  Dt

    """

    DtD = np.dot(mat.T, mat)
    if calm :
        avdiag = np.average( np.diag(DtD) )
        for i in range(40,60):
            DtD[i,i] += avdiag*0.05
    lam, u = diagonalisation(DtD)
    for i in range(lam.shape[0]-nfilt):
        u[:,i] /= np.sqrt(lam[i])
    if( nfilt>0 ):
        u[:,-nfilt:] = 0
    return np.dot(np.dot(u, u.T), mat.T)








def diagonalisation(A):
    """
    Input arg:
    <A> : square, symmetric matrix

    Ouput:
    tuple with sorted eigenvalues, and corresponding eigenvectors

    This function is equivalent to np.linalg.eigh() (=diagonalisation
    of a symmetric matrix), plus sorting the eigenvalues.
    """
    # diagonalization of matrix A
    # One will find matrix M so that :
    #      M.T * A * M = diagonal
    (eigenvalue, M) = np.linalg.eigh(A)
    # Now we need to sort the eigenvalues and the eigenvectors, it not already done...
    sortindex = np.argsort(eigenvalue)[::-1]
    eigenvalue = eigenvalue[sortindex]
    M = M[:,sortindex]
    return (eigenvalue, M)










def filterOutPiston( modes, pistonMode, pistonProj ):
    """
    Input args :
    <modes>      2D array of floats, is an array containing N modes, it is an array
60xN
    <pistonMode> 1D array of floats, is the piston mode defined on actuator voltages
(60 components)
    <delta>      2D array of floats, is the geometric cov matrix of the DM (60x60),
fixed.

    The function will suppress piston from all the modes of matrix <modes>.
    """

    # compute the scalar product between pistonMode and (delta.pistonMode)
    # this is a projection of the piston mode on itself
    pnorm = np.dot(pistonMode, pistonProj)
    # compute the scalar products between each of the modes and (delta.pistonMode)
    # this is a projection of each mode on the piston
    proj  = np.dot(modes.T, pistonProj)
    # normalize by the value of the projection of piston on itself
    proj /= pnorm
    # now we know the amount of piston in each mode, we need to subtract it ..
    (m, n) = modes.shape
    fmodes = np.zeros((m,n)) # create the output, filtered modes
    for i in range(n):    # for each mode ....
        fmodes[:,i] = modes[:,i] - pistonMode*proj[i]
    return fmodes






def KLbasisSlopesSpace(fullNoll, Z2S, miaHO, miaHO1):
    """
    Input args:
    <fullNoll> : matrix (120, 120), given as an input (covariance
                 matrix of zernike modes)
    <Z2S> : matrix zernike to slopes
    <miaHO> : HODM interaction matrix
    <miaHO1> : filtered pseudo inverse of miaHO

    Outputs:
    The procedure returns
    - the eigenvalues of the KL modes of the
    measurements projected onto the commanded subspace of miaHO1
    - the modes
    """
    n = fullNoll.shape[0]   # number of zernike in the noll matrix
    ns, nz = Z2S.shape      # nber of slopes and zernike
    # see what's the common number of zernike of both
    nz = min(n,nz)
    # cut matrices to the right size
    fullNoll = fullNoll[0:nz, 0:nz]
    Z2S  = Z2S[:, 0:nz]
    # compute covariance matrix of slopes Css = Z2S * Czz * Z2S.T
    covSlopes = np.dot( Z2S, np.dot(fullNoll, Z2S.T))

    proj = np.dot(miaHO, miaHO1)
    covSlopes = np.dot( proj, np.dot(covSlopes, proj.T))

    turb, eigss = diagonalisation(covSlopes)
    return turb, eigss








def createModalBasis( miaHO, miaTT, delta, fullNoll, Z2S, nfilt ):
    """

    Input arguments :
    <miaHO>   is the measured interaction matrix of DM (120x60, i.e. 120 lines, 60
columns)
    <miaTT>   is the measured interaction matrix of TT (120x2, i.e. 120 lines, 2
columns)
    <delta>   is the geometric covariance matrix of actuators, it is computed
elsewhere.
               It is a square, symmetric matrix 60x60
    <fullNoll> 2D array, covariance matrix of Zernike modes in atmospheric turbulence,
               computed from Noll (1976) article.
    <Z2S> 2D array, matrix zernike to slopes
    <nfilt>   scalar int, to be adjusted during AITs. Equal to 9 by default.

    Outputs arguments :
    <TT2HO>      : matrix of the tilts of the DM (will be used in sparta)
    <HO2TT>      : projection matrix from HODM voltages onto the 2 TT voltages
    <M2V>        : matrix of system modes (will be used in the modal control optim
and elsewhere)
    <pistonMode> : piston mode (the mode that makes the DM flat, but pistonned)
    <pistonProj> : the vector that extracts the value of the piston from a set of
voltages of the HODM.
    <S2M>        : modal control matrix, that will be used in particular in the modal
                    control optimization.
    <AWFbasis>   : basis of modes for anti windup filter
    <SMAbasis>   : basis of modes for saturation management

     note: see ESO doc
      VLT-SPE-ESO-16100-4856_iss1_SPARTA-Light Reference Implementation Technical
Specifications

    """
    # get number of slopes <ns> and number of actuators of DM <na>
    (ns, na) = miaHO.shape

    # generalized inverse of miaHO
    # nfilt = 5      # to be adjusted during AITs
    miaHO1 = pseudoinv( miaHO, nfilt)     # miaHO1 has transposed shape wrt miaHO

    # Computation of set of DM voltages that reproduce a pure tilt
    # WARNING: those tilts may contain some piston
    TT2HO = np.dot(miaHO1, miaTT)

    # Computation of PISTON mode as the least visible mode with highest variance
    # method 1: obsolete
    # lam, eigenvect = simultaneousDiagonalisation( delta[0:60,0:60],np.dot(miaHO.T, miaHO) )
    # pistonMode = eigenvect[:,-1]
    # method 2 (preferred)
    lam, eigenvect = diagonalisation( np.dot(miaHO.T, miaHO) )
    pistonMode = eigenvect[:,-1]
    pistonProj = np.dot(delta,pistonMode)

    # Removing pistonMode from TT2HO
    TT2HO = filterOutPiston( TT2HO, pistonMode, pistonProj )

    # computing the reverse projection HO2TT
    HO2TT = np.dot(TT2HO, pseudoinv(TT2HO,0))

    # computing modes as KL in the measurement space
    turb, eigss = KLbasisSlopesSpace(fullNoll, Z2S, miaHO, miaHO1)
    nuseful = na-nfilt
    modes = np.dot(miaHO1, eigss )[:, 2:nuseful]

    # make all those modes orthogonal to piston
    modes = filterOutPiston( modes, pistonMode, pistonProj)

    # Normalize the basis to the same value as TT2HO.
    var = np.sum(modes*np.dot(delta, modes), axis=0)
    varTT2HO = np.average( np.sum(TT2HO*np.dot(delta, TT2HO), axis=0) )
    modes = np.dot( modes, np.diag(np.sqrt(varTT2HO/var)) )

    M2V = np.zeros((na+2, nuseful))
    M2V[:na,2: ] = modes
    M2V[na+0, 0] = 1.
    M2V[na+1, 1] = 1.

    # Gather the 2 interaction matrices miaHO and miaTT into a single matrix mia (tilts at the end)
    mia = np.zeros((ns, na+2))
    mia[:,0:na]    = miaHO
    mia[:,na:na+2] = miaTT

    # Now computing command matrix
    M2S = np.dot(mia, M2V)      # this is the Modal Interaction Matrix
    S2M = pseudoinv( M2S, 0 )   # this is the Modal Control Matrix
    #   ................. AWF ......................
    # compute the projection basis of a voltage vector onto the commanded space
    #
    M2VHO = np.zeros((na, nuseful))
    M2VHO[:,0:2] = TT2HO
    M2VHO[:,2: ] = M2V[:na,2: ]
    projAWF = np.identity(na) - np.dot(M2VHO, pseudoinv( M2VHO, 0 ))
    lam, mod = diagonalisation(projAWF)
    AWFbasis = np.dot(mod[:,0:nfilt], np.diag(1/np.sqrt(lam[0:nfilt])))

    # some debug stuff (useless in the final version)
    debugmode = False
    if debugmode==True:
        plt.figure(1)
        plt.clf()
        plt.subplot(3,1,1)
        mc = computeSystemControlMatrix( S2M, M2V, TT2HO )
        plt.imshow(mc, interpolation='none')
        covSlopes = np.dot( Z2S, np.dot(fullNoll, Z2S.T))
        covact = np.dot( mc, np.dot(covSlopes, mc.T))
        plt.subplot(3,1,2)
        plt.plot(np.sqrt(np.diag(covact)))
        covnoise = np.dot( mc, mc.T)
        plt.plot(np.sqrt(np.diag(covnoise)))
        plt.plot(np.sqrt(np.diag(covact)/np.diag(covnoise)))
        iii = np.identity(136)
        turbfull, eigss = KLbasisSlopesSpace(fullNoll, Z2S, iii, iii)
        plt.subplot(3,1,3)
        vartot = np.sum(turbfull)
        plt.plot(100 - np.cumsum(turb)/vartot*100)
        print 'Residu = %f' % (1-np.sum(turb)/vartot)

    #   ................. SMA ......................
    SMAbasis = createSMAbasis(delta, pistonMode, pistonProj)


    return TT2HO, HO2TT, M2V, pistonMode, pistonProj, S2M, AWFbasis, SMAbasis







def createSMAbasis(delta, pistonMode, pistonProj):
    """
    Input args:
    <delta>   is the geometric covariance matrix of actuators, it is computed
elsewhere.
               It is a square, symmetric matrix 60x60
    <pistonMode> : piston mode (will be used in sparta)

    This will create a basis orthogonal to piston
    with last modes having large voltages and only small phase variance.

    """
    m = filterOutPiston( np.identity(60), pistonMode, pistonProj )
    lam, mo = diagonalisation( np.dot(m.T, np.dot(delta,m)) )
    mo = np.dot(m, mo)

    SMAbasis = np.zeros(delta.shape)
    SMAbasis[:,0] = pistonMode
    SMAbasis[:,1:] = mo[:,:-1]

    return SMAbasis









def pseudoOpenLoopReconstruction(slopesdata, voltsdataHO, voltsdataTT, miaHO, miaTT,
Fs, latency):
    """

    Input arguments :
    <slopesdata>  : set of N frames of slopes (i.e. centroids), (120xN)
    <voltsdataHO> : set of N frames of volts (i.e. DM commands) synchronized with
slopes, (60xN)
    <voltsdataTT> : set of N frames of volts (i.e. TT commands) synchronized with
slopes, (2xN)
    <miaHO>       : measured interaction matrix (120x60, i.e. 120 lines, 60 columns)
    <miaTT>       : measured interaction matrix (120x2, i.e. 120 lines, 2 columns)
    <Fs>          : sampling frequency (Hertz)
    <latency>     : latency, i.e. time between the start of WFS integration and time
when the
                  actuator reaches 50% of its command, in seconds.

    Returned value :
    open loop slopes data, as if they had been acquired in open-loop.

    This procedure applies both in open-loop or closed-loop.
    In the particular case of open-loop, all voltages voltsdata are 0, and the
function will
    just return <slopesdata>.

    The method to reconstruct pseudo openloop slopes is to subtract the
    action of the DM (projected into measurement space (slopes)) from
    the residuals. The DM action has to be subtracted with the proper
    time shift according to the system latency. This time-shift
    operation algorithm is evaluated in the time domain in our function
    (one usually make use of Fourier domain in other algorithms cited
    in the literature, but we do not make this choice here because we
    think it is not convenient).

    """
    # get number of slopes <ns> and number of actuators of DM <na>
    (ns, na) = miaHO.shape

    # slopes that were "added" by the action of the DM
    sv = np.dot(miaHO, voltsdataHO.T)
    sv = sv + np.dot(miaTT, voltsdataTT.T)

    # those slopes <sv> need now to be subtracted from the residual slopes
    # with the proper delay.
    # But first one needs to process latency a bit.

    # latency expressed in frames (i.e. about 1.35 frames for gravity).
    latency_frame = Fs * latency
    # latency needs to be split into integer and fractional part,
    # i.e. 1.35 = 1 + 0.35
    i_latency_frame = int(latency_frame)               # integer part, i.e. 1
    f_latency_frame = latency_frame - i_latency_frame  # franctional part, i.e. 0.35

    # now <sv> will be shifted in time
    (nrec, ns) = sv.shape
    svs = np.zeros((nrec, ns))  # just to allocate memory space same size as <sv>
    for i in range(nrec):
        j = i-i_latency_frame   # 1 frame before
        k = j-1                 # 2 frames before ...
        if(j<0): j=0            # except we don't have data before start of file !..
        if(k<0): k=0            # idem
        svs[:,i] = sv[:,j] * (1-f_latency_frame) + sv[:,k] * f_latency_frame

    # subtraction of time-shifted DM action from residual slopes data
    openloopSlopes = slopesdata.T - svs

    return openloopSlopes










def hcor(freq, Fs, latency, G, BP):
    """
    Input arguments :
    <freq> is a 1D array of frequencies (usually 1024 points ranging from Fs/2048 to
Fs/2).
    <Fs> is the sampling frequency
    <latency> is the latency, in seconds, between the BEGINNING of the integration
and the
              start of the command.
    <G> is a scalar, it's the loop gain
    <BP> is a scalar, the cutoff frequency of the DM (seen as a 1st order filter)

    On output, retunrs the square modulus of the correction transfer function of the
    system
    """
    Te=1./Fs                         # sampling period
    p = 1j*2*np.pi*freq              # Laplace variable
    Hint = 1./(1-np.exp(-p*Te))      # numeric integrator
    Hccd = (1.-np.exp(-p*Te))/(p*Te) # zero-order hold with 1/2 frame delay
    Hdac = Hccd                      # well, same.
    tdelay = latency - Te            # time between END of the integratino and start
of command
    Hret = np.exp(-p*tdelay)         # latency transfer function
    # transfer func of the DM, as a 1st order filter
    Hmir = 1./(1. + 1j*freq/BP)
    # open-loop transfer function
    Hbo = Hint * Hccd * Hdac * Hret * Hmir
    # correction transfer function
    Hcor   = 1./np.abs(1+Hbo*G)**2
    return Hcor







def modalControlOptimization( slopesdata, voltsdataHO, voltsdataTT, miaHO, miaTT,
S2M, M2V, gmax, Fs, latency, BP ):
    """
    Input arguments :
    <slopesdata> : set of N frames of slopes (i.e. centroids), (120xN)
    <voltsdataHO> : set of N frames of volts (i.e. DM commands) synchronized with
slopes, (60xN)
    <voltsdataTT> : set of N frames of volts (i.e. TT commands) synchronized with
slopes, (2xN)
    <miaHO>      : measured interaction matrix (120x60, i.e. 120 lines, 60 columns)
    <miaTT>      : measured interaction matrix (120x2, i.e. 120 lines, 2 columns)
    <S2M>        : modal control matrix
    <M2V>        : matrix of system modes
    <gmax>       : scalar floating point value, maximum gain.
    <Fs>         : sampling frequency (Hertz)
    <latency>    : latency, in seconds, i.e. time between the start of WFS
integration and time
                  when the actuator reaches 50% of its command.
    <BP>         : scalar, cutoff frequency of the DM (seen as a 1st order filter)

    Returned value :
    Array of optimal gains on the 50 modes.

    """
    # determines how many frames are recorded in the data set,
    # will be used later (could also be passed as an argument..)
    (nrec, ns) = slopesdata.shape

    # Reconstruction of pseudo open-loop slopes
    slopes = pseudoOpenLoopReconstruction(slopesdata, voltsdataHO, voltsdataTT,
miaHO, miaTT, Fs, latency)

    # conversion of slopes into modal coefficients
    modes = np.dot(S2M, slopes)
    (nmod, na) = modes.shape   # just to have the number of modes here

    # Produces an array of gains
    ngain = 15       # value to be adjusted during AITs
    gmin = 0.0       # TBC, maybe we'll keep gmin=0.001 to allow for static
aberration compensation, if any. TBC during AIT.
    G = np.linspace(np.sqrt(gmin), np.sqrt(gmax), ngain)**2   # 1D array from gmin
to gmax in ngain points

    npfft = nrec//2    # number of useful points in the FFT
    # create a 1D array of frequency ranging from Fs/nrec to Fs/2.0 in npfft points
    freq = np.linspace(Fs/nrec, Fs/2.0, npfft)
    optimumGain = np.zeros(nmod)  # memory alloc for 1D array of optimum gains

    # Initialisation of all transter functions computed for all the gains
    Hcor = np.zeros((ngain, npfft))
    for j in range(ngain):
        Hcor[j,:] = hcor(freq, Fs, latency, G[j], BP)

    # Fourier transform of modes coefficients.
    # Compute the square modulus of Fourier transform along time, for each mode.
    # Then multiplies by transfer function for different gains, sum everything
    # over frequency to get the error, and selects the minimum value.
    for i in range(nmod):
        # square modulus of Fast Fourier Transform
        tmp = np.abs(np.fft.fft( modes[i,:] ))**2
        # take only half of the points, and reject 1st point
        fftmodes = (2.0/nrec) * tmp[1:npfft+1]

        for j in range(ngain):
            phaseError = np.sum(Hcor[j,:] * fftmodes)

            # Initializes things at startup
            if(j==0):
                minPhaseError = phaseError
                jmin = j

            # Search for the minimum value, and the associated gain
            if( phaseError<minPhaseError ):
                minPhaseError = phaseError
                jmin = j


        optimumGain[i] = G[jmin]

    return optimumGain






def computeSystemControlMatrix( S2M, M2V, TT2HO ):
    """

    Input parameters :
    <S2M>      is the modal control matrix as computed by function
computeModalControlMatrix()
    <M2V>      is the matrix of the modes computed by the function createModalBasis()
    <TT2HO>    is the matrix of the DM-tilt modes (60,2) computed by the same function
               than <M2V>.

    Returned parameter : the command matrix of the system.
    """
    # matrix multiply  modes * S2M
    # In a "normal system",
    mc = np.dot(M2V, S2M)

    # transforming the tip-tilt correction matrix into a DM-tilt correction matrix
    tiltHO = np.dot(TT2HO, mc[60:62,:])
    # adding this DM-tilt correction matrix to the DM one
    mc[0:60,:] += tiltHO

    return mc







def computeSystemOptimizedControlMatrix( S2M, gainvector, M2V, TT2HO ):
    """

    Input parameters :
    <S2M>        is the modal control matrix as computed by function
computeModalControlMatrix()
    <gainvector> a 1D vector of modal gains, supposed to be provided by function
modalControlOptimization()
    <M2V>        is the matrix of the modes computed by the function createModalBasis()
    <TT2HO>      is the matrix of the DM-tilt modes (60,2) computed by the same
function
                 than <M2V>.

    Returned parameter : the command matrix of the system.
    """
    # let's allocate memory space
    tmp = np.zeros( S2M.shape )
    # let's multiply by the gains (dimensions of gainvector
    # and S2M must comply)
    ngain = gainvector.shape[0]
    for i in range(ngain):
        tmp[i,:] = S2M[i,:] * gainvector[i]

    # matrix multiply  M2V * (gains*S2M)
    mc = np.dot(M2V, tmp)

    # Now we need to take into account the tiptilt handling operated by
    # the CTR component of SpartaLight. We need to add to the DM part
    # of the matrix the compensation of tiptilt by the DM actuators.
    # transformation of the tip-tilt correction matrix into a DM-tilt correction matrix
    tiltHO = np.dot(TT2HO, mc[60:62,:])
    # adding this DM-tilt correction matrix to the DM one
    mc[0:60,:] += tiltHO

    return mc






"""


 _____         _
|_   _|__  ___| |_ ___
  | |/ _ \/ __| __/ __|
  | |  __/\__ \ |_\__ \
  |_|\___||___/\__|___/



 ____       _
|  _ \  ___| |__  _   _  __ _
| | | |/ _ \ '_ \| | | |/ _` |
| |_| |  __/ |_) | |_| | (_| |
|____/ \___|_.__/ \__,_|\__, |
                        |___/


"""




def plm(fi, mat, window, transpose=False):
    m = np.dot(fi, mat)
    n = mat.shape[1]
    k = 1 + int(np.sqrt(n-1))
    plt.figure(window)
    plt.clf()
    for i in range(n):
        plt.subplot(k,k,i+1)
        if( transpose ):
            plt.imshow(m[:,:,i].T)
        else:
            plt.imshow(m[:,:,i])







def modesFromKL(miaHO, miaHO1, turb, eigss, order=False):
    """
    A plotting routine for gendron ...
    Just ignore it.
    """
    modes = np.dot(miaHO1, eigss )
    realmes = np.dot(miaHO, modes )
    residu = np.sum((realmes - eigss)**2, axis=0)

    Watt = np.sum(modes[40:,:]**2,axis=0)
    merit = turb/Watt**2/np.exp(residu/0.2)
    if order==True :
        sortindex = np.argsort(residu)#[::-1]
        residu = residu[sortindex]
        turb   = turb[sortindex]
        Watt   = Watt[sortindex]
        merit  = merit[sortindex]
        modes  = modes[:,sortindex]
    plt.figure(1)
    plt.subplot(2,2,1)
    plt.plot(residu)
    plt.subplot(2,2,2)
    plt.plot( Watt )
    plt.subplot(2,2,3)
    plt.plot( np.sum(modes,axis=0))
    plt.subplot(2,2,4)
    plt.plot( turb )
    plt.yscale('log')

    return modes





def pseudoinv_calm( mat, nfilt, calm=1 ):
    """
    le gros tests du rico .....
    """

    DtD = np.dot(mat.T, mat)
    A = np.identity(60)
    A[40:60,40:60] *= calm
    lam, u = simultaneousDiagonalisation(A, DtD)
    u = u[:,0:60-nfilt]
    Du = np.dot(mat, u)
    return np.dot(u, pseudoinv( Du, 0))




def simultaneousDiagonalisation(A, B):
    """
    Input args:
    <A> a square matrix, symmetric, positive definite.
    <B> a square matrix, symmetric, positive definite.

    The procedure will diagonalize A and B, and find a matrix P
    in a way that:
    Pt.A.P = Identity
    Pt.B.P = diagonal matrix

    Usage:
    This procedure will be called by setting
    A = Delta (the geometric covariance matrix)
    B = D.T*D (with D the interaction matrix)
    """
    # diagonalization of matrix A
    # One will find matrix M so that :
    #      M.T * A * M = diagonal
    (eigenvalue, M) = diagonalisation(A)
    # Now the eigenvectors (columns of M) will be divided by the square root
    # of the eigenvalues, i.e.
    #      M  =  M / diag(sqrt(eigenvalues))
    for i in range(eigenvalue.shape[0]):
        M[:,i] /= np.sqrt(eigenvalue[i])
    # Then one compute the matrix
    #      MtBM  =  M.T * B * M
    MtBM = np.dot(M.T, np.dot(B,M))
    # and we diagonalize this last matrix to find a matrix N so that
    #      N.T * MtBM * N = diagonal
    (eigenvalue, N) = diagonalisation(MtBM)
    # finally, the modes are  M*N
    P = np.dot(M,N)
    # Proof:
    # If we compute the 2 following quantities :
    #     P.T * A * P = N.T * M.T * A * M * N
    #                 = N.T * 1/sqrt(eig) * M'.T * A * M' * 1/sqrt(eig) * N
    #                 = N.T * N
    #                 = identity
    #     P.T * B * P = N.T * M.T * B * M * N
    #                 = N.T * MtBM * N
    #                 = eigenvalues
    # we see that the basis <P> turns both matrices into a diagonal form.
    #

    return (eigenvalue, P)





"""


These comments are the list of commands you should run in a PYTHON
terminal to be able to use these functions.





# GENDRON DEBUG INIT !!!!
path = '/home/egendron/HOMEMAC/Projets/GRAVITY/PistonAntiWindup/'
v2z = pyfits.getdata(path+'IF_DM4/v2z_DM4.fits')
#z2v = pyfits.getdata(path+'IF_DM4/z2v_DM4.fits')
#screen = pyfits.getdata(path+'Simuls/screen2048.fits')
#sBuf = pyfits.getdata(path+'Simuls/sBuf2048.fits')
fi = pyfits.getdata(path + 'IF_DM4/if_macao_64x64_microns_400V.fits')
# adding the piston to the influence functions ...
pupmsk = (fi[0,:,:]!=0)
for i in range(60):
    fi[i,:,:] += pupmsk * v2z[0,i]/3.50

delta = np.dot(fi.reshape((60,64*64)),fi.reshape((60,64*64)).T)

"""

# Reads useful basic data : an interaction matrix (TT and HO) and a matrix
# <delta> are stored somewhere on the disk, pre-computed during
# AITs.
# The HO matrix is stored in FITS as NAXIS1=60, NAXIS2=136.
# In python the HO matrix is read as (136,60)
mypath = '/home/egendron/Projets/GRAVITY/SVN_Packages/'
miaHO = pyfits.getdata(mypath+'HORecnCalibrat.RESULT_IM.fits')
miaTT = pyfits.getdata(mypath+'TTRecnCalibrat.RESULT.IM.fits')
delta = pyfits.getdata(mypath+'delta_MACAO.fits')

#miaHO = pyfits.getdata(mypath+'simu_mia_HO.fits').T
#miaTT = pyfits.getdata(mypath+'simu_mia_TT.fits').T
#delta1 = pyfits.getdata(mypath+'simu_delta.fits')[0:60,0:60]



Z2S = pyfits.getdata(mypath+'Z2S_136s_119z.fits').T
#Z2S = pyfits.getdata(mypath+'Z2S_120s_119z_xxxyyy.fits').T
#norder=11
#Jmax = ((norder+1)*(norder+2))/2  # Jmax = 55
#S2Z = np.linalg.pinv(Z2S[:,0:Jmax-1])       # Inverts the Z2S it to get a slope-2-zernike matrix


#r,th = gral_zer.mkxy(64,0)
#cubzer = np.zeros((64,64,Jmax-1))
#pup = r<1.001
#for i in range(Jmax-1):
#    cubzer[:,:,i] = gral_zer.zer(r,th,i+2) * pup


fullNoll = pyfits.getdata(mypath+'ZernikeCovar_Noll.fits')
#noll = fullNoll[0:Jmax-1,0:Jmax-1]

# Creates modal basis and outputs <TT2HO> (matrix of the tilts of the
# DM (will be used in sparta)), <M2V> (matrix of system modes
# (will be used in the modal control optim and elsewhere)), <pistonMode>
# (piston mode (will be used in sparta)), and also <S2M> the modal
# control matrix, i.e. a control matrix that multiplies to slopes
# vectors and produces the modes coefficients.
(TT2HO, HO2TT, M2V, pistonMode, pistonProj, S2M, AWFbasis, SMAbasis) =
createModalBasis( miaHO, miaTT, delta, fullNoll, Z2S, 9)


mc = computeSystemControlMatrix( S2M, M2V, TT2HO )


"""

# Reads (or acquires..) a set of slopes and voltages
#s = pyfits.getdata(mypath+'simu_slopeset.fits')
s = np.random.randn(2048, 136)
vHO = s[:,0:60] * 0.0   # volts are just 0 here
vTT = s[:,0:2] * 0.0   # volts are just 0 here
(nframes, ns) = s.shape
s = s + np.random.randn( nframes, ns )*40   # add random numbers on slopes to
simulate noise


Fs = 200.0        # Sampling frequency (Hz), please replace by yours.
latency = 0.008   # Latency (s), please replace by yours. If you don't know, keep
0.008 s.
BP = 700.         # MACAO DM first resonnant frequency (Hz), should be ok.


# If you're able to run the modal control optimization, do this :
#
# Perform modal control optimization: optiGain is the array of the modal gains.
optiGain = modalControlOptimization( s, vHO, vTT, miaHO, miaTT, S2M, M2V, 0.8, Fs,
latency, BP )


# If you're not yet in a state that allows you to do the modal
# control optimization, instead of the previous line, do that:
gDM = 0.2
gTT = 0.3
nmod = M2V.shape[1]
optiGain = np.zeros( nmod )  # produces an array of gain
optiGain[0:2] = gTT
optiGain[2:]  = gDM

# Compute the optimized control matrix to be fed into the system

mc = computeSystemOptimizedControlMatrix( S2M, optiGain, M2V, TT2HO )


# Now .. write mc into a fits file, load it into the system, TT2HO also, and .. try,
test, etc.


"""
