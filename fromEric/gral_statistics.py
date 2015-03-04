# -*- coding: utf-8 -*-
"""
Created on Fri Dec 19 15:32:29 2014

@author: egendron



How to launch the procedures in this file:

# Before doing anything: prepare yourself.
import pyfits

# define where the files are located... (you'll need to update this line
# when running script on your machine !!)
path = '/home/egendron/Projets/GRAVITY/SVN_Packages/'


# First, get a Zernike reconstruction matrix made....
Z2S = pyfits.getdata(path+'Z2S_136s_119z.fits').T
norder=9
Jmax = ((norder+1)*(norder+2))/2
Z2S = Z2S[:,0:Jmax-1]
S2Z = np.linalg.pinv(Z2S)       # Inverts the Z2S it to get a slope-2-zernike matrix


# Define a few parameters ...........................
pixelSizeArcsec = 0.51  # arcseconds
DiamVLT = 8.00          # meters
Fs = 200.0              # Hetz
latency = 1.4           # frames

# read the interaction matrices ....
miaHO = pyfits.getdata(path+'HORecnCalibrat.RESULT_IM.fits')
miaTT = pyfits.getdata(path+'TTRecnCalibrat.RESULT.IM.fits')

# and read some data .....
a = pyfits.getdata(path+'LoopData_2.fits')
s = a.field(4)    # slopes
v = a.field(5)    # volts on HO DM
vtt = a.field(6)  # volts on TT


# READY !!!...  Real work can begin.
#
#

userExample(s, v, vtt, miaHO, miaTT, Fs, latency, S2Z, pixelSizeArcsec, DiamVLT, 2)


# That's all.

"""




import numpy as np
import scipy.optimize
import math

# This is just for plotting ...
import matplotlib.pyplot as plt



def gral_pseudoOpenLoopReconstruction(slopesdata, voltsdataHO, voltsdataTT, miaHO, miaTT, Fs, latency):
    """
    Input arguments :
    <slopesdata>  : set of N frames of slopes (i.e. centroids), shape (N,120)
    <voltsdataHO> : set of N frames of volts (i.e. DM commands) synchronized with slopes, shape (N, 60)
    <voltsdataTT> : set of N frames of volts (i.e. DM+TT commands) synchronized with slopes, shape (N, 2)
    <miaHO>       : measured interaction matrix (120x60, i.e. 120 lines, 60 columns)
    <miaTT>       : measured interaction matrix (120x2, i.e. 120 lines, 2 columns)
    <Fs>          : sampling frequency (Hertz)
    <latency>     : latency, i.e. time between the start of WFS integration and time when the
                    actuator reaches 50% of its command, in seconds.
    Returned value :
    open loop slopes data, as if they had been acquired in open-loop.

    This procedure applies both in open-loop or closed-loop.
    In the particular case of open-loop, all voltages voltsdata are 0, and the function will
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
    # determines how many frames are recorded in the data set,
    # will be used later (could also be passed as an argument..)
    nrec = slopesdata.shape[0]

    # slopes that were "added" by the action of the DM
    sv = np.dot(voltsdataHO, miaHO.T)
    sv += np.dot(voltsdataTT, miaTT.T)

    # those slopes <sv> need now to be subtracted from the residual slopes
    # with the proper delay.
    # But first one needs to process latency a bit.

    # latency expressed in frames (i.e. about 1.35 frames for gravity).
    latency_frame = Fs * latency
    # latency needs to be split into integer and fractional part,
    # i.e. 1.35 = 1 + 0.35
    i_latency_frame = int(latency_frame)                # integer part, i.e. 1
    f_latency_frame = latency_frame - i_latency_frame   # franctional part, i.e. 0.35

    # now <sv> will be shifted in time
    svs = np.zeros(sv.shape)   # just to allocate memory space same size as <sv>
    for i in range(nrec):
        j = i-i_latency_frame    # 1 frame before
        k = j-1                  # 2 frames before ...
        if(j<0): j=0             # except we don't have data before start of file !..
        if(k<0): k=0             # idem
        svs[i,:] = sv[j,:] * (1-f_latency_frame) + sv[k,:] * f_latency_frame

    # subtraction of time-shifted DM action from residual slopes data
    openloopSlopes = slopesdata - svs

    return openloopSlopes









def conversion_pixels_nmrms(pixelSizeArcsec, DiamVLT):
    """
    Input args:
    <pixelSizeArcsec> : scalar floating point = pixel size of WFS in arcsec = 0.51 (Pixel size of WFS is 0.51'')
    <DiamVLT>         : scalar floating point = VLT diameter in meters = 8.0 (VLT is 8m diameter)

    Output:
    scalar float value, that allows to transform zernike coeffs into nm rms.
    """
    # scaling factor to transform "nanometers rms" to "WFS pixels"
    RASC = 180*3600/3.14159265358979324    # number of arcsec in 1 radian
    unitfactor = 1e9 * DiamVLT * (pixelSizeArcsec/RASC) / 2.0
    return unitfactor




def varianceZernike_nm2( slopes , S2Z, pixelSizeArcsec, DiamVLT):
    """
    Input args:
    <slopes> : 2D array, array of slopes, shape (nframes, 136)
    <S2Z>    : 2D array, zernike reconstruction matrix
    <pixelSizeArcsec> : scalar floating point = pixel size of WFS in arcsec = 0.51 (Pixel size of WFS is 0.51'')
    <DiamVLT>         : scalar floating point = VLT diameter in meters = 8.0 (VLT is 8m diameter)

    Compute variance of zernike polynomials from an array of slopes <slopes>.

    """
    # zernike reconstruction, dimensions (nbzer, nbframes)
    szer = np.dot(S2Z, slopes.T)
    # varz = rms(szer)**2
    varz = np.var(szer, axis=1)
    # scale variances to get nanometres^2
    unitfactor = conversion_pixels_nmrms(pixelSizeArcsec, DiamVLT)
    varz *= unitfactor**2       # in nm^2
    return varz








def noise_inPix2(slopes, zero=True):
    """
    Input args:
    <slopes> : 2D array of slopes, shape is (nbframes, nbslopes=136)

    Output:
    A 1D array of nbslopes values (136 values), that are the noise variance on the 136 slopes
    The noise is expressed in pixels^2

    """
    npt = slopes.shape[0]          # number of frames in the slope set
    if( npt<10 ):             # meaningless result anyway if there are less than 10 frames
        return 0.0
    # Time-average of each slopes
    sa = np.average(slopes, axis=0)
    # remove average value of each slopes
    slopes_0 = slopes - sa
    # Variance
    a0 = np.average(slopes_0**2, axis=0)
    # Correlation with a shift of 1 frame
    a1 = np.average(slopes_0[0:npt-1,:]*slopes_0[1:npt,:], axis=0)
    # Correlation with a shift of 2 frames
    a2 = np.average(slopes_0[0:npt-2,:]*slopes_0[2:npt,:], axis=0)
    # parabolic fit
    varNoise = a0 - (4*a1-a2)/3.
    # negative values (theoretically impossible) are set to 0.00, by default.
    if(zero):
        varNoise[varNoise<0]=0.00

    return varNoise




def propagateNoiseOnZernike(varNoise, S2Z, pixelSizeArcsec, DiamVLT):
    """
    Input args:
    <varNoise> : 1D array of float, this is the output of function noise_inPix2()
    <S2Z> is a 2D array of float, it is the matrix "slope-to-zernike"
    <pixelSizeArcsec> : scalar float, pixel size in arcseconds, 0.51'' for gravity (TBC in AIT)
    <DiamVLT> : scalar float, diameter of VLT (meters) = 8.00

    Output:
    Returns an 1D array of floats, the noise variance, in nm^2 (propagated on the zernike).

    The input data <varNoise> is the array of the noise variance (in pix^2) for each slope.
    The noise variance is computed in pixels^2 on the slopes, and
    then propagated through the Zernike reconstruction matrix
    S2Z, and expressed in nm^2.


    """
    # scale variances to get nanometres^2
    unitfactor = conversion_pixels_nmrms(pixelSizeArcsec, DiamVLT)
    # propagation of each individual subaperture noise onto zernike coefficients
    zernoiseVar = np.dot( (S2Z**2), varNoise) * unitfactor**2
    return zernoiseVar       # returns the full zernike variances, nm^2





def meritFunction(fittedParams, data, i0, i1):
    """

    Input args:
    <fittedParams> is a 1D array of 2 parameters (D/r0 and D/L0), that will be
                   adjusted by the fitting procedure.
    <data>         is the 1D array of the data (floats, the measured zernike spectrum to be fitted).
    <i0> and <i1>  define the range of zernike indexes where the spectrum has to
                   be evaluated.

    Output: a 1D array of floats, same shape as <data>.

    This function computes the quadratic error between the
    model zernikeSpectrum_L0() and the <data>.
    The function zernikeKolmoSpectrum() returns the natural log of the spectrum of the Zernike
    variances, in the range of indexes i0 to i1.
    This function returns the difference between this spectrum and the data.

    CALLED FROM:
    This function is called by the procedure scipy.optimize.leastsq() in function getr0L0()

    """
    # Here we need to put an abs() around fitted parameters, because the fitting
    # procedure may set them to some negative value, while searching for a minimum. This negative value
    # would cause a floating-point interrupt --> hence the abs(), to prevent this.
    # It results is that the minimisation procedure may find, on output, a negative r0
    # and/or a negative L0, but no worry : we will just have to take the abs() of r0 and L0
    # at the very end, to get their right values.
    D_r0 = abs(fittedParams[0])
    D_L0 = abs(fittedParams[1])
    err = (zernikeSpectrum_L0(i0, i1, D_r0, D_L0, log=True, aliasing=False) - data)**2
    return err






def zernikeSpectrum_L0(i0, i1, dro, dlo, log=False, aliasing=False):
    """
    Input args:
    <i0> is the number of the 1st zernike for which the spectrum has to be computed
    <i1> is the number of the last zernike
    <dro> is the value of D/r0
    <dlo> is the value of D/L0
    <log> : boolean. If log==True the natural log of the values are returned.
    <aliasing> : boolean. If aliasing==True the spectrum is scaled to take aliasing
    into account. The influence of the aliasing of the 9x9 hartmann is simulated by
    multiplying the theoretical spectrum by coefficients (that have been calculated
    elsewhere) that simulate the aliasing impact.

    Output:
    A 1D array of floats, length (i1-i0+1).

    The function returns the theoretical spectrum of Zernike coeffs,
    with a given r0 and outer scale. It uses a formula given by Rodolph Conan in his phD thesis.
    The returned vector ranges from Zernike i0 to i1, i.e. there are i1-i0+1 coeffs.

    The function works for large L0>D/4, i.e. for 0<= (D/L0) < 4, it will diverge for L0<D/4.

    CALLED FROM:
    Function is called by meritFunction()
    """
    # The expression of the spectrum is computed thanks to a Taylor expansion series.
    # kmax is the number of terms of the taylor expansion.
    # kmax is set to 50, which is already a pretty large value .. it can be increased
    # if higher accuracy is required, e.g. for computing spectra for small L0<D/4,
    # however the success is not guaranteed.
    kmax=50
    # Number of coeffs computed
    ncoeffs = i1-i0+1
    # allocate memory to store results
    result=np.zeros(ncoeffs)

    pidf0=np.pi*dlo  # product pi*D*fo = pi*D/L0
    #  multiplicative factor  1.1628911592006026 :
    #  fact=2*gamma_R(11./6.)
    #  fact/=pi^(3./2.)
    #  fact*=((24./5)*gamma_R(6./5))^(5./6)
    dro53_116 = dro**(5./3) * 1.1628911592006026
    # Compute radial orders n of Zernike i0 and i1.
    nmin = int( (-1.+np.sqrt(8*(i0-1)+1))/2.)
    nmax = int( (-1.+np.sqrt(8*(i1-1)+1))/2.)

    # The convergence of the series is not garanteed for values of D/L0<4 (i.e. pi*D/L0<4*pi)
    # Function may return something stupid when D/L0<4.
    # A test was previously setup so that the function returns 0.
    # Now (EG, 21 feb 2014) the function limits the L0 values (clipping) to L0>D/4
    if pidf0>np.pi*4:
        pidf0=np.pi*4

    # Computation of the spectrum.
    # The variable <inn> is looping over radial orders. As all the zernike modes with the
    # same radial order have the same variance, the computation will be done only once per radial order,
    # and the variance will be copied onto all the elements of the array with the same
    # radial order. The Zernike numbers of a given radial order are determined by
    # J0 and J1.
    # Then J0 and J1 are limited to the range i0 to i1 (because this is where coeffs should be
    # computed)
    for inn in range(nmin,nmax+1):
        # computation of theoretical zernike variances of given radial order <inn>
        s = series_sum_diag(inn, kmax, pidf0) * (inn+1) * dro53_116
        # determine ranges of zernike index with this rad order
        J0 = (inn*(inn+1)) / 2 + 1   # first zernike index of given radial order inn
        J1 = ((inn+1)*(inn+2)) / 2   # last zernike index of given radial order inn
        # convert the zernike range of index into array index
        a0 = max(J0,i0) - i0         # array index of first coeff of rad ord inn
        a1 = min(J1,i1) - i0 + 1     # array index +1 of last coeff of rad ord inn
        # copy the computed value in the result array
        result[a0:a1] = s

    # Unbias from aliasing, if required.
    if( aliasing==True ):
        zerfactor = [1.00127,1.00127,1.0089,1.00421,1.01433,1.05697,1.05697,1.0464,1.0464,1.07974,1.12548,1.03053,1.07697,1.05695,1.0647,1.0647,1.13497,1.13497,1.12401,1.12401,1.46191,1.32155,1.30285,0.752463,1.25776,0.911428,0.994586,2.69474,2.69474,2.18867,2.18867,1.3627,1.3627,1.39158,1.39158]
        result = result * zerfactor

    # take natural log(), if required.
    if( log==True ):
        result = np.log(np.abs(result))

    return result







def series_sum_diag(n,kmax,pidf0):
    """
    Function coded by Fabrice Vidal. Derived from Yorick package "noll.i" written by E. Gendron.

    Input args:
    <n> : scalar integer. radial order of a Zernike mode
    <kmax> : scalar integer.
    <pidf0> : scalar float equal to pi*D/L0

    This function computes a sum of terms composed with Gamma
    functions, for the computation of zernike variances with finite
    L0. The formula is extracted from R. Conan thesis (equation 4.16
    page 116).
    The formula has been modified ton make use of the relation
    Gamma(k+x)=(k-1+x).Gamma(k-1+x) and thus avoid to compute the Gamma
    functions at each iteration.

    CALLED FROM:
    this function is called by zernikeSpectrum_L0()
  """

    n2 = 2*n

    # fn = factorial(n)
    fn = 1.00
    # compute (2+n1+n2)! using the previous results again
    for i in range(2,2+n2+1):
        fn*=i
    # computation of all gamma() coefficients
    # UDO: if you don't have a C-library to compute gamma functions, I can provide
    # you with a simple algorithm to do that. Just email me.
    pregamma = [math.gamma(i) for i in [1.5 + n, 5./6 - n, n-5./6, n+23./6]]

    # Initialisation
    uk = pregamma[0] * pregamma[1] / ((n+1) * fn)   # u0

    #  0.6494396447776356 = gamma_R(7/3.)*gamma_R(11./6)/gamma_R(17./6)
    vk = pregamma[2]*0.6494396447776356 / pregamma[3]

    pidf0_n253 = pidf0**(n2-5/3.)
    pidf0_2 = pidf0*pidf0
    pidf0_2k = 1.00                # preparing pidf0^(2*k), for k=0
    fk = 1.0                       # preparing k!, for k=0

    s = (uk * pidf0_n253 + vk)     # s(k=0)

    for k in range(1, kmax+1):
        fk *= k                    # k!
        pidf0_2k *= pidf0_2        # computation of  pidf0^(2*k)

        uk *= ((0.5+n+k)*(k+n))/((2+n2+k)*(1+n+k)*(5./6-n-k))
        vk *= (k+4./3)*(k+5./6)/((n-5./6-k)*(n+17./6+k)*(11./6+k))

        tmps = (pidf0_n253 * uk + vk) * pidf0_2k / fk

        if(k%2):
            s -= tmps
        else:
            s += tmps


    return s







def getr0L0( varz, noiseOnZernike, pixelSizeArcsec, DiamVLT, i0=4, i1=55):
    """
    <varz>            : 1D array of openloop zernike variances, or reconstructed
                        pseudo-openloop zernike variances, coming from the output
                        of varianceZernike_nm2(), expressed in (nm rms)^2
    <noiseOnZernike>  : 1D array of float. It is the output of propagateNoiseOnZernike()
    <pixelSizeArcsec> : pixel size in arcseconds, 0.51'' for gravity (TBC in AIT)
    <DiamVLT>         : scalar float, diameter of VLT (meters) = 8.00
    <i0> and <i1>     : scalar integers, with default values 4 and 55 for Gravity-CIAO.

    Outputs:
    This function computes r0 and L0 (both scalar floats), both in meters.
    It is doing a non-linear least-square fit of the measured Zernike spectrum
    against a modelled spectrum with r0 and L0 as free parameters.

    """
    # Subtraction of propagated noise from the measured Zernike variances
    #
    # It may happen that for very HIGH noise conditions, the variance of the Zernike is
    # completely dominated by noise. Then, it may happen that the estimation of the variance
    # becomes negative after subtracting the propagated noise.
    # Those negative values need to be eliminated: we select them, and replace them by the
    # initial variance, divided by 100
    n = varz.shape[0]
    varz_noNoise = np.zeros(varz.shape)
    for i in range(n):
        if varz[i]>noiseOnZernike[i]:
            varz_noNoise[i] = varz[i] - noiseOnZernike[i]
        else:
            varz_noNoise[i] = varz[i] / 100.0

    # conversion from (nm^2) to (rd rms)^2
    lambda_vis = 500.0                    # r0 will be expressed at lambda=500 nm.
    nm2rd =  (2*np.pi/lambda_vis)**2      # conversion coeff from nm to radians (at lambda=500nm)
    varz_noNoise *= nm2rd**2              # now varz_noNoise is in (radians^2)

    # taking the natural log of the data, starting at Zernike index i0 (included) and finishing at i1 (included).
    # As the S2Z produces zernike coefficients starting at index=2 (tip), one need to subtract 2.
    data = np.log(varz_noNoise[i0-2:i1-1])

    # Now we will fit the coefficients (D/r0) and (D/L0): initial guess will be r0=12cm and L0=20m (just because it makes sense)
    initguess_Dr0L0 = [8.0/0.12, 8.0/20.0]
    # Non-linear least-square fit using a Levernberg-Marquardt algorithm
    # General information: http://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
    # GPL C routines : http://users.ics.forth.gr/~lourakis/levmar/
    Dr0L0Fitted = scipy.optimize.leastsq(meritFunction, initguess_Dr0L0, args=(data, i0, i1))[0]
    # We have put an abs() around the fitted parameters, because the fitting
    # procedure may potentially set them to a negative value while searching for a minimum. This negative value
    # would cause a floating-point interrupt --> hence the abs().
    # It results from this that the minimisation procedure may find, on output, a negative r0
    # and/or a negative L0, but no worry : we just have to take the abs() of r0 and L0
    # at the end, to get their right values.
    Dr0L0Fitted = np.abs(Dr0L0Fitted)

    r0 = DiamVLT/Dr0L0Fitted[0]
    L0 = DiamVLT/Dr0L0Fitted[1]

    return r0, L0







def getTau0(slopes, noiseOnSlopes, samplingFrequency):
    """
    Input args:
    <slopes>           : 2D array of floats. Slopes, with a shape (Nframes, Nslopes=136).
    <noiseOnSlopes>    : 1D array of floats (shape (Nslopes=136)), which is the
                         output of noise_inPix2(slopes, zero=True)
    <sampligFrequency> : scalar float. Sampling frequency, in Hertz.

    Output:
    Scalar value, tau0, in seconds.

    The function computes the correlation time of slopes.
    It computes the autocorrelation for each slope, it subtracts the noise variance
    from the central point, and sums all the autocorrelation to get an average.
    Then, it searches for the time required to reach the half-maximum.
    A series of 5 points bracketing the half-maximum are selected, and a linear regression
    is done in order to locally fit the autocorrelation curve nicely around those 5 points,
    and then find a proper interpolation.
    """
    # number of frames of the array <slopes>
    n, nslopes = slopes.shape
    # Just a little check ... if there are too few values, the result is meaningless
    # and the rest of the program will crash anyway
    if n<10:
        print 'Too few values, tau0 cannot be computed.'
        return 0.00
    # allocate memory for intermediate result
    autocorr = np.zeros(n)
    # compute autocorrelations of slopes, and sum everything. Autocorrelation is computed by
    # taking the Fourier transform, then square of the modulus, then Fourier transform again,
    # and ignore imaginary part (keep real part only)
    for i in range(nslopes):
        # compute time-average of slope variation
        sa = np.average(slopes[:,i])
        autocorr += np.fft.fft(np.abs(np.fft.fft(slopes[:,i]-sa))**2).real / n**2
        autocorr[0] -= noiseOnSlopes[i]

    # Now normalizing the max of the autocorrelation to 1.00
    if autocorr[0]==0 or autocorr[0]<0:
        print 'Error. This should not happen .. data corrupted or all null ? Please check'
        return 0.00

    # Search for the first point just below half-maximum, starting at index k=2
    n /= 2
    k=2
    while (autocorr[k]>(0.5*autocorr[0]) and (k<n)):
        k += 1              # search first point below 0.5

    # linear regression on 5 points around the half-max point (k-2, k-1, k, k+1, k+2)
    if( k>3 and (k+2)<n ):
        y = autocorr[k-2:k+3]   # extract those 5 points
        x = np.arange(k-2,k+3)  # create a vector x=[k-2, k-1, k, k+1, k+2]
        # Formula of the linear regression of y versus x
        # tmp = (np.average(x*x)-np.average(x)^2.) : tmp is always equal to 2 in our case, whatever the value of k.
        tmp = 2.0
        a = (np.average(y*x)-k*np.average(y)) / tmp
        b = np.average(y) - a*k
        if  a<0 :
            # success...
            k = (0.5-b)/a      # solve equation a*k+b = 0.5
        else:
            # big fail ...
            print 'Error. Unable to determine tau0, weird autocorrelation ... Please check data.'
            print k, a, b
            #return 0.00
    else:
        print 'Correlation time tau0 is too short. Cannot measure it properly.'
        print 'Were those data *really* openloop ones ? .. please check.'
        print 'Otherwise ask Gendron for a better algorithm.'
        return 0.00

    tau0 = k / samplingFrequency   # tau0 in seconds

    return tau0







"""
         _____ _    _   ___  __
        |  ___| |  | | | \ \/ /
        | |_  | |  | | | |\  /
        |  _| | |__| |_| |/  \
        |_|   |_____\___//_/\_\

"""


def computeFluxAverage( fluxBuffer ):
    """
    Input args:
    <fluxBuffer>    : a 2D array of flux, as acquired by sparta

    Output:
    scalar value, the average flux (in ADU)

    """

    # really average everything in space and time
    avgFlux = np.average(fluxBuffer)

    return avgFlux












"""
 ____  __  __     ____    _  _____ _   _ ____      _  _____ ___ ___  _   _
|  _ \|  \/  |   / ___|  / \|_   _| | | |  _ \    / \|_   _|_ _/ _ \| \ | |
| | | | |\/| |   \___ \ / _ \ | | | | | | |_) |  / _ \ | |  | | | | |  \| |
| |_| | |  | |    ___) / ___ \| | | |_| |  _ <  / ___ \| |  | | |_| | |\  |
|____/|_|  |_|   |____/_/   \_\_|  \___/|_| \_\/_/   \_\_| |___\___/|_| \_|





"""

def checkDmSaturation(saturationValue, volts):
    """
    Input args:
    <saturationValue> : scalar float, value of saturation
    <volts> : 2D array of shape (N, 60)

    The MACAO DM is made of 5 circular rings of electrodes. The saturation of electrodes
    will be checked for 2 lists of electrodes : the outer ring with 20 electrodes (41 to 60)
    and all the other ones (inner ones, 40 electrodes, 1 to 40).

    WARING NOTE: this function has been written, checked for functionnality, but
    never tested against real data.

    """
    # get dimensions of array volts
    nframes, nact = volts.shape
    if( nact!=60 ):
        print "AAAH !!! stop here .. the rest of this procedure will not work anyway !!"
        return 0
    # allocte memory
    count = np.zeros(nact)
    # count number of saturations per electrode
    for i in range(nact):
        count[i] = np.sum( np.abs(volts[:,i])>=saturationValue )
    # counting saturations for "inner electrodes" and "outer ones" ...
    # sorry for the rough hardcoding here. If you have better idea, please do so ... !!
    countInner = np.sum(count[0:40]) / nframes * 100.0
    countOuter = np.sum(count[40:60]) / nframes * 100.0

    return (countInner, countOuter)









"""



 _____ ___  _   _ ____  ___ _____ ____
|  ___/ _ \| | | |  _ \|_ _| ____|  _ \
| |_ | | | | | | | |_) || ||  _| | |_) |
|  _|| |_| | |_| |  _ < | || |___|  _ <
|_|   \___/ \___/|_| \_\___|_____|_| \_\

    _    _   _    _    _  __   ______ ___ ____
   / \  | \ | |  / \  | | \ \ / / ___|_ _/ ___|
  / _ \ |  \| | / _ \ | |  \ V /\___ \| |\___ \
 / ___ \| |\  |/ ___ \| |___| |  ___) | | ___) |
/_/   \_\_| \_/_/   \_\_____|_| |____/___|____/




These routines are made for vizualization of the data, by the user.




"""






def FFThz(signal, samplingFreq):
    """
    Input args:
    <signal>       : a 1D array of floats
    <samplingFreq> : scalar float, sampling frequency in Hz.

    Output:
    Power spectral density (PSD) of the input signal. This function
    takes care of the units of the output PSD. If signal has units 'u',
    the PSD is expressed in units 'u^2 per Hz'.
    Warning: only the half of the frequencies, plus excluding null
    frequency, i.e.(from Fs/n to Fs/2), are returned.

    """
    n = signal.size
    d = np.abs(np.fft.fft(signal))[1:n/2+1]       # we drop the point at f=0, and points at f>Fs/2
    d = d**2/(samplingFreq*n/2)                   # proper normalization to get units 'u^2/Hz'
    d[n/2-1]/=2                                   # the 'nyquist point' (=last point) needs to be divided by 2
    return d




def createFrequencyAxis(npt, samplingFreq):
    """
    Input args:
    <npt> : scalar integer, number of points.
    <samplingFrequency> : scalar float, sampling frequency in Hz.

    Output:
    a 1D array of floats, which are the frequency axis that is matching the PSD samples
    created by function FFThz().
    First element is (samplingFreq/npt), last element is (samplingFreq/2). There are (npt/2) elements in the array.


    """
    # return np.linspace(0,samplingFreq,npt+1)[1:npt/2+1]          # works fine, but heavy
    # return np.linspace(samplingFreq/npt,samplingFreq/2,npt/2)    # works fine, better optimized
    return (1+np.arange(npt/2)) * samplingFreq/npt                 # best formulation for a more direct C-translation ...





def doflat( signal ):
    """
    Input args:
    <signal>       : a 1D array of floats

    Output:
    1D array of floats, same length as <signal>

    """
    n = signal.size
    # check length of signal is not stupid .. and allows avoir the
    # floating-point error of dividing by 0 if n==1 (which should
    # never happen, by the way.. but who knows ... )
    if n<3:
        print "Warning: you're attempting to compute an FFT on a signal of length %d. Please check ..!" % n
        return signal
    a = signal - (np.arange(n)*(signal[n-1]-signal[0])/(n-1) + signal[0])
    return a






def dowin( signal ):
    """
    Input args:
      <signal>       : a 1D array of floats

    Output:
      1D array of floats, same length as <signal>

    The function multiplies the signal by a window with a shape in window=sin()^2,
    with window[edges]=0, and window[middle]=1.
    """
    n = signal.size
    # check length of signal is not stupid .. and allows avoir the
    # floating-point error of dividing by 0 if n==1 (which should
    # never happen, by the way.. but who knows ... )
    if n<3:
        print "Warning: you're attempting to compute an FFT on a signal of length %d. Please check ..!" % n
        return signal
    # application of a window
    a = signal * np.sin(np.arange(n)*(np.pi/(n-1)))**2
    # The signal is multiplied by sin()^2, which means that the square of the
    # signal is multiplied by sin()^4. The integral from 0 to 1 of sin(pi.t)^4
    # is 3/8, that's why we need to divide by this value in order to re-normalize the total
    # variance. However this factor applies on the PSD, while here we're dealing with
    # the signal: so we just need to compensate by a factor of sqrt(8/3)=1.63299316185545
    a *= 1.63299316185545
    return a







def buildPSD( signal, samplingFreq, flat=False, win=True, K=1024 ):
    """
    Input args:
      <signal>       : a 1D array of floats,
      <samplingFreq> : scalar float, sampling frequency in Hz.

    Output:
      Returns a 1D array of floats, the PSD of <signal>.

    The final PSD is obtained using Welsh method, cutting the data set into
    subsets of K=1024 frames (by default), overlapping by K/2=512 samples,
    and windowed with a sin²(πt/N), and making the PSD for each and
    averaging all of them.

    The advantage is
      - a greater accuracy of the spectrum (less 'noise'), thus better identification of spectrum features (noise, peaks, ..)

    drawbacks are
      - a (slightly) coarser resolution in frequency,
      - absence of points at very low frequencies


    """
    # length of signal
    N = signal.size
    # if the signal is shorter than 1024 samples, then we do not attempt to cut it .. it's too short anyway.
    if( N<=K ):
        K = N

    # allocate memory for result. The output DSP will have K/2 points.
    psd = np.zeros(K/2)
    k = 0
    # this is where we're potentially gonna cut our signal
    ind = range(0, N, K/2)  # with N=2048 and K=1024,  ind=[0, 512, 1024, 1536]
    # let's cut, and compute our PSDs
    for i in ind:
        # this test is to ensure the subset we're gonna cut does actually exists ...
        if( i+K<=N ):
            tmp = signal[i:i+K]   # extract K points from the signal
            # prevent signal to 'jump' between edges, to avoid introducing fake high frequencies
            if( flat==True ):
                tmp = doflat(tmp)
            # apply a sin^2 windowing
            if( win==True ):
                tmp = dowin(tmp)
            # compute power spectral density
            psd += FFThz ( tmp, samplingFreq )
            # count how many PSD have been computed
            k += 1
    # divide by the number of additions (the test is just here to avoid a floating-point
    # exception, although this should never happen here: k MUST at least be equal to 1)
    if k>0:
        psd /= k
    return psd





def buildCumulativePSD( psd, samplingFreq ):
    """
    Input args:
    <psd>          : a 1D array of floats given by buildPSD()
    <samplingFreq> : scalar float, sampling frequency in Hz.

    Output:
    Returns a 1D array of floats, the square root of the frequency-integrated PSD.

    """
    cpsd = np.cumsum(psd)
    n = psd.size
    cpsd *= samplingFreq/2/n
    cpsd = np.sqrt(cpsd)
    return cpsd






def estimateFittingError_CIAO(r0):
    """
    Input args:
    <r0>         : scalar floating point = value of r0 (gieen at 500nm), in meters.

    """
    pitchMACAO = 1.0   # crude approximation.. will be refined during AITs
    # Computation of the fitting error, in rd^2 at 500 nm.
    fittingErr = 0.23 * (pitchMACAO/r0)**(5./3.)
    # Now we convert the rd^2 to nm rms
    lambda_vis = 500.0    #
    fittingErr *= (lambda_vis/2/np.pi)**2
    return fittingErr





def SQRT(x):
    """
    <x> : scalar float
    The function returns the square root of x when positive, and -sqrt(-x) when x<0.
    """
    if x<0:
        return -np.sqrt(-x)
    else:
        return np.sqrt(x)




def computeStrehl(lambdaIR, varTotal, varLow, r0, DiamTel):
    """
    Input args:
    <lambdaIR> : scalar float, IR wavelength (in nanometers)
    <varTotal> : scalar float, total variance (in nm^2) of the residual wavefront
    <varLow>   : scalar float, variance (in nm^2) of low order modes (tip & tilt)
    <r0>       : scalar float, r0 (in meters) at 500 nm
    <DiamTel>  : scalar float, diameter of the telescope (in meters)

    The procedure computes the Strehl ratio using the approximation given by
    (Parenti & Sasiela, JOSA 11, No 1, 1994).
    This approximation assumes that the PSF is an airy pattern with a SR given by the
    Marechal approximation (exp(-sigma^2)) from the variance of high order modes, convolved
    by the jitter of tip and tilt that reduce the previous Strehl by an amount 1/(1+sigma_tt^2),
    plus a term due to the halo in 1/(1+(D/r0IR)^2) applicable to all the energy which is not in
    the previous airy pattern, i.e. (1-SR).
    When the variance of high order modes is large, the expression tends towards
    1/(1+(DiamTel/r0_IR)**2) which is the expression given by (Yura, JOSA 63, 1973).

    """
    # conversion factor nm^2 to radians^2
    nm2_to_rd2 = (2*np.pi/lambdaIR)**2
    # computation of variance of high orders
    varHigh = varTotal - varLow
    # Strehl Ratio High Orders
    SRHO = np.exp(-varHigh * nm2_to_rd2)
    # scaling r0 in the IR
    lambdaVis = 500.    # r0 is given at 500 nm
    r0_IR = r0 * (lambdaIR/lambdaVis)**(6./5)
    # Strehl computation
    SR = SRHO / (1 + varLow * nm2_to_rd2) + (1-SRHO)/(1+(DiamTel/r0_IR)**2)
    # translate this in percent
    SR *= 100.0
    return SR







def userExample(s, v, vtt, miaHO, miaTT, Fs, latency, S2Z, pixelSizeArcsec, DiamVLT, izer):
    """
    <s>           : set of N frames of slopes (i.e. centroids), shape (N,120)
    <v>           : set of N frames of volts (i.e. DM commands) synchronized with slopes, shape (N, 60)
    <vtt>         : set of N frames of volts (i.e. DM+TT commands) synchronized with slopes, shape (N, 2)
    <miaHO>       : measured interaction matrix (120x60, i.e. 120 lines, 60 columns)
    <miaTT>       : measured interaction matrix (120x2, i.e. 120 lines, 2 columns)
    <Fs>          : scalar floating point = sampling frequency (Hertz)
    <latency>     : scalar floating point = latency (in seconds), i.e. time between the start of WFS
                    integration and time when the actuator reaches 50% of its command, in seconds.
    <S2Z>         : 2D array, zernike reconstruction matrix
    <pixelSizeArcsec> : scalar floating point = pixel size of WFS in arcsec = 0.51 (Pixel size of WFS is 0.51'')
    <DiamVLT>     : scalar floating point = VLT diameter in meters = 8.0 (VLT is 8m diameter)
    <izer>        : scalar integer, a zernike number

    This function is just an example that shows how to use all the
    functions in this file.

    """

    # Do the pseudo openloop reconstruction
    POLslopes = gral_pseudoOpenLoopReconstruction(s, v, vtt, miaHO, miaTT, Fs, latency)

    # compute the variance of the Zernike
    POLvarz = varianceZernike_nm2( POLslopes , S2Z, pixelSizeArcsec, DiamVLT)
    varz    = varianceZernike_nm2( s         , S2Z, pixelSizeArcsec, DiamVLT)

    # compute the noise on the POLslopes ..
    noiseOnSlopes = noise_inPix2(POLslopes, zero=True)

    # propagate this noise on Zernike modes
    noiseOnZernike = propagateNoiseOnZernike(noiseOnSlopes, S2Z, pixelSizeArcsec, DiamVLT)

    # Now compute r0 and L0:
    r0, L0 = getr0L0( POLvarz, noiseOnZernike, pixelSizeArcsec, DiamVLT, i0=4, i1=55)

    # and compute correlation time
    tau0   = getTau0(POLslopes, noiseOnSlopes, Fs)

    # check saturation of DM electrodes
    saturationValue = 9.5   # saturation occurs after 9.50 volts
    countInner, countOuter = checkDmSaturation(saturationValue, v)

    # Display some useful statistical data:
    print "r0 at 500 nm          %.2f cm" % (r0*100.0)
    print "Seeing                %.2f ''" % (0.103/r0)
    print "Outer scale           %.1f m" % (L0)
    print "Correlation time      %.1f ms" % (tau0*1000.0)
    print "DM saturation         %.1f/%.1f (inner/outer)" % (countInner, countOuter)
    # Noise ...
    avgNoiseOnSlopes = np.average(noiseOnSlopes)
    print "Noise on slopes       %.2f pixels rms (%.2f arcsec rms)" % (np.sqrt(avgNoiseOnSlopes),np.sqrt(avgNoiseOnSlopes)*pixelSizeArcsec)
    wfsNoiseLevel = np.sum(noiseOnZernike)
    print "WFS Noise level       %.0f nm rms" % (np.sqrt(wfsNoiseLevel))
    # Estimation of tip tilt residuals
    tiltResidualError = (varz[0]+varz[1]) - (noiseOnZernike[0]+noiseOnZernike[1])  # expressed in nm^2, for 2 axis (tip and tilt summed)
    RASC = 180*3600/np.pi   # conversion factor radian to arcsec
    tiltResidualError_arcs = (RASC*4e-9/DiamVLT)**2 *  tiltResidualError / 2.0     # translation in arcsec^2, for 1 axis only
    print "Tilt residual error   %.2f ''" % (SQRT(tiltResidualError_arcs))
    print "Tilt residual error   %.1f %% of Airy disc" % (100.0*SQRT(tiltResidualError_arcs) / (RASC*2200.e-9/DiamVLT))
    # Fitting error
    fittingError = estimateFittingError_CIAO(r0)
    print "Fitting error         %.0f nm rms" % (np.sqrt(fittingError))
    # Total errors
    wfsRawTotalError = np.sum(varz)
    print "WFS total raw error   %.0f nm rms (incl. measured noise)" % (np.sqrt(wfsRawTotalError))
    print "WF error from WFS     %.0f nm rms (same w/o noise)" % (SQRT(wfsRawTotalError - wfsNoiseLevel))
    CiaoError = wfsRawTotalError - wfsNoiseLevel + fittingError * 1.30
    print "CIAO error            %.0f nm rms (full WF error: BW, noise, fitting, aliasing)" % (SQRT(CiaoError))
    lambdaIR = [1250.0, 1650.0, 2200.0]    # wavelength band J, H, K in nanometers
    for lam in lambdaIR:
        print "SR at %d nm         %.0f %%" % (int(lam), computeStrehl(lam, CiaoError, tiltResidualError, r0, DiamVLT))


    # Plot graphs of Zernike spectrum (i.e. Zernike variances) of closed-loop and pseudo openloop data
    plt.figure(1)
    plt.clf()
    plt.subplot(2,2,1)
    zindex = np.arange(varz.size) + 2   # +2 because the 1st zernike is Z2=tip, then Z3=tilt, etc.
    plt.plot(zindex, np.sqrt(varz), color='red')
    plt.plot(zindex, np.sqrt(POLvarz))
    plt.xlabel('Zernike number')
    plt.ylabel('Stdev (nm rms)')
    plt.xlim(zindex[0],zindex[-1])
    plt.yscale('log')
    plt.grid(True)

    plt.subplot(2,2,3)
    plt.plot(zindex, np.sqrt(np.cumsum(varz)), color='red')
    plt.plot(zindex, np.sqrt(np.cumsum(POLvarz)))
    plt.xlabel('Zernike number')
    plt.ylabel('Cumulative stdev (nm rms)')
    plt.yscale('linear')
    plt.xlim(zindex[0],zindex[-1])
    plt.grid(True)

    # We will need to convert slopes (in pixels) to Zernike (in nm rms) and for that
    # purpose we need a conversion factor
    unitfactor = conversion_pixels_nmrms(pixelSizeArcsec, DiamVLT)

    # Computing power spectral densities of a given Zernike mode
    # let's choose, for instance, Zernike number 4  (4 is defocus)
    # .. err .. no, I prefer to pass <izer> as an argument
    K = 512   # why not
    f = createFrequencyAxis(K, Fs)       # compute frequency axis, first.
    signal = np.dot(POLslopes, S2Z[izer-2,:]) * unitfactor   # getting temporal signal of Zernike. The [izer-2] is required because the S2Z starts at Z_2 (tip).
    psd = buildPSD( signal, Fs, flat=True, win=True, K=K )
    signalCL = np.dot(s, S2Z[izer-2,:]) * unitfactor
    psdCL = buildPSD( signalCL, Fs, flat=True, win=True, K=K )

    # now trying to plot things ...
    plt.subplot(2,2,2)
    plt.plot(f, psd)
    plt.plot(f, psdCL, color='red')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power spectral density $(nm^2/Hz)$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(f[0],f[-1])
    plt.grid(True)


    plt.subplot(2,2,4)
    plt.plot(f, buildCumulativePSD(psd, Fs))
    plt.plot(f, buildCumulativePSD(psdCL, Fs), color='red')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Cumulative PSD (nm rms)')
    plt.xscale('log')
    plt.grid(True)
    plt.xlim(f[0],f[-1])

    plt.show(block=False)
