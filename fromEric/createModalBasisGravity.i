

/*


  These comments are the list of commands you should run in a Yorick
  terminal to be able to use these functions.
  

// Reads useful basic data : an interaction matrix and a matrix
// <delta> that is stored somewhere on the disk, pre-computed during
// AITs.

miaHO = fits_read("simu_mia_HO.fits");
miaTT = fits_read("simu_mia_TT.fits");
delta = fits_read("simu_delta.fits");



// Creates modal basis and outputs <TT2HO> (matrix of the tilts of the
// DM (will be used in sparta)), <M2V> (matrix of system modes
// (will be used in the modal control optim and elsewhere)), <piston>
// (piston mode (will be used in sparta)), and also <S2M> the modal
// control matrix, i.e. a control matrix that multiplies to slopes
// vectors and produces the modes coefficients.

gral_createModalBasis, miaHO, miaTT, delta, TT2HO, M2V, piston, S2M;




// Reads (or acquires..) a set of slopes and voltages
s = fits_read("simu_slopeset.fits");
v = s(1:62,) * 0.0;   // volts are just 0 here
s += random_n( dimsof(s) )*40;   // add random numbers on slopes to simulate noise


Fs = 200.0;        // Sampling frequency (Hz), please replace by yours.
latency = 0.008;   // Latency (s), please replace by yours. If you don't know, keep 0.008 s.
BP = 700.;         // MACAO DM first resonnant frequency (Hz), should be ok.


// If you're able to run the modal control optimization, do this :
//
// Perform modal control optimization: optiGain is the array of the modal gains.
optiGain = gral_modalControlOptimization( s, v, miaHO, miaTT, S2M, M2V, 0.8, Fs, latency, BP );


// If you're not yet in a state that allows you to do the modal
// control optimization, instead of the previous line, do that:
gDM = 0.2;
gTT = 0.3;
optiGain = array(gDM, numberof(M2V(1,)));  // produces an array of gain
optiGain(-1:0) = gTT;


// Compute the optimized control matrix to be fed into the system

mc = gral_computeSystemOptimizedControlMatrix( S2M, optiGain, M2V, TT2HO );


Now .. write mc into a fits file, load it into the system, TT2HO also, and .. try, test, etc.


 */


func gral_generalized_inverse(mat , nfilt, ploteigenval=)
/* DOCUMENT gral_generalized_inverse( mat, nfilt, ploteigenval=)
   
   <mat>   : matrix to be inverted
   <nfilt> : number of modes to be filtered out
   <ploteigenval> : optional flag, for plotting eigenvalues for the user (because
                    users are usually keen on having eigenvalues displayed).

   The function computes the generalized inverse of a matrix, using a
   singular value decomposition.
   The returned matrix has a transposed shape wrt input matrix <mat>.

   The way this routine is coded here may be dependent on the choice
   of either row-major (also said C-contiguous) or column-major
   (Fortran contiguous) (see http://en.wikipedia.org/wiki/Row-major_order)
   that is made by the chosen langage.
   See for instance:
   Numerical Recipes (Press, et. al. Cambridge University Press 1988)
   has a good discussion of how to use the SVD -- see section 2.9.
   or just
   http://en.wikipedia.org/wiki/Singular_value_decomposition
   
   The package NumPy provides a pseudo-inverse calculation through its
   functions matrix.I and linalg.pinv; its pinv uses the SVD-based
   algorithm. SciPy adds a function scipy.linalg.pinv that uses a
   least-squares solver. High quality implementations of SVD, QR, and
   back substitution are available in standard libraries, such as
   LAPACK.
   
*/
{
  local u;
  local vt;
  // singular value decomposition. Produces 2 matrices (u and vt) and eigenvalues (s).
  // The SVdec allows us to rewrite the matrix mat as
  // mat = u . s . vt
  s = SVdec(mat,u,vt);
  
  // plot things, just for the purpose of debugging this code. Otherwise useless.
  if( ploteigenval==1 ) {
    fma; logxy,0,1;
    na = dimsof(mat)(0);
    plg, [lam(0),lam(1)], [1,1]*(na-nfilt+0.5), color="red", marks=0;
    plg, lam, marks=0; plmk, msize=0.3, lam, marker=4;
  }
  
  // inversion of eigenvalues, letting to 0 those filtered out.
  if( nfilt>0 ) {
    s1 = s;
    // this is to avoid any 'floating point interrupt', as some of the
    // filtered eigenvals may be equal to 0.
    // Syntax warning: in Yorick, arr(-3:0) means "the last 4 elements of the
    // array <arr>. Index '0' stands for the last element of the array.
    s1(1-nfilt:0) = 1.0;
    // inversion of all eigenvalues
    s1 = 1./s1;
    // setting back to 0 all the 
    s1(1-nfilt:0) = 0.0;
  } else {
    s1 = 1.0 / s;
  }
  // computing the inverse now. The inverse is given by
  // 1/mat = transpose(vt) . s1 . transpose(u)
  m1 = vt(+,) * (u*s1(-,))(,+);
  return m1;
}






func gral_identite(n)
/* DOCUMENT  idmat = gral_identite(n)
     Produces a nxn identity matrix
     
   SEE ALSO:
 */
{
  // memory allocation, filled with 0.00
  a = array(0.0,n,n);
  // Set diagonal to 1
  for(i=1;i<=n;i++) a(i,i) = 1.0;
  return a;
}







func gral_filterOutPiston( modes, piston, delta )
/* DOCUMENT fmodes = gral_filterOutPiston( modes, piston, delta );
   Input arguments :
     <modes> is an array containing N modes, it is an array 60xN
     <piston> is the piston mode defined on actuator voltages (60 components)
     <delta> is the geometric cov matrix of the DM (60x60), fixed.

   The function will suppress piston from all the modes of matrix <modes>.
   
   SEE ALSO:
 */
{
  pnorm = piston(+) * (delta(,+)*piston(+))(+);
  proj = modes(+,) * (delta(,+)*piston(+))(+);
  proj /= pnorm;
  fmodes = modes - piston(,-)*proj(-,);
  return fmodes;
}





func gral_createModalBasis( miaHO, miaTT, delta,   &TT2HO, &M2V, &piston, &S2M )
/* DOCUMENT  gral_createModalBasis, miaHO, miaTT, delta,   TT2HO, M2V, piston, S2M;
     
     Input arguments :
     <miaHO>   is the measured interaction matrix of DM (120x60, i.e. 120 lines, 60 columns)
     <miaTT>   is the measured interaction matrix of TT (120x2, i.e. 120 lines, 2 columns)
     <delta>   is the geometric covariance matrix of actuators, it is computed elsewhere.
               It is a square, symmetric matrix 60x60

     Outputs arguments :
     <TT2HO>    : matrix of the tilts of the DM (will be used in sparta)
     <M2V>      : matrix of system modes (will be used in the modal control optim and elsewhere)
     <piston>   : piston mode (will be used in sparta)
     <S2M>      : modal control matrix, that will be used in particular in the modal
                  control optimization.
                  

     note: see ESO doc
      VLT-SPE-ESO-16100-4856_iss1_SPARTA-Light Reference Implementation Technical Specifications

   SEE ALSO:
 */
{
  // get number of slopes <ns> and number of actuators of DM <na>
  ns = dimsof(miaHO)(2);   // number of slopes measured by the Shack-Hartmann, ns should be 120
  na = dimsof(miaHO)(3);   // number of actuators of HO DM, na should be 60
  

  // generalized inverse of miaHO
  nfilt = 8;   // to be adjusted during AITs
  miaHO1 = gral_generalized_inverse(miaHO , nfilt);       // miaHO1 has transposed shape wrt miaHO
  //lam = SVdec(miaHO,U,VT);     

  // Computation of set of DM voltages that reproduce a pure tilt
  TT2HO = miaHO1(,+)*miaTT(+,);
  
  // Computation of associated measurements
  Mtilt = miaHO(,+) * TT2HO(+,);
  
  
  // Compute tilt-filtering matrix
  // The tilt which is filtered out is the DM-TT.
  // Mtilt.inverse(transpose(Mtilt).Mtilt).transp(Mtilt)
  tmp = Mtilt(,+) * (LUsolve(Mtilt(+,)*Mtilt(+,))(,+) * Mtilt(,+))(+,);
  // now the DMtilt-filtering matrix
  filterOutTilt = (gral_identite(ns) - tmp);     
  
  // now filter out any tilt from the interaction matrix
  miaHOf = filterOutTilt(,+) * miaHO(+,);

  // Diagonalizing matrix (transpose(miaHOf).miaHOf) and getting eigenvectors
  // Two eigenvalues=0 appear here, due to 2 filtered modes
  lam = SVdec(miaHOf(+,)*miaHOf(+,), miaHOf_eigenvect);


  // Now removing piston
  // I decide that piston is mode number 58 (as 59 and 60 are the 2 filtered tilts).
  // I've checked this is extremely close to reality ..
  // The SVdec has produced modes that are orthogonal to piston in the
  // measurement space, but not in the phase space. Now, one has to
  // remove the piston phase component to all other modes.
  piston = miaHOf_eigenvect(,58);   // i'm choosing this here .. but piston could be something else. To be tuned during AITs.
  DMmodes = gral_filterOutPiston( miaHOf_eigenvect, piston, delta(1:60,1:60) );
  
  // One will 'cut' useless modes, to keep only the 'visible' ones
  // The basis is reduced to 50 modes
  // 50 to be adjusted during AITs .. may be different for the 4 systems ...
  nuseful = 50;
  DMmodes = DMmodes(,1:nuseful);
  
  // Then one will define a matrix including 62 actuators (DM and tip-tilt),
  // and (nuseful+2) modes. The last "+2" modes are the tip-tilt of the tip-tilt gimbal mount.
  M2V = array(0.0, 62, nuseful+2);
  M2V(1:60, 1:nuseful) = DMmodes;
  M2V(61,nuseful+1) = 1.0;
  M2V(62,nuseful+2) = 1.0;

  
  // Gather the 2 interaction matrices miaHO and miaTT into a single matrix mia (tilts at the end)
  mia = array(0.0, ns, na+2);
  mia(,1:na) = miaHO;
  mia(,na+1:na+2) = miaTT;
  
  // Now computing command matrix
  mim = mia(,+) * M2V(+,);   // this is the Modal Interaction Matrix
  S2M = gral_generalized_inverse(mim,0);     // this is the Modal Control Matrix
}









func gral_pseudoOpenLoopReconstruction(slopesdata, voltsdata, miaHO, miaTT, Fs, latency)
/* DOCUMENT openloopSlopes = gral_pseudoOpenLoopReconstruction(slopesdata, voltsdata, miaHO, miaTT, Fs, latency);

     Input arguments :
   <slopesdata> : set of N frames of slopes (i.e. centroids), (120xN)
   <voltsdata>  : set of N frames of volts (i.e. DM+TT commands) synchronized with slopes, (62xN)
   <miaHO>      : measured interaction matrix (120x60, i.e. 120 lines, 60 columns)
   <miaTT>      : measured interaction matrix (120x2, i.e. 120 lines, 2 columns)
   <Fs>         : sampling frequency (Hertz)
   <latency>    : latency, i.e. time between the start of WFS integration and time when the
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
   
   SEE ALSO:
 */
{
  // get number of slopes <ns> and number of actuators of DM <na>
  ns = dimsof(miaHO)(2);   // number of slopes measured by the Shack-Hartmann, ns should be 120
  na = dimsof(miaHO)(3);   // number of actuators of HO DM, na should be 60
  
  // Gather the 2 interaction matrices miaHO and miaTT into a single matrix mia (tilts at the end)
  mia = array(0.0, ns, na+2);
  mia(,1:na) = miaHO;
  mia(,na+1:na+2) = miaTT;
  
  // determines how many frames are recorded in the data set,
  // will be used later (could also be passed as an argument..)
  nrec = numberof(slopesdata(1,));
  
  // slopes that were "added" by the action of the DM
  sv = mia(,+) * voltsdata(+,);

  // those slopes <sv> need now to be subtracted from the residual slopes
  // with the proper delay.
  // But first one needs to process latency a bit.
  
  // latency expressed in frames (i.e. about 1.35 frames for gravity).
  latency_frame = Fs * latency;                      
  // latency needs to be split into integer and fractional part,
  // i.e. 1.35 = 1 + 0.35
  i_latency_frame = long(latency_frame);              // integer part, i.e. 1
  f_latency_frame = latency_frame - i_latency_frame;  // franctional part, i.e. 0.35

  // now <sv> will be shifted in time
  svs = sv;   // just to allocate memory space same size as <sv>;
  for(i=1; i<=nrec; i++) {
    j = i-i_latency_frame;  // 1 frame before
    k = j-1;                // 2 frames before ...
    if(j<1) j=1;            // except we don't have data before start of file !..
    if(k<1) k=1;            // idem
    svs(,i) = sv(,j) * (1-f_latency_frame) + sv(,k) * f_latency_frame;
  }

  // subtraction of time-shifted DM action from residual slopes data
  openloopSlopes = slopesdata - svs;
  
  return openloopSlopes;
}




func gral_hcor(freq, Fs, latency, G, BP)
/* DOCUMENT   H = gral_hcor(freq,Fs,latency,G,BP);

   Input arguments :
   <freq> is a 1D array of frequencies (usually 1024 points ranging from Fs/2048 to Fs/2).
   <Fs> is the sampling frequency
   <latency> is the latency, in seconds, between the BEGINNING of the integration and the
             start of the command.
   <G> is a scalar, it's the loop gain
   <BP> is a scalar, the cutoff frequency of the DM (seen as a 1st order filter)

   On output, retunrs the square modulus of the correction transfer function of the system
   
*/
{
  Te=1./Fs;                         // sampling period
  p = 1i*2*pi*freq;                 // Laplace variable
  Hint = 1./(1-exp(-p*Te));         // numeric integrator
  Hccd = (1.-exp(-p*Te))/(p*Te);    // zero-order hold with 1/2 frame delay
  Hdac = Hccd;                      // well, same.
  tdelay = latency - Te;            // time between END of the integratino and start of command
  Hret = exp(-p*tdelay);            // latency transfer function
  // transfer func of the DM, as a 1st order filter
  Hmir = 1./(1. + 1i*freq/BP);
  // open-loop transfer function
  Hbo = Hint * Hccd * Hdac * Hret * Hmir;
  // correction transfer function
  Hcor   = 1./abs(1+Hbo*G)^2;
  return Hcor;
}






func gral_modalControlOptimization( slopesdata, voltsdata, miaHO, miaTT, S2M, M2V, gmax, Fs, latency, BP )
/* DOCUMENT optiGain = gral_modalControlOptimization( slopesdata, voltsdata, miaHO, miaTT, S2M, M2V, gmax, Fs, latency, BP )
     
     Input arguments :
   <slopesdata> : set of N frames of slopes (i.e. centroids), (120xN)
   <voltsdata>  : set of N frames of volts (i.e. DM commands) synchronized with slopes, (62xN)
   <miaHO>      : measured interaction matrix (120x60, i.e. 120 lines, 60 columns)
   <miaTT>      : measured interaction matrix (120x2, i.e. 120 lines, 2 columns)
   <S2M>        : modal control matrix
   <M2V>        : matrix of system modes
   <gmax>       : scalar floating point value, maximum gain.
   <Fs>         : sampling frequency (Hertz)
   <latency>    : latency, in seconds, i.e. time between the start of WFS integration and time
                  when the actuator reaches 50% of its command.
   <BP>         : scalar, cutoff frequency of the DM (seen as a 1st order filter)
                  
     Returned value :
   Array of optimal gains on the 50 modes.
   
   SEE ALSO:
 */
{
  // determines how many frames are recorded in the data set,
  // will be used later (could also be passed as an argument..)
  nrec = numberof(slopesdata(1,));
  
  // Reconstruction of pseudo open-loop slopes
  slopes = gral_pseudoOpenLoopReconstruction(slopesdata, voltsdata, miaHO, miaTT, Fs, latency);

  // conversion of slopes into modal coefficients
  modes = S2M(,+) * slopes(+,);
  nmod = numberof(modes(,1));   // just to have the number of modes here

  // Produces an array of gains
  ngain = 15;       // value to be adjusted during AITs
  gmin = 0.0;       // TBC, maybe we'll keep gmin=0.001 to allow for static aberration compensation, if any. TBC during AIT.
  G = span(sqrt(gmin), sqrt(gmax), ngain)^2;   // 1D array from gmin to gmax in ngain points
  
  // Fourier transform of modes coefficients.
  // Compute the square modulus of Fourier transform along time, for each mode.
  // Then multiplies by transfer function for different gains, sum everything
  // over frequency to get the error, and selects the minimum value.
  
  npfft = nrec/2;   // number of useful points in the FFT
  // create a 1D array of frequency ranging from Fs/nrec to Fs/2.0 in npfft points
  freq = span(Fs/nrec, Fs/2.0, npfft);
  optimumGain = array(0.0, nmod);  // memory alloc for 1D array of optimum gains
  for(i=1; i<=nmod; i++) {
    // square modulus of Fast Fourier Transform
    tmp = abs(fft( modes(i,) ))^2;
    // take only half of the points, and reject 1st point
    fftmodes = (2.0/nrec) * tmp(2:npfft+1);

    for(j=1; j<=ngain; j++) {
      phaseError = sum(gral_hcor(freq, Fs, latency, G(j), BP) * fftmodes);

      // Initializes things at startup
      if(j==1) {
        minPhaseError = phaseError;
        jmin = j;
      }
      
      // Search for the minimum value, and the associated gain
      if( phaseError<minPhaseError ) {
        minPhaseError = phaseError;
        jmin = j;
      }
    }
    optimumGain(i) = G(jmin);
  }
  return optimumGain;
}





func gral_computeSystemControlMatrix( S2M, M2V, TT2HO )
/* DOCUMENT mc = gral_computeSystemControlMatrix( S2M, M2V, TT2HO )
     
     Input parameters :
     <S2M>      is the modal control matrix as computed by function computeModalControlMatrix()
     <M2V>      is the matrix of the modes computed by the function gral_createModalBasis()
     <TT2HO>    is the matrix of the DM-tilt modes (60,2) computed by the same function
                than <M2V>.

     Returned parameter : the command matrix of the system.
     
   SEE ALSO:
 */
{
  // matrix multiply  modes * S2M
  // In a "normal system",
  mc = M2V(,+) * S2M(+,);

  // transforming the tip-tilt correction matrix into a DM-tilt correction matrix
  tiltHO = TT2HO(,+) * mc(61:62,)(+,);
  // adding this DM-tilt correction matrix to the DM one
  mc(1:60,) += tiltHO;
  
  return mc;
}




func gral_computeSystemOptimizedControlMatrix( S2M, gainvector, M2V, TT2HO )
/* DOCUMENT mc = gral_computeSystemOptimizedControlMatrix( S2M, gainvector, M2V, TT2HO )
     
     Input parameters :
     <S2M>        is the modal control matrix as computed by function computeModalControlMatrix()
     <gainvector> a 1D vector of modal gains, supposed to be provided by function gral_modalControlOptimization()
     <M2V>        is the matrix of the modes computed by the function gral_createModalBasis()
     <TT2HO>      is the matrix of the DM-tilt modes (60,2) computed by the same function
                  than <M2V>.

     Returned parameter : the command matrix of the system.
     
   SEE ALSO:
 */
{
  // let's allocate memory space and copy the data from the S2M
  tmp = S2M;
  // let's multiply by the gains (dimensions of gainvector and S2M must comply)
  ngain = numberof(gainvector);
  for(i=1; i<=ngain; i++) {
    tmp(i,) *= gainvector(i);
  }
  // matrix multiply  modes * (gains*S2M)
  mc = M2V(,+) * tmp(+,);

  
  // Now we need to take into account the tiptilt handling operated by
  // the CTR component of SpartaLight. We need to add to the DM part
  // of the matrix the compensation of tiptilt by the DM actuators.
  // transformation of the tip-tilt correction matrix into a DM-tilt correction matrix
  tiltHO = TT2HO(,+) * mc(61:62,)(+,);
  // adding this DM-tilt correction matrix to the DM one
  mc(1:60,) += tiltHO;
  
  return mc;
}



