# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 16:52:46 2014

Mix of various routines around Zernike modes

@author: tristanbuey
"""


# Load libraries
import numpy as np;






def gammX(n, i0):
    """
    gammX(n, i0)
    Computes Noll matrix for derivatives of Zernike Gx.
    The matrix applies on a Zernike vector z, and produces
    the Zernike decomposition z' of x-derivative :
    z' = Gx . z

    n = number of zernike coefficients on the input vector z.
    i0 = zernike index of the first coefficient of z (1=piston, 2=tip, ...)

    It results that Gx is a matrix with a size (n+i0-1, n).
    """
    gg = np.zeros((i0+n-1,n));
    # variable i will span Zernike indexes, starting at i0.
    for i in range (i0, i0+n):
        # variable j spans Zernike indexes, starting at piston
        # and stopping at i
        for j in range(1, i+1):
            gg[j-1,i-i0] = gamX(i,j);
    return gg;




def gammY(n, i0):
    """
    gammY(n, i0)
    Computes Noll matrix for derivatives of Zernike Gy.
    The matrix applies on a Zernike vector z, and produces
    the Zernike decomposition z' of y-derivative :
    z' = Gy . z

    n = number of zernike coefficients on the input vector z.
    i0 = zernike index of the first coefficient of z (1=piston, 2=tip, ...)

    It results that Gy is a matrix with a size (n+i0-1, n).
    """
    gg = np.zeros((i0+n-1,n));
    # variable i will span Zernike indexes, starting at i0.
    for i in range(i0, i0+n):
        # variable j spans Zernike indexes, starting at piston
        # and stopping at i
        for j in range(1, i+1):
            gg[j-1,i-i0] = gamY(i,j);
    return gg;



"""
A lot of sub-functions to calculate the Noll matrix
"""

def pair(number):
    return number % 2 == 0;

def impair(num):
   return num % 2 != 0;



def nm(i):
    """
    For a given Zernike mode of index <i>, returns the
    radial and azimutal orders (n,m)
    """
    n = int( (-1.+np.sqrt(8*(i-1)+1))/2.);
    p = (i-(n*(n+1))/2);
    k = n%2;
    m = int((p+k)/2)*2 - k;
    return (n,m);


def gamY(i,j):
    """
    Input arguments:
    2 scalar int i and j, that are indexes of Zernike modes.

    Returns the coefficient of the derivative matrix of (Noll R.J., 1976)
    The algorithm coded below is a python translation of the series of
    rules that Noll has enounced in his article of 1976, for derivating
    Zernike.
    Warning: Unfortunately Noll had made a little error in his rules,
    that has been corrected in this program.
    """
    # determine radial and azimutal orders of Zernike number i
    ni,mi = nm(i);
    # idem for j
    n,m = nm(j);

    # Noll's rules :
    if(mi==(m-1) or mi==(m+1)):
        if(m==0 or mi==0):
            if((m==0 and impair(i)) or (mi==0 and impair(j))):
                return np.sqrt(2*(n+1)*(ni+1));
            else:
                return 0.00;
        else:
            if(impair(i+j)):
                if((mi==m+1 and impair(j)) or (mi==m-1 and pair(j))):
                    return -np.sqrt((n+1)*(ni+1));
                else:
                    return np.sqrt((n+1)*(ni+1));
            else:
                return 0.0;
    else:
        return 0.0;

    return;



def gamX(i,j):
    """
    Input arguments:
    2 scalar int i and j, that are indexes of Zernike modes.

    Returns the coefficient of the derivative matrix of (Noll R.J., 1976)
    The algorithm coded below is a python translation of the series of
    rules that Noll has enounced in his article of 1976, for derivating
    Zernike.
    Warning: Unfortunately Noll had made a little error in his rules,
    that has been corrected in this program.
    """
    # determine radial and azimutal orders of Zernike number i
    ni,mi = nm(i);
    # idem for j
    n,m = nm(j);

    # Noll's rules :
    if(mi==m-1 or mi==m+1):
        if(m==0 or mi==0):
            if((m==0 and pair(i)) or (mi==0 and pair(j))):
                return np.sqrt(2*(n+1)*(ni+1));
            else:
                return 0.00;
        else:
            if( (j+i)%2==0 ):
                return np.sqrt((n+1)*(ni+1));
            else:
                return 0.00;
    else:
        return 0.0;

    return;






def polyfute(m,n):
    """
    Les coefs des poly de zer sont des K_mn(s).
    Le coeff K_mn(s) pondère r^(n-2s)
    Il y a la relation de recurrence
    K_mn(s+1) =  K_mn(s) * ((n+m)/2-s)*((n-m)/2-s)/(s+1)/(n-s)
    Il y a aussi
    K_mn(0) = n! / ((n+m)/2)! / ((n-m)/2)!
    """

    a = np.zeros(n+1)

    # Calcul de K_mn(0)
    st = 2                           # start index for dividing by ((n-m)/2)!
    coef = 1.00
    for i in range((n+m)/2+1, n+1):
        if( st<=((n-m)/2) and i%st==0 ) :
            j = i/st
            st = st+1
            coef = coef*j
        else:
            coef = coef*i

    # division by ((n-m)/2)! (has already been partially done)
    for i in range(st,(n-m)/2+1):
        coef = coef / i

    a[n] = round(coef);   # pour K_nm(0)

    for i in range(1,(n-m)/2+1):
        coef = -coef * ((n+m)/2-i+1)*((n-m)/2-i+1);
        coef = coef / i;
        coef = coef / (n-i+1);
        a[n-2*i] = round(coef)

    return a






def evaluate_poly(n,m,a,r):
    """
     evaluate_poly(n,m,a,r)
     n is the radial order
     m is the azimutal order
     a[] is the list of coefficient, with a(i+1) the coeff of r^i
     r is the variable of the polynomial
    """
    if n>1 :
        r2 = r*r

    p = a[n]
    for i in range(n-2,m-1,-2):
        p = p*r2 + a[i]

    if(m==0): return p
    elif(m==1): p*=r
    elif(m==2): p*=r2
    else: p = p * r**m

    return p



def zer(r,t,i):
    """
    Computes Zernike polynom of index i, at point (r,t)

    The algo is using
    1) a recursive function to compute coefficients of the polynom, so that
    there is no factorial function of large numbers to invoke (risk of
    roundoff errors, plus takes more exec time)

    2) a smarter way to compute polynomial expressions such as
     ax^3+bx^2+cx+d = x(x(ax+b)+c)+d
    to avoid roundoff errors and minimize number of operations

    """

    if(i==1):
        return np.ones_like(r+t)

    # calcul de n et m a partir de i
    n = int( (-1.+np.sqrt(8*(i-1)+1))/2.)
    p = (i-(n*(n+1))/2);
    k = n%2;
    m = int((p+k)/2)*2 - k;

    a = polyfute(m,n)

    Z = evaluate_poly(n,m,a,r) * np.sqrt(n+1);
    if( m!=0 ):
        Z *= np.sqrt(2);
        if( i%2 ):
            Z *= np.sin(m*t)
        else:
            Z *= np.cos(m*t)
    return Z





def mkxy(npt, center, xy=0):
    """
    Defines a meshgrid of npt X npt points, and express
    it in polar coordinates.
    Returns a tuple (r,theta).

    When center==1, the center of the x,y map is 'Fourier-centered'.
    When optional argument xy==1, the returned tuple is (x, y)
    When optional argument xy==2, the returned tuple is (r, theta, x, y)
    """
    # generate an array of coordinates
    if center==1:
        x = np.linspace(-1,1,npt+1)[0:npt]
    else:
        x = np.linspace(-1,1,npt)
    # Warning: meshgrid swaps the 2 axis x and y.
    # That's why it's required to re-swap them here
    y,x = np.meshgrid(x,x)
    if( xy==1 ):
        return (x,y)
    # generates a map of the distance of subapertures to pupil center
    r = np.sqrt(x**2 + y**2)
    # generates a map of the azimut angle of subapertures
    theta = np.arctan2(y,x)
    if( xy==2 ):
        return r,theta,x,y
    else:
        return r,theta