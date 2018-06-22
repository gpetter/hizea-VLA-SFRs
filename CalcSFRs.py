import numpy as np
from sympy import *


# Calculate a luminosity, star formation rate, and uncertainties given a flux density
# Using equation relating synchrotron emission to star formation rate given in Murphy et. al (2011)
# Also using Condon & Matthews (2018) to calculate luminosity distance
def calc_params(flux, flux_error, redshift, redshift_error):
    # Defining symbols (sympy)
    # redshift
    z = Symbol('z')
    # redshift uncertainty
    zunc = Symbol('zunc')
    # frequency
    nu = Symbol('nu')
    nuunc = Symbol('nuunc')
    # alpha value
    al = Symbol('al')
    alunc = Symbol('alunc')
    # flux
    f = Symbol('f')
    f_unc = Symbol('f_unc')
    # Hubble constant
    Ho = Symbol('Ho')
    Ho_unc = Symbol('Ho_unc')

    # define speed of light and Hubble distance
    c = 299792.458  # km/s
    Dho = c/Ho


    # Define symbolic formluas for desired quantities
    # convenience definition from Murphy et al. paper
    a = 1/(1+z)
    # Comoving distance formula
    Dc = (Dho/(a/(1-a)+0.2278+0.2070*(1-a)/(0.785+a)-0.0158*(1-a)/((0.312+a)**2)))*3.08567758128*(10**24)  # cm
    # luminosity distance formula
    Dl = (1+z)*Dc  # cm
    # spectral luminosity distance
    Dl_nu = Dl*((1+z)**(-(al+1)/2))
    # inverse square law to get luminosity
    Lumform =(4*np.pi*f*(10**-23)*Dl_nu**2)  # ergs/s
    # SFR formula in solar masses/yr (Murphy et. al)
    SFRform = (6.64e-29*(nu**(-al))*Lumform)
    # luminosity uncertainty formula - simple error propagation
    Lum_unc = ((diff(Lumform, z)*zunc)**2+(diff(Lumform, f)*f_unc)**2+(diff(Lumform, Ho)*Ho_unc)**2)**0.5
    # SFR uncertainty formula
    SFRuncertainty = ((diff(SFRform, z)*zunc)**2+(diff(SFRform, nu)*nuunc)**2+
                      (diff(SFRform, al)*alunc)**2+(diff(SFRform, f)*f_unc)**2+(diff(SFRform, Ho)*Ho_unc)**2)**.5

    # Define constants
    Hubble = 70
    Hubble_unc = 2
    freqs = 1.51976491105  # GHz
    freqsigs = 0.01
    alphas = -0.7
    alphasig = 0.05

    output = []

    # substitute values in from arguments to calculate
    output.append(Lumform.subs({f: flux, z: redshift, al: alphas, Ho: Hubble}))
    output.append(Lum_unc.subs({f: flux, z: redshift, f_unc: flux_error, zunc: redshift_error, al: alphasig,
                                Ho: Hubble, Ho_unc: Hubble_unc}))
    output.append(SFRform.subs({nu: freqs, al: alphas, f: flux, z: redshift, Ho: Hubble}))

    output.append(SFRuncertainty.subs({nu: freqs, al: alphas, f: flux, z: redshift,
                zunc: redshift_error, nuunc: freqsigs, alunc: alphasig, f_unc: flux_error, Ho: Hubble,
                Ho_unc: Hubble_unc}))
    return(output)


