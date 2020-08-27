#
# Author: Grayson Petter
# Date: 6/28/18
# Description: Given a flux density in Jy and a redshift of a galaxy observed at a frequency nu, this code will
# return a luminosity, star formation rate, and uncertainties.
# This code uses equations given in Murphy et. al (2011), and Condon & Matthews (2018).
#

import numpy as np
from sympy import *

allowthermal = True


def tabatabaei(nu, al, z):
    thermal_frac = 1./(1.+13.*((nu*(1+z)) ** (0.1 + al)))
    return(thermal_frac)

# Calculate a luminosity, star formation rate, and uncertainties given a flux density
# Using equation relating synchrotron emission to star formation rate given in Murphy et. al (2011)
# Also using Condon & Matthews (2018) to calculate spectral luminosity distance
def calc_params(flux, flux_error, redshift, redshift_error):

    # Defining symbols (sympy)
    z = Symbol('z')     # redshift
    zunc = Symbol('zunc')   # redshift uncertainty
    nu = Symbol('nu')   # frequency
    nuunc = Symbol('nuunc')     # frequency uncertainty
    al = Symbol('al')   # non-thermal spectral index alpha
    alunc = Symbol('alunc')     # alpha uncertainty
    f = Symbol('f')     # flux
    f_unc = Symbol('f_unc')     # flux uncertainty
    Ho = Symbol('Ho')   # hubble constant
    Ho_unc = Symbol('Ho_unc')   # hubble uncertainty

    # define speed of light and Hubble distance
    c = 299792.458  # km/s
    Dho = c/Ho

    # Define symbolic formluas for desired quantities
    # convenience definition from Murphy et al. paper
    a = 1/(1+z)
    # Comoving distance formula
    # 3E24 factor converts between cm and Mpc
    Dc = (Dho/(a/(1-a)+0.2278+0.2070*(1-a)/(0.785+a)-0.0158*(1-a)/((0.312+a)**2)))*3.08567758128*(10**24)  # cm
    # luminosity distance formula
    Dl = (1+z)*Dc  # cm
    # spectral luminosity distance
    Dl_nu = Dl*((1+z)**(-(al+1)/2))
    # inverse square law to get luminosity
    Lumform =(4*np.pi*f*(10**-23)*Dl_nu**2)  # ergs/s
    # SFR formula in solar masses/yr (Murphy et. al)
    kroupa_to_salpeter = 1.5
    #SFRform = kroupa_to_salpeter*(6.64e-29*(nu**(-al))*Lumform)
    if allowthermal:
        L_NT = Lumform/(1+1/13.*((1+z)*nu)**(-0.1-al))
        SFRform = 6.64e-29 * (nu ** (-al)) * L_NT
    else:
        SFRform = 6.64e-29 * (nu ** (-al)) * Lumform
    # luminosity uncertainty formula - simple error propagation
    Lum_stat_unc = ((diff(Lumform, f)*f_unc)**2)**0.5
    Lum_syst_unc = ((diff(Lumform, z)*zunc)**2+(diff(Lumform, Ho)*Ho_unc)**2)**0.5
    # SFR uncertainty formula
    SFR_stat_uncertainty = ((diff(SFRform, f)*f_unc)**2)**.5
    SFR_syst_uncertainty = ((diff(SFRform, z) * zunc) ** 2 + (diff(SFRform, Ho) * Ho_unc)**2 + (diff(SFRform, al) * alunc)**2) ** .5

    # Define constants
    Hubble = 70     # km/s/Mpc
    Hubble_unc = 2
    freqs = 1.51976491105  # GHz
    freqsigs = 0
    alphas = -0.8
    alphasig = 0.05

    output = []

    # substitute in values into symbolic expressions
    # SFR values
    SF = SFRform.subs({nu: freqs, al: alphas, f: flux, z: redshift, Ho: Hubble})
    SF_stat = SFR_stat_uncertainty.subs({nu: freqs, al: alphas, f: flux, z: redshift,
                                            f_unc: flux_error, Ho: Hubble})
    SF_syst = SFR_syst_uncertainty.subs({nu: freqs, al: alphas, f: flux, z: redshift, zunc: redshift_error,
                                         alunc: alphasig, f_unc: flux_error, Ho: Hubble, Ho_unc: Hubble_unc})
    # luminosity values
    Lum = Lumform.subs({f: flux, z: redshift, al: alphas, Ho: Hubble})
    Lum_stat = Lum_stat_unc.subs({f: flux, z: redshift, f_unc: flux_error, Ho: Hubble, al: alphas})
    Lum_syst = Lum_syst_unc.subs({f: flux, z: redshift, f_unc: flux_error, zunc: redshift_error,
                                  Ho: Hubble, Ho_unc: Hubble_unc})

    output.append(Lum)
    output.append(Lum_stat)
    output.append(SF)
    output.append(SF_stat)

    return output




