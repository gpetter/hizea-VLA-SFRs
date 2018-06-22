import astropy
import os
import numpy as np
from photutils import CircularAperture
from photutils import EllipticalAperture
from photutils import EllipticalAnnulus
from photutils import CircularAnnulus
from photutils import aperture_photometry
from astropy.io import fits
from astropy.table import Table
import MeasureFluxes
import CalcSFRs

t = Table.read('VLAsample.csv')

a = np.empty(len(t))
a[:] = np.nan
t['Flux'], t['Flux_error'], t['Luminosity'], t['Luminosity_error'], t['SFR'], t['SFR_error'] = a, a, a, a, a, a

os.chdir('/users/gpetter/DATA')

names = os.listdir('data_v1')
os.chdir('data_v1')

for name in names:
    idx = np.where(t['Name'] == name)[0]
    flux_measured = MeasureFluxes.photometry(name)
    t['Flux'][idx] = flux_measured[0]
    t['Flux_error'][idx] = flux_measured[1]
    params_measured = CalcSFRs.calc_params(flux_measured[0], flux_measured[1], t['Z'][idx], 0.001)
    t['Luminosity'][idx] = params_measured[0]
    t['Luminosity_error'][idx] = params_measured[1]
    t['SFR'][idx] = params_measured[2]
    t['SFR_error'][idx] = params_measured[3]

print(t)

os.chdir('/users/gpetter/PycharmProjects/untitled')