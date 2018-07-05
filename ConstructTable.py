#
# Author: Grayson Petter
# Date: 6/29/18
# Description: Code that calls on other functions to do aperture photometry and calculate a star formation rate
# for each galaxy, compiling all of the relevant data into an astropy table ready for science
#


import os
import numpy as np
from astropy.table import Table
import CalcSFRs
CalcSFRs = reload(CalcSFRs)
import MeasureFluxes
MeasureFluxes = reload(MeasureFluxes)
import SigFigs
reload(SigFigs)
import ReturnImfitPars
reload(ReturnImfitPars)

sig_fig_toggle = True
get_imfits = True


# Read in data table given to me by collaboration
t = Table.read('VLAsample.csv')

# Append columns to table for new data
a = np.empty(len(t))
a[:] = np.nan
b = np.zeros(len(t))
t['data'], t['21 cm Flux'], t['21 cm Flux Error'], t['Luminosity'], t['Luminosity Error (stat.)'], \
t['21 cm SFR'], t['21 cm SFR Error (stat.)'], t['detect_aper'], t['detect_pix'], t['RMS'], \
t['Npixperbeam'], t['Nbeams'], t['MaxValue'], t['Max/noise'], \
t['Flux/error'] = b, a, a, a, a, a, a, a, a, a, a, a, a, a, a

if get_imfits:
    t['imfit flux'], t['imfit flux_err'], t['imfit max'], t['imfit max_err'], t['ratio'] = a, a, a, a, a

# Defining units of each column
t['RA (J2000)'].unit = 'deg'
t['Dec (J2000)'].unit = 'deg'
t['IR SFR'].unit = 'solMass/yr'
t['21 cm Flux'].unit = 'Jy'
t['21 cm Flux Error'].unit = 'Jy'
t['Luminosity'].unit = 'erg/s'
t['Luminosity Error (stat.)'].unit = 'erg/s'
t['21 cm SFR'].unit = 'solMass/yr'
t['21 cm SFR Error (stat.)'].unit = 'solMass/yr'
t['RMS'].unit = 'Jy/beam'
t['MaxValue'].unit = 'Jy/beam'

os.chdir('/users/gpetter/DATA')
names = os.listdir('data_v1')
os.chdir('data_v1')


# Go to each galaxy, call photometry and calculate SFRs scripts, and add the outputs to the table
for name in names:

    idx = np.where(t['Name'] == name)[0]
    t['data'][idx] = True

    flux_measured = MeasureFluxes.photometry(name)

    sig_fig_fluxes = SigFigs.sig_figs(False, 2, flux_measured[0], flux_measured[1])

    if sig_fig_toggle:
        t['21 cm Flux'][idx] = sig_fig_fluxes[0]
        t['21 cm Flux Error'][idx] = sig_fig_fluxes[1]
        t['RMS'][idx] = '%f' % float(('%.' + '%sg' % 2) % float(flux_measured[2]))
    else:
        t['21 cm Flux'][idx] = flux_measured[0]
        t['21 cm Flux Error'][idx] = flux_measured[1]
        t['RMS'][idx] = flux_measured[2]

    t['Npixperbeam'][idx] = flux_measured[3]
    t['Nbeams'][idx] = flux_measured[4]
    t['MaxValue'][idx] = flux_measured[5]
    t['Max/noise'][idx] = flux_measured[6]
    t['Flux/error'][idx] = flux_measured[7]

    if flux_measured[0] > 3*flux_measured[1]:
        t['detect_aper'][idx] = True
    else: t['detect_aper'][idx] = False

    if flux_measured[6] > 3:
        t['detect_pix'][idx] = True
    else: t['detect_pix'][idx] = False

    params_measured = CalcSFRs.calc_params(flux_measured[0], flux_measured[1], t['Z'][idx], 0.001)
    if sig_fig_toggle:
        Lum_sig_fig = SigFigs.sig_figs(True, 2, params_measured[0], params_measured[1])
        t['Luminosity'][idx] = Lum_sig_fig[0]
        t['Luminosity Error (stat.)'][idx] = Lum_sig_fig[1]
        SFR_sig_fig = SigFigs.sig_figs(False, 2, params_measured[2], params_measured[3])
        t['21 cm SFR'][idx] = SFR_sig_fig[0]
        t['21 cm SFR Error (stat.)'][idx] = SFR_sig_fig[1]
    else:
        t['Luminosity'][idx] = params_measured[0]
        t['Luminosity Error (stat.)'][idx] = params_measured[1]
        t['21 cm SFR'][idx] = params_measured[2]
        t['21 cm SFR Error (stat.)'][idx] = params_measured[3]

    if get_imfits:
        try:
            imfitpars = ReturnImfitPars.get_imfit(name)
            t['imfit flux'][idx] = imfitpars[0]
            t['imfit flux_err'][idx] = imfitpars[1]
            t['imfit max'][idx] = imfitpars[2]
            t['imfit max_err'][idx] = imfitpars[3]
            t['ratio'][idx] = float(imfitpars[1])/float(flux_measured[1])
        except:
            os.chdir('..')



# Filter table to contain only sources with associated data presently
t_data = t[np.where(t['data'])[0]]

# Filter table to contain only detections
t_detect = t_data[np.where(np.multiply(t_data['detect_pix'], t_data['detect_aper']))[0]]

print(t_data)

os.chdir('/users/gpetter/PycharmProjects/untitled')
t_detect.write('detected_table.csv', format='csv', overwrite=True)
t_data.write('table.csv', format='csv', overwrite=True)
