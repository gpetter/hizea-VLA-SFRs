#
# Author: Grayson Petter
# Date: 6/29/18
# Description: Code that calls on other functions to do aperture photometry and calculate a star formation rate
# for each galaxy, compiling all of the relevant data into an astropy table ready for science
#


import os
import numpy as np
from astropy.table import Table, Column
from astropy.io import fits
import CalcSFRs
CalcSFRs = reload(CalcSFRs)
import MeasureFluxes
MeasureFluxes = reload(MeasureFluxes)
import SigFigs
reload(SigFigs)
import ReturnImfitPars
reload(ReturnImfitPars)
import GetGalaxyList
reload(GetGalaxyList)
#import WISE
#reload(WISE)
import Templates
reload(Templates)

#######################################################################
# parameters
# Use sig figs when displaying fluxes, luminosities, and SFRs in the table (calculations are still done with floats)
sig_fig_toggle = True

# Toggle whether to use aperture photometry or to retrieve imfit results for flux measurements
get_imfits = True

# Name of data table containing source list with RAs, Decs, IR SFRs, etc. Might need to delete columns by hand prior
table_name = 'VLAsample.csv'

test_temps = True
#######################################################################

start_dir = os.getcwd()

# Read in data table given to me by collaboration
t = Table.read(table_name)

# Append columns to table for new data
a = np.empty(len(t))
a[:] = np.nan
b = np.zeros(len(t))
t['data'], t['21 cm Flux'], t['21 cm Flux Error'], t['Luminosity'], t['Luminosity Error (stat.)'], \
t['21 cm SFR'], t['21 cm SFR Error (stat.)'],  t['RMS'], t['IR Error syst.'] = b, a, a, a, a, a, a, a, a

# Optional toggle to retrieve data from imfit logs
if get_imfits:
    t['imfit max'], t['imfit max_err'], t['aper flux err/imfit err'], t['geo_mean/rms'], t['detect'] = a, a, a, a, a
    t['imfit max'].unit = 'Jy/beam'
    t['imfit max_err'].unit = 'Jy/beam'
else:
    t['detect_aper'], t['detect_pix'], t['Max/noise aper'], t['Flux/error aper'], \
    t['Npixperbeam'], t['Nbeams'], t['MaxValue aper'] = a, a, a, a, a, a, a
    t['MaxValue aper'].unit = 'Jy/beam'

# Defining units of each column
t['RA (J2000)'].unit = 'deg'
t['Dec (J2000)'].unit = 'deg'
t['IR SFR'].unit = 'solMass/yr'
t['IR Luminosity'].unit = 'solLum'
t['21 cm Flux'].unit = 'Jy'
t['21 cm Flux Error'].unit = 'Jy'
t['Luminosity'].unit = 'erg/(s*Hz)'
t['Luminosity Error (stat.)'].unit = 'erg/(s*Hz)'
t['21 cm SFR'].unit = 'solMass/yr'
t['21 cm SFR Error (stat.)'].unit = 'solMass/yr'
t['RMS'].unit = 'Jy/beam'


names = GetGalaxyList.return_galaxy_list()

# Go to each galaxy, call photometry and calculate SFRs scripts, and add the outputs to the table
for name in names:

    

    # If the name exists in my directory, set data flag to true, get index in table
    idx = np.where(t['Name'] == name)[0]
    t['data'][idx] = True


    if test_temps:
	    #short_name = name[:5]
	    #wisefluxes = WISE.mag_to_flux(short_name)
	    #twelve_lum = (CalcSFRs.calc_params(wisefluxes[0], 0, t['Z'][idx], 0)[0])/(10**7)
	    #twent_lum = (CalcSFRs.calc_params(wisefluxes[1], 0, t['Z'][idx], 0)[0])/(10**7)
	    lum_syst = Templates.test_templates(t['Z'][idx])
	    t['IR Error syst.'][idx] = lum_syst

    # Call photometry script, returning flux and error, as well as other parameters
    flux_measured = MeasureFluxes.photometry(name, True)

    Flux = 0
    Flux_error = 0
    img_rms = flux_measured[2]

    if sig_fig_toggle:
        t['RMS'][idx] = '%f' % float(('%.' + '%sg' % 2) % float(img_rms))
    else:
        t['RMS'][idx] = img_rms

    





    # Retrieve parameters derived by imfit to compare to our estimates
    if get_imfits:
        try:
            imfitpars = ReturnImfitPars.get_imfit(name)
            Flux = float(imfitpars[0])
            Flux_error = float(imfitpars[1])
            if sig_fig_toggle:
                sig_im_flux = SigFigs.sig_figs(False, 2, Flux, Flux_error)
                t['21 cm Flux'][idx] = sig_im_flux[0]
                t['21 cm Flux Error'][idx] = sig_im_flux[1]
            else:
                t['21 cm Flux'][idx] = Flux
                t['21 cm Flux Error'][idx] = Flux_error
            t['imfit max'][idx] = imfitpars[2]
            t['imfit max_err'][idx] = imfitpars[3]
            t['aper flux err/imfit err'][idx] = float(flux_measured[1]) / float(imfitpars[1])
            t['geo_mean/rms'][idx] = np.sqrt(Flux * float(imfitpars[2])) / float(t['RMS'][idx])
            if t['geo_mean/rms'][idx] > 3 and t['21 cm Flux'][idx] > 0:
                t['detect'][idx] = True
            else:
                t['detect'][idx] = False
        except:
            print('failed')
            os.chdir('..')
    else:
        # Set elements in table equal to results from photometry. Apply sig figs if desired
        Flux = flux_measured[0]
        Flux_error = flux_measured[1]
        t['Npixperbeam'][idx] = flux_measured[3]
        t['Nbeams'][idx] = flux_measured[4]
        t['MaxValue aper'][idx] = flux_measured[5]
        if sig_fig_toggle:
            sig_fig_fluxes = SigFigs.sig_figs(False, 2, Flux, Flux_error)
            t['21 cm Flux'][idx] = sig_fig_fluxes[0]
            t['21 cm Flux Error'][idx] = sig_fig_fluxes[1]

        else:
            t['21 cm Flux'][idx] = Flux
            t['21 cm Flux Error'][idx] = Flux_error


    with open((name + '/text/max.txt'), 'w') as f_max:
        f_max.write('%s' % float(flux_measured[5]))



    if not get_imfits:

        t['Max/noise aper'][idx] = flux_measured[6]
        t['Flux/error aper'][idx] = flux_measured[7]

        # If flux is a 3 sigma result, count as detection
        if flux_measured[0] > 3*flux_measured[1]:
            t['detect_aper'][idx] = True
        else: t['detect_aper'][idx] = False

        # If maximum pixel value is a 3 sigma result, count as detection
        if flux_measured[6] > 3:
            t['detect_pix'][idx] = True
        else: t['detect_pix'][idx] = False
        detection = t['detect_aper'][idx] and t['detect_pix'][idx]
    else:
        detection = t['detect'][idx]



    # If both previous detection criteria are met, the galaxy is deemed a detection, calculate a SFR and Luminosity
    if detection:
        params_measured = CalcSFRs.calc_params(Flux, Flux_error, t['Z'][idx], 0.001)
        os.chdir(name)
        with open('text/detect.txt', 'w') as f:
            f.write('%s' % 1)
        os.chdir('..')

    # If galaxy is not detected, calculate a SFR and luminosity with 3x the image rms to give an upper limit
    else:
        params_measured = CalcSFRs.calc_params(3*img_rms, img_rms, t['Z'][idx], 0.001)
        os.chdir(name)
        with open('text/detect.txt', 'w') as f:
            f.write('%s' % 0)
        os.chdir('..')


    # Enter calculated SFRs and luminosities to table, with sig figs if desired
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

    



# Filter table to contain only sources with associated data presently
t_data = t[np.where(t['data'])[0]]


print(t_data['Name', 'Z', 'IR SFR', 'IR SFR Err', '21 cm Flux', '21 cm Flux Error', 'Luminosity', 'Luminosity Error (stat.)', '21 cm SFR', '21 cm SFR Error (stat.)'])


os.chdir(start_dir)
t_data.write('table.csv', format='csv', overwrite=True)
t_data['Name', 'Z', 'IR SFR', 'IR SFR Err', '21 cm Flux', '21 cm Flux Error', 'Luminosity', 'Luminosity Error (stat.)', '21 cm SFR', '21 cm SFR Error (stat.)'].write('textable', format='aastex', overwrite=True)

names = GetGalaxyList.return_galaxy_list()

t_obs = Table([t_data['Name']])

c = np.zeros(len(t_obs))

t_obs['rms'], t_obs['Bmaj'], t_obs['Bmin'], t_obs['PA'] = c, c, c, c

t_obs['rms'].unit = 'uJy/beam'
t_obs['Bmaj'].unit = 'arcsec'
t_obs['Bmin'].unit = 'arcsec'
t_obs['PA'].unit = 'deg'


for name in names:
	
	os.chdir(name)
	idx = np.where(t_obs['Name'] == name)[0]
	with open('text/stdev.txt', 'r') as f:
		std_dev = float(f.readline())  # Jy/beam

	imgname = '%s.cutout.pbcor.fits' % name
	hdu = fits.open(imgname)
	bmaj = hdu[0].header['bmaj']
	bmin = hdu[0].header['bmin']
	angle = hdu[0].header['bpa']
	
	if sig_fig_toggle:
		std_dev = round((std_dev*(10**6)), 1)
		bmaj = round((bmaj*3600), 1)
		bmin = round((bmin*3600), 1)
		angle = round(angle, 1)		



	t_obs['rms'][idx] = std_dev
	t_obs['Bmaj'][idx] = bmaj
	t_obs['Bmin'][idx] = bmin
	t_obs['PA'][idx] = angle

	os.chdir('..')
	
os.chdir(start_dir)
t_obs.write('obs_table.csv', format='csv', overwrite=True)
t_obs.write('obs_tex', format='aastex', overwrite=True)
	
