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
import glob
import pandas as pd
import pickle
CalcSFRs = reload(CalcSFRs)
import MeasureFluxes
MeasureFluxes = reload(MeasureFluxes)
import ReturnImfitPars
reload(ReturnImfitPars)
import GetGalaxyList
reload(GetGalaxyList)
import WISE
reload(WISE)
import Templates
reload(Templates)

#######################################################################
# parameters

# Toggle whether to use aperture photometry or to retrieve imfit results for flux measurements
get_imfits = True

# Name of data table containing source list with RAs, Decs, IR SFRs, etc. Might need to delete columns by hand prior
table_name = 'VLAsample.csv'

# test_temps = False

wise_colors = True

template_colors = True
#######################################################################

start_dir = os.getcwd()

# Read in data table given to me by collaboration
t = Table.read(table_name)

# Append columns to table for new data
a = np.empty(len(t))
a[:] = np.nan
b = np.zeros(len(t))
t['data'], t['21 cm Flux'], t['21 cm Flux Error'], t['Luminosity'], t['Luminosity Error (stat.)'], \
t['21 cm SFR'], t['21 cm SFR Error (stat.)'],  t['RMS'], t['w1-w2'], \
t['w3-w4'], t['w2-w3'], t['MySFR'], t['MySFR Err'], \
t['q'], t['W3'], t['W3_err'], t['W4'], t['W4_err'] = b, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a

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
t['W3'].unit = 'Jy'
t['W4'].unit = 'Jy'
t['W3_err'].unit = 'Jy'
t['W4_err'].unit = 'Jy'


names = GetGalaxyList.return_galaxy_list()
templateSFRratio = []

# Go to each galaxy, call photometry and calculate SFRs scripts, and add the outputs to the table
for name in names:

	# If the name exists in my directory, set data flag to true, get index in table
	idx = np.where(t['Name'] == name)[0]
	t['data'][idx] = True


	"""if test_temps:
		#short_name = name[:5]
		#wisefluxes = WISE.mag_to_flux(short_name)
		#twelve_lum = (CalcSFRs.calc_params(wisefluxes[0], 0, t['Z'][idx], 0)[0])/(10**7)
		#twent_lum = (CalcSFRs.calc_params(wisefluxes[1], 0, t['Z'][idx], 0)[0])/(10**7)
		lum_syst = Templates.test_templates(t['Z'][idx])
		#t['IR Error syst.'][idx] = lum_syst
		#t['IR Error Chary'][idx] = Templates.chary_elbaz(t['Z'][idx])"""

	WISEfluxes = WISE.mag_to_flux(name[:5])
	t['W3'][idx] = round(WISEfluxes[0], 6)
	t['W4'][idx] = round(WISEfluxes[1], 5)
	t['W3_err'][idx] = round(WISEfluxes[2], 7)
	t['W4_err'][idx] = round(WISEfluxes[3], 6)

	if wise_colors:
		colors = WISE.colors(name[:5])
	t['w1-w2'][idx] = colors[0]
	t['w3-w4'][idx] = colors[1]
	t['w2-w3'][idx] = colors[2]

	# Call photometry script, returning flux and error, as well as other parameters
	flux_measured = MeasureFluxes.photometry(name, True)

	Flux = 0
	Flux_error = 0
	img_rms = flux_measured[2]

	t['RMS'][idx] = round(img_rms, 8)

	# Retrieve parameters derived by imfit to compare to our estimates
	if get_imfits:
		try:
			imfitpars = ReturnImfitPars.get_imfit(name)
			Flux = float(imfitpars[0])
			Flux_error = float(imfitpars[1])

			t['imfit max'][idx] = imfitpars[2]
			t['imfit max_err'][idx] = imfitpars[3]
			t['aper flux err/imfit err'][idx] = float(flux_measured[1]) / float(imfitpars[1])
			t['geo_mean/rms'][idx] = np.sqrt(Flux * float(imfitpars[2])) / img_rms

			if t['geo_mean/rms'][idx] > 3 and Flux > 0:
				t['detect'][idx] = True
			else:
				t['detect'][idx] = False


			if t['detect'][idx]:
				t['21 cm Flux'][idx] = round(Flux, 6)
				t['21 cm Flux Error'][idx] = round(Flux_error, 7)
			else:
				t['21 cm Flux'][idx] = round(3*img_rms, 6)
				t['21 cm Flux Error'][idx] = np.nan

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
		params_measured = CalcSFRs.calc_params(Flux, Flux_error, t['Z'][idx], t['Z_err'][idx])
		os.chdir(name)
		with open('text/detect.txt', 'w') as f:
			f.write('%s' % 1)
		os.chdir('..')

	# If galaxy is not detected, calculate a SFR and luminosity with 3x the image rms to give an upper limit
	else:
		params_measured = CalcSFRs.calc_params(3*img_rms, np.nan, t['Z'][idx], t['Z_err'][idx])
		os.chdir(name)
		with open('text/detect.txt', 'w') as f:
			f.write('%s' % 0)
		os.chdir('..')

	# Enter calculated SFRs and luminosities to table
	lumin = float(params_measured[0])
	t['Luminosity'][idx] = round(params_measured[0], -28)
	t['Luminosity Error (stat.)'][idx] = round(params_measured[1], -28)
	t['21 cm SFR'][idx] = round(params_measured[2], 1)
	t['21 cm SFR Error (stat.)'][idx] = round(params_measured[3], 1)

	newSFR = Templates.IR_SFRs(t['Z'][idx], name[:5])
	t['MySFR'][idx] = round(float(newSFR[0]), 1)
	t['MySFR Err'][idx] = round(float(newSFR[1]), 1)
	sfr_to_lum = 1/(3.88e-37)
	irlum = float(newSFR[0])*sfr_to_lum
	# extrapolates our 1.519GHz observation to 1.4 GHz using alphaNT = -0.8
	lum_14 = lumin*1.0674
	q = np.log10((irlum/3.75e12)/(lum_14/1.e7))
	t['q'][idx] = round(q, 3)







# Filter table to contain only sources with associated data presently
t_data = t[np.where(t['data'])[0]]

print(t_data['Name', 'Z', 'IR SFR', 'IR SFR Err', '21 cm Flux', '21 cm Flux Error', 'Luminosity', 'Luminosity Error (stat.)', '21 cm SFR', '21 cm SFR Error (stat.)'])

t_non_agn = t_data[np.where(t_data['21 cm SFR'] < 1000.)[0]]
print(np.median(t_non_agn['21 cm SFR']))
print(np.mean(t_non_agn['21 cm SFR']))

os.chdir(start_dir)
t_data.write('table.csv', format='csv', overwrite=True)
# write out table in tex format
t_data['Name', 'Z', '21 cm Flux', '21 cm Flux Error', 'Luminosity', 'Luminosity Error (stat.)', '21 cm SFR', '21 cm SFR Error (stat.)', 'W3', 'W3_err', 'W4', 'W4_err', 'MySFR', 'MySFR Err'].write('textable', format='aastex', overwrite=True)



####################################################
# make observation stats table
####################################################
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
	
