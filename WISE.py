#
# Author: Grayson Petter


import os
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
import pandas as pd
import CalcSFRs
CalcSFRs = reload(CalcSFRs)

fc_w_three = [0.9169, 0.9393, 1.0088, 1.1344, 0.9373, 1., 1.1081, 1.2687]
fc_w_four = [0.9905, 0.9934, 1.0013, 1.0142, 0.9926, 1., 1.013, 1.0319]
w_three = 29.045
w_four = 8.284




def mag_to_flux(name):
	
	pat = 'unWISE/%s.fits' % (name)
	t = Table.read(pat)

	print('1-2')
	print(float(t['w1_mag'])-float(t['w2_mag']))
	print('2-3')
	print(float(t['w2_mag'])-float(t['w3_mag']))
	print('3-4')
	print(float(t['w1_mag'])-float(t['w2_mag']))
	
	print('choose alpha')
	alpha = input()
	
	color_corr = [fc_w_three[alpha], fc_w_four[alpha]]

	flux_three = (w_three/color_corr[0])*(10**(-(float(t['w3_mag'])/2.5)))
	flux_four = (w_four/color_corr[1])*(10**(-(float(t['w4_mag'])/2.5)))
	

	return(flux_three, flux_four)

test = mag_to_flux('J1107')
print(test)
print(((CalcSFRs.calc_params(test[0], 0, 0.467, 0))[0])/(10**7))
print(((CalcSFRs.calc_params(test[1], 0, 0.467, 0))[0])/(10**7))
