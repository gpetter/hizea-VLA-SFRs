#
# Author: Grayson Petter


import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astropy.modeling import models
from astropy import units as u
import math
#from specutils.spectra import Spectrum1D
#from specutils.fitting import fit_generic_continuum

import pandas as pd
import glob

c = 299792458.

# bins of 1 GHz
d_nu = 1000000000.


templates = glob.glob('Comprehensive_library/*.txt')




def redshift_spectrum(z, template):

	t = pd.read_csv(template, sep='   ', engine='python')
	t.columns = ['lambda_e', 'L', 'dL']

	wavelengths = np.array(t.iloc[:, 0])
	Lums = np.array(t.iloc[:, 1])

	shifted_len = []
	shifted_Lum = []

	for y_emit in wavelengths:
		shifted_len.append(float(y_emit)*(1+z))
	#for l_emit in Lums:
		#shifted_Lum.append(float(l_emit)/(1+z))
	shifted_len = np.array(shifted_len)
	#shifted_Lum = np.array(shifted_Lum)
	twelve_mu = (np.abs(shifted_len - 12)).argmin()
	twenty_two_mu = (np.abs(shifted_len -22)).argmin()

	return(shifted_len, Lums, Lums[twelve_mu], Lums[twenty_two_mu])



def interpolate_spec(shifted_spec):

	# convert wavelengths to frequencies
	nus = (10**6)*c/(shifted_spec[0])

	# reverse lists so frequencies go from low to high for convenience
	reversed_nus = np.flipud(nus)
	reversed_lums = np.flipud(shifted_spec[1])

	# find smallest factor of d_nu Hz greater than the smallest frequency in the list
	start = (reversed_nus[0]+int(d_nu))-(reversed_nus[0] % int(d_nu))

	# number of d_nu Hz intervals in entire template
	span = reversed_nus[len(reversed_nus)-1]-reversed_nus[0]
	chunks = int(math.floor(span/d_nu))

	new_nus, new_lums = [], []
	current_nu = start

	# linearly interpolate to frequencies in d_nu Hz steps
	for x in range(chunks):
		new_nus.append(current_nu)
		new_lums.append(np.interp(current_nu, reversed_nus, reversed_lums))
		current_nu += d_nu
	return(new_nus, new_lums)



# integrate spectrum using trapezoid method (rectangle below 
def integrate_spectrum(freqs, Ls):
	tot_ir = 0.
	for x in range(len(freqs)-1):
		if Ls[x+1] > Ls[x]:
			rect = d_nu*Ls[x]
			tri = (d_nu*(Ls[x+1]-Ls[x]))/2.
			
		else:
			rect = d_nu*Ls[x+1]
			tri = (d_nu*(Ls[x]-Ls[x+1]))/2.

		tot_ir = tot_ir + (rect+tri)


	return(tot_ir)

def test_templates(zz):
	ratio_list = []
	for x in range(len(templates)):	
	

		shifted_spectrum = redshift_spectrum(zz, templates[x])
		twenty_two_lum = shifted_spectrum[3]
		interped_spectrum = interpolate_spec(shifted_spectrum)
		total_ir = integrate_spectrum(interped_spectrum[0], interped_spectrum[1])
		ratio_list.append(total_ir/twenty_two_lum)
	ratio_list = np.array(ratio_list)
	stdev = np.std(ratio_list)
	range_ratios = np.ptp(ratio_list)
	print(stdev, range_ratios)








"""plt.clf()
plt.cla()
plt.close()
plt.figure(121)
plt.plot((c/thing[0]), thing[1])
plt.show()
plt.clf()
plt.cla()
plt.close()"""


