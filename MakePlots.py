import matplotlib.pyplot as plt
import matplotlib.colors as matcol
import matplotlib as mpl
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
from astropy.table import Table
import os
import aplpy
import glob
import pandas as pd
import matplotlib.transforms as tf
from astropy.table import Table
from astropy.io import fits
from astropy.stats import biweight_midvariance
from astropy.cosmology import FlatLambdaCDM
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils import Background2D
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.utils import median_survival_times
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
import scipy.stats as st
import astropy.units as u
from mpl_toolkits.axes_grid1 import make_axes_locatable
import CalcSFRs
reload(CalcSFRs)
import WISE
reload(WISE)
import GetGalaxyList
reload(GetGalaxyList)
import Templates
reload(Templates)

kmf = KaplanMeierFitter(alpha=0.16)

####################################################################################
# parameters

# Keep this true
use_imfit = True
# Use plt.annotate to place the names of the galaxies next to points
label_points = False

marker_size = 10
####################################################################################





# Read in astropy table
t = Table.read('table.csv')
# read in Smolcic 2017 COSMOS data
t_survey_raw = Table.read('smolcic.fit')

if use_imfit:
    correlation = t['detect']
else:
    correlation = np.multiply(t['detect_pix'], t['detect_aper'])

# Table containing detections only
t_detect = t[np.where(correlation == 1)[0]]

# Filter table into 3 separate tables, one for detections, non-detections, and flagged sources (AGN)
t_nondetect = t[np.where(correlation == 0)[0]]
t_ok = t_detect[np.where(t_detect['21 cm SFR'] < 1000)[0]]
t_bad = t_detect[np.where(t_detect['21 cm SFR'] > 1000)[0]]

t_color = t[np.where(~np.isnan(t['w3-w4']))[0]]

t_shift = t_color[np.where(t_color['IR SFR'] > (1.5 * t_color['21 cm SFR']))[0]]
t_else = t_color[np.where((t_color['IR SFR'] < (1.5 * t_color['21 cm SFR'])) & (t_color['21 cm SFR'] < 1000))[0]]

# survey data to plot on top of SFR comparison plot
surv_redshift = np.array(t_survey_raw['zbest'])
surv_rad_lum = np.array(t_survey_raw['logL21cm'])
surv_ir_lum = np.array(t_survey_raw['logLTIRSF'])

# multiplying columns, so we can eliminate data points which don't have all 3 (redshift, radio luminosity, and IR luminosity)
multed = np.multiply(np.multiply(surv_redshift, surv_rad_lum), surv_ir_lum)
t_survey_data = t_survey_raw[np.where(multed > 0)[0]]

# applying redshift cut similar to our sample
surv_z = np.array(t_survey_data['zbest'])
t_survey = t_survey_data[np.where((surv_z > 0.4) & (surv_z < 0.75))]


# Do two weighted linear fits to detections only. Fix the second fit to have a y-intercept of zero
def linear_fit(x_vals, y_vals, x_err, y_err):
    # Do first fit with just y errors
    tmp_fit = np.polyfit(x_vals, y_vals, 1, w=1. / y_err)
    tmp_fit_fn = np.poly1d(tmp_fit)

    # Determine total error from previous fit
    err_tot = np.sqrt((y_err) ** 2 + (tmp_fit_fn[1] * x_err) ** 2)

    # Fit again with total error
    fit, residuals, _, _, _ = np.polyfit(x_vals, y_vals, 1, w=1 / err_tot, full=True)
    fit_fn = np.poly1d(fit)
    print(fit_fn)

    # Define linear function for scipy.curve_fit to fit
    def func(x, a, b):
        return a + b * x

    # Fit while holding y intercept to zero, use y errors as weights
    tmp_popt_cons, _ = curve_fit(func, x_vals, y_vals, bounds=([0, -10], [1e-10, 10]), sigma=y_err)

    # Use fitted function to calculate the total error as a result of x and y errors
    new_err_tot = np.sqrt((y_err) ** 2 + (tmp_popt_cons[1] * x_err) ** 2)

    # Use the total error to fit the data again
    popt_cons, _ = curve_fit(func, x_vals, y_vals, bounds=([0, -10], [1e-10, 10]), sigma=new_err_tot)

    held_to_zero = np.poly1d([popt_cons[1], 0])
    print(held_to_zero)

    """reduced_chi_squared = residuals[0] / (float(len(x_vals) - 2))
    print(reduced_chi_squared)"""

    return [fit_fn, held_to_zero]

sfr_to_use = 'MySFR'
err_to_use = 'MySFR Err'


# main plot for paper
# comparing SFRs derived from radio & IR, plus extra axes for comparing luminosities in radio & IR
def plot_all_SFRs():

    x_axis_lim = 600
    y_axis_lim = 600

    # Making arrays of SFRs and uncertainties for non-AGN detections
    irsfrok = np.array(t_ok[sfr_to_use]).astype(float)
    irok_uncertainty = np.array(t_ok[err_to_use])
    ##############################
    radiosfr_ok = np.array(t_ok['21 cm SFR']).astype(float)
    sfr_ok_uncertainty = np.array(t_ok['21 cm SFR Error (stat.)']).astype(float)
    names = t_ok['Name']
    labels = np.random.randint(0, 3, size=len(irsfrok))



    # and for non-detections
    ir_non_detect = np.array(t_nondetect[sfr_to_use])
    ir_non_unc = np.array(t_nondetect[err_to_use])
    radio_non = np.array(t_nondetect['21 cm SFR'])
    radio_non_unc = np.array(t_nondetect['21 cm SFR Error (stat.)'])
    non_names = t_nondetect['Name']

    # Get indices of sources which lie below proposed detection limit
    # detect_x_lims = (irsfrok < 30)
    # non_detect_x_lims = (ir_non_detect < 30)

    # Generate one-to-one line
    one_to_one = np.poly1d([1, 0])
    # Do the linear fits to non-AGN detections only
    #fits = linear_fit(irsfrok, radiosfr_ok, irok_uncertainty, sfr_ok_uncertainty)

    # Generate subplots for the broken axis effect, ax2 for majority of data and ax for AGN
    #fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 10), dpi=300, gridspec_kw={'height_ratios': [1, 4]})
    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
    plt.style.use('default')
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"




    # Plot all data with different markers. For non-detections, make y-axis upper-limit arrows. For sources below
    # proposed detection limit, make x-axis upper limit arrows.

    ok = ax.errorbar(irsfrok, radiosfr_ok, yerr=sfr_ok_uncertainty, xerr=irok_uncertainty, fmt='o', ecolor='k',
                      markeredgecolor='b', markerfacecolor='b', capsize=3, ms=marker_size, elinewidth=2.)

    non_detect = ax.errorbar(ir_non_detect, radio_non, yerr=2 * radio_non / 20, xerr=ir_non_unc, fmt='o', ecolor='k',
                              c='gold', capsize=3, uplims=True, marker='o', ms=marker_size, elinewidth=2.)

    # Plot the linear fits, the one-to-one line, and the detection limit dashed lines

    # fit_line = ax2.plot(np.linspace(0, x_axis_lim), fits[0](np.linspace(0, x_axis_lim)), 'g')
    one_to_one_line = ax.plot(np.linspace(0, x_axis_lim), one_to_one(np.linspace(0, x_axis_lim)), 'k', linestyle='--')
    fixed_line = ax.plot(np.linspace(0, x_axis_lim), np.poly1d([1. / 2.5, 0])(np.linspace(0, x_axis_lim)), color='k', alpha=0.7, linestyle=':')


    # copy of bottom axis to show luminosities
    tempax = ax.twinx()
    ax2 = tempax.twiny()


    # Titles
    # plt.suptitle('Star Formation Rate Comparison', y=0.91)
    fig.text(0.04, 0.5, '$\mathrm{SFR_{1.5 GHz} \ (M_{\odot} yr^{-1}})$', va='center', rotation='vertical', fontsize=24)
    fig.text(0.97, 0.5, '$\mathrm{L_{1.5 GHz}}$ (W Hz $^{-1}$)', va='center', rotation='vertical', fontsize=24)
    ax.set_xlabel('$\mathrm{SFR_{IR}} \ (\mathrm{M}_{\odot} \mathrm{yr}^{-1})$', fontsize=24)
    ax2.set_xlabel('$\mathrm{L_{IR}}$ (W)', fontsize=24)

    # put equations of linear fits on plot
    # plt.annotate('%s' % fits[0], (45, 180))
    # plt.annotate('%s' % fits[1], (50, 80))

    # ratio between radio luminosity (erg/s/Hz) and SFR. Multiplying SFR bounds by this number translates them to lumiosity bounds (W/Hz)
    radio_conversion = (t_ok['Luminosity'][0]/1e7)/t_ok['21 cm SFR'][0]

    low_sfr_lim = 10.

    # Log scales, axis limits
    #ax.set_ylim(flagRadioSFR[0] - 10000, flagRadioSFR[0] + 10000)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_yscale('log')


    ax.set_ylim(low_sfr_lim, y_axis_lim)
    ax2.set_ylim(low_sfr_lim * radio_conversion, y_axis_lim * radio_conversion)

    # multiplying an IR SFR limit by this number converts to a luminosity in Watts
    # 3.88e-37 is murphy2011 factor
    ir_conversion = 1./(3.88e-37)

    ax.set_xlim(low_sfr_lim, x_axis_lim)
    ax2.set_xlim(low_sfr_lim * ir_conversion, x_axis_lim * ir_conversion)

    """# star forming population defined by any galaxy with SFG = True
    SFG = t_survey[np.where(np.array(t_survey['SFG']) == 'T')[0]]
    # AGN population defined by any galaxy with EITHER Xray, MIR, or SED flagged AGN
    AGN = t_survey[np.where((np.array(t_survey['XrayAGN']) == 'T') | (np.array(t_survey['MIRAGN']) == 'T') | (
    np.array(t_survey['SEDAGN']) == 'T'))[0]]

    # plot contours of density of points on plot
    SFG_plot = sns.kdeplot(((10 ** SFG['logLTIRSF']) * (l_solar)), (10 ** SFG['logL21cm']), ax=ax3, n_levels=4,
                           cmap='Blues', alpha=0.3)
    AGN_plot = sns.kdeplot(((10 ** AGN['logLTIRSF']) * (l_solar)), (10 ** AGN['logL21cm']), ax=ax3, n_levels=4,
                           cmap='Reds', alpha=0.3)
    agn_legend = Line2D([0], [0], marker='o', color='w', label='Scatter',
                          markerfacecolor='w', markeredgecolor='darkred')
    SFG_legend = Line2D([0], [0], marker='o', color='w', label='Scatter',
                          markerfacecolor='w', markeredgecolor='midnightblue')"""

    # Legend

    legs = plt.legend((ok, non_detect),
                      ('Detections', '$3\sigma$ Upper Limits'), prop={'size': 18},
                      loc=2)





    ax.tick_params(axis='both', which='both', labelsize=16)
    ax2.tick_params(axis='both', which='both', labelsize=16)
    tempax.yaxis.set_tick_params(labelsize=16)

    ax.text(14.5, 16, '$\mathrm{SFR}_{\mathrm{IR}} = \mathrm{SFR}_{\mathrm{1.5 GHz}}$', rotation=47, rotation_mode='anchor', fontsize=16)

    ax.text(39.6, 13, r'$\mathrm{SFR}_{\mathrm{IR}} = 2.5 \times \mathrm{SFR}_{\mathrm{1.5 GHz}}$', rotation=46,
             rotation_mode='anchor', fontsize=16)



    plt.savefig('SFRs.pdf', overwrite=True, bbox_inches='tight', dpi=400)
    plt.clf()
    plt.cla()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()

# deprecated version of lum vs speed, comparing radio lum and outflow v linearly
"""def lum_v_speed():


	lum_ok = np.array(t_ok['Luminosity'])
	lum_uncertainty_ok = np.array(t_ok['Luminosity Error (stat.)'])
	z_ok = np.array(t_ok['v_out'])
	z_uncertainty_ok = 0
	names_ok = t_ok['Name']

	lum_bad = np.array(t_bad['Luminosity'])
	lum_uncertainty_bad = np.array(t_bad['Luminosity Error (stat.)'])
	z_bad = np.array(t_bad['v_out'])
	z_uncertainty_bad = 0
	names_bad = t_bad['Name']

	lum_non = np.array(t_nondetect['Luminosity'])
	lum_uncertainty_non = np.array(t_nondetect['Luminosity Error (stat.)'])
	z_non = np.array(t_nondetect['v_out'])
	z_uncertainty_non = 0
	names_non = t_nondetect['Name']

	fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 10), dpi=300, gridspec_kw={'height_ratios': [1, 4]})

	ok = ax2.errorbar(z_ok, lum_ok, yerr=lum_uncertainty_ok, xerr=z_uncertainty_ok,
					  fmt='o', ecolor='k', capsize=2, c='b', ms=marker_size)
	bad = ax.errorbar(z_bad, lum_bad, yerr=lum_uncertainty_bad, xerr=z_uncertainty_bad,
					   fmt='o', ecolor='k', capsize=2, c='r', marker='x', ms=marker_size)
	non = ax2.errorbar(z_non, lum_non, yerr=lum_non/5, xerr=z_uncertainty_non,
					   fmt='o', ecolor='k', capsize=2, c='gold', uplims=True, marker='v', ms=marker_size)

	ax.legend((ok, bad, non), ('Detections', 'AGN', 'Non-Detection Upper Limits'))

	#plt.suptitle('Luminosity vs. Wind Speed', y=0.92)
	fig.text(0.05, 0.5, '1.5 GHz Luminosity (erg s$^{-1}$ Hz$^{-1}$)', va='center', rotation='vertical', fontsize=18)

	plt.yscale('log')
	plt.xlabel('Outflow Velocity (km s$^{-1}$)', fontsize=18)

	# Hack to make the diagonal hashes on broken axis
	d = .015  # how big to make the diagonal lines in axes coordinates
	# arguments to pass to plot, just so we don't keep repeating them
	kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
	ax.plot((-d, +d), (-3.5 * d, 3.5 * d), **kwargs)  # top-left diagonal
	ax.plot((1 - d, 1 + d), (-3.5 * d, +3.5 * d), **kwargs)  # top-right diagonal
	# plt.gca().set_aspect('equal', adjustable='box')

	kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
	ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
	ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

	# hide the spines between ax and ax2
	ax.spines['bottom'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax.xaxis.tick_top()
	ax.tick_params(labeltop='off')  # don't put tick labels at the top
	ax2.xaxis.tick_bottom()


	if label_points:
		for x in range(len(lum_ok)):
			plt.annotate((names_ok[x].split('.')[0])[:5], (z_ok[x], lum_ok[x]), xytext=(0, 2), textcoords='offset points',
						 ha='right', va='bottom')
		for x in range(len(lum_bad)):
			plt.annotate((names_bad[x].split('.')[0])[:5], (z_bad[x], lum_bad[x]), xytext=(0, 2), textcoords='offset points',
						 ha='right', va='bottom')

		for x in range(len(lum_non)):
			plt.annotate((names_non[x].split('.')[0])[:5], (z_non[x], lum_non[x]), xytext=(0, 2), textcoords='offset points',
						 ha='right', va='bottom')
	plt.savefig('lum_vs_speed.png', overwrite=True, bbox_inches=0)
	plt.clf()
	plt.cla()
	plt.close()
	plt.clf()
	plt.cla()
	plt.close()
	plt.clf()
	plt.cla()
	plt.close()"""


# compare sfr derived from radio vs outflow speed (exclude J0827, the radio loud AGN)
def sfr_v_speed():
    t_fine = t[np.where((t['21 cm SFR'] < 1000))[0]]

    radsfrs = np.log10(np.array(t_fine['21 cm SFR'])/(np.pi * np.square(np.array(t_fine['Re']))))
    irsfrs = np.log10(np.array(t_fine['MySFR'])/(np.pi * np.square(np.array(t_fine['Re']))))
    vsa = np.array(t_fine['v_out'])
    idxs1 = np.where(vsa > 0)[0]
    radsfrs = radsfrs[idxs1]
    vsa = vsa[idxs1]
    irsfrs = irsfrs[idxs1]


    # detections

    SFR_ok = np.log10(np.array(t_ok['21 cm SFR'])/(np.pi * np.square(np.array(t_ok['Re']))))

    SFR_uncertainty_ok = (np.array(t_ok['21 cm SFR Error (stat.)'])/np.array(t_ok['21 cm SFR']))/np.log(10)
    irsfr_ok = np.log10(np.array(t_ok['MySFR'])/(np.pi * np.square(np.array(t_ok['Re']))))
    irsfr_unc_ok = (np.array(t_ok['MySFR Err'])/np.array(t_ok['MySFR']))/np.log(10)
    v_ok = np.array(t_ok['v_out'])
    idxs = np.where(v_ok > 0)[0]
    SFR_ok = SFR_ok[idxs]
    v_ok = v_ok[idxs]
    SFR_uncertainty_ok = SFR_uncertainty_ok[idxs]
    irsfr_ok = irsfr_ok[idxs]
    irsfr_unc_ok = irsfr_unc_ok[idxs]
    v_uncertainty_ok = 0
    names_ok = t_ok['Name']

    # non-detections

    SFR_non = np.log10(np.array(t_nondetect['21 cm SFR'])/(np.pi * np.square(np.array(t_nondetect['Re']))))

    SFR_uncertainty_non = np.array(t_nondetect['21 cm SFR Error (stat.)'])/np.array(t_nondetect['21 cm SFR'])/np.log(10)
    irsfr_non = np.log10(np.array(t_nondetect['MySFR'])/(np.pi * np.square(np.array(t_nondetect['Re']))))
    irsfr_unc_non = np.array(t_nondetect['MySFR Err'])/np.array(t_nondetect['MySFR'])/np.log(10)
    v_non = np.array(t_nondetect['v_out'])
    idxs2 = np.where(v_non > 0)[0]
    SFR_non = SFR_non[idxs2]
    v_non = v_non[idxs2]
    irsfr_non = irsfr_non[idxs2]
    irsfr_unc_non = irsfr_unc_non[idxs2]
    v_uncertainty_non = 0
    names_non = t_nondetect['Name']



    """all_vs = np.log10(np.concatenate((v_ok, v_non)))
    all_IR = np.concatenate([irsfr_ok, irsfr_non])
    all_rad = np.concatenate([SFR_ok, SFR_non])
    print(list(radsfrs), list(np.log10(vsa)))"""



    ir_coeff = st.spearmanr(list(irsfrs), list(np.log10(vsa)))
    print(ir_coeff[1])



    avgerr = np.average(irsfr_unc_ok)
    avgerr2 = np.average(irsfr_unc_non)
    avgerrfinal = (avgerr+avgerr2)/2

    fig = plt.figure(45, (10, 5))

    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2, sharey=ax1)

    plt.setp(ax2.get_yticklabels(), visible=False)

    plt.style.use('default')
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"


    ok = ax1.scatter(SFR_ok, v_ok, c='b', s=5*marker_size, zorder=3, label='Radio Detections')
    non = ax1.scatter(SFR_non, v_non, c='gold', marker='o', s=5*marker_size, zorder=3, label=r'Radio 3$\sigma$ Upper Limits')
    ax1.errorbar(SFR_non, v_non, xerr=0.1, fmt='o', ecolor='k', xuplims=True, marker='None', zorder=2)


    ax1.axvline(3.47, c='k', linestyle='--')
    ax1.text(3.3, 500, 'Eddington Limit', rotation=90)
    ax1.text(0.9, 2200, r'$\rho = %s$' % 0.46)

    ax2.axvline(3.47, c='k', linestyle='--')
    ax2.text(3.3, 500, 'Eddington Limit', rotation=90)
    ax2.text(0.9, 2200, r'$\rho = %s$' % round(ir_coeff[0], 2))

    ir = ax2.scatter(irsfr_ok, v_ok, c='r', s=5*marker_size, zorder=1, label='IR Detections')
    ax2.scatter(irsfr_non, v_non, c='r', s=5*marker_size, zorder=1)

    ax1.errorbar(1, 1000, xerr=avgerrfinal, ecolor='gray', capsize=2)
    ax2.errorbar(1, 1000, xerr=avgerrfinal, ecolor='gray', capsize=2)


    plt.yscale('log')
    ax1.set_ylabel(r'Average Mg$_{\mathrm{II}}$ Outflow Velocity (km s$^{-1})$', fontsize=12)

    ax1.set_xlabel(r'log$\left(\mathrm{\frac{\Sigma_{SFR_{1.5 GHz}}}{M_{\odot} yr^{-1} kpc^{-2}}} \right)$', fontsize=12)
    ax2.set_xlabel(r'log$\left(\mathrm{\frac{\Sigma_{SFR_{IR}}}{M_{\odot} yr^{-1} kpc^{-2}}} \right)$', fontsize=12)

    #plt.legend(loc='upper left', prop={'size': 8})

    if label_points:
        for x in range(len(SFR_ok)):
            plt.annotate((names_ok[x].split('.')[0])[:5], (v_ok[x], SFR_ok[x]), xytext=(0, 2),
                         textcoords='offset points',
                         ha='right', va='bottom')

        for x in range(len(SFR_non)):
            plt.annotate((names_non[x].split('.')[0])[:5], (v_non[x], SFR_non[x]), xytext=(0, 2),
                         textcoords='offset points',
                         ha='right', va='bottom')

    ax2.tick_params(axis='y', which='both', left=False, labelleft=False)
    plt.subplots_adjust(wspace=0.05)
    plt.savefig('sigma_vs_speed.pdf', overwrite=True, bbox_inches='tight', dpi=400)
    plt.clf()
    plt.cla()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close('all')


def plot_lum_vs_z():
    es_to_W = 10 ** 7

    lum_ok = np.array(t_ok['Luminosity']) / es_to_W
    lum_uncertainty_ok = np.array(t_ok['Luminosity Error (stat.)']) / es_to_W
    z_ok = np.array(t_ok['Z'])
    z_uncertainty_ok = 0
    names_ok = t_ok['Name']

    lum_bad = np.array(t_bad['Luminosity']) / es_to_W
    lum_uncertainty_bad = np.array(t_bad['Luminosity Error (stat.)']) / es_to_W
    z_bad = np.array(t_bad['Z'])
    z_uncertainty_bad = 0
    names_bad = t_bad['Name']

    lum_non = np.array(t_nondetect['Luminosity']) / es_to_W
    lum_uncertainty_non = np.array(t_nondetect['Luminosity Error (stat.)']) / es_to_W
    z_non = np.array(t_nondetect['Z'])
    z_uncertainty_non = 0
    names_non = t_nondetect['Name']

    print(st.spearmanr(list(z_non), list(np.log10(lum_non))))

    plt.style.use('default')
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"

    fig = plt.figure(0, (8, 8))

    ok = plt.errorbar(z_ok, lum_ok, yerr=lum_uncertainty_ok, xerr=z_uncertainty_ok,
                      fmt='o', ecolor='k', capsize=2, c='b', ms=marker_size)
    bad = plt.errorbar(z_bad, lum_bad, yerr=lum_uncertainty_bad, xerr=z_uncertainty_bad,
                      fmt='o', ecolor='k', capsize=2, c='r', marker='x', ms=marker_size)
    non = plt.errorbar(z_non, lum_non, yerr=lum_non / 5, xerr=z_uncertainty_non,
                       fmt='o', ecolor='k', capsize=2, c='gold', uplims=True, marker='o', ms=marker_size)

    plt.legend((ok, bad, non), ('Detections', 'Radio AGN', '$3\sigma$ Upper Limits'), loc=2, fontsize=15)

    fig.text(0, 0.5, '$\mathrm{L_{1.5GHz}} \ (\mathrm{W \ Hz}^{-1}$)', va='center', rotation='vertical', fontsize=25)

    plt.yscale('log')
    plt.xlabel('$z$', fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=16)

    plt.savefig('lum_vs_z.pdf', overwrite=True, bbox_inches='tight', dpi=400)
    plt.clf()
    plt.cla()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()





def Lum_hist():
    plt.clf()
    plt.cla()
    plt.close()
    plt.style.use('default')
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"
    lums_ok = list(np.array(t_ok['Luminosity']) / (10000000.))
    lums_non = list(np.array(t_nondetect['Luminosity']) / (10000000.))
    plt.clf()
    plt.cla()
    plt.close()
    plt.close()
    plt.close()
    plt.close()
    #fig = plt.figure(24, (4,3))
    colors = ['b', 'gold']
    plt.hist([lums_ok, lums_non], bins=20, stacked=True, color=colors, label=['Detections', '$\mathrm{3 \sigma \ Upper \ Limits}$'])
    plt.legend()
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
    plt.xlabel('$\mathrm{L_{1.5 GHz} \ (W \ Hz}^{-1}$)', fontsize=18)
    plt.ylabel('Number', fontsize=18)
    plt.savefig('lum_hist.png', overwrite=True, bbox_inches='tight', dpi=300)
    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()



def wise_colors():

    """"# w1-w2 vs w3-w4
    y_shift = np.array(t_shift['w1-w2']).astype(float)
    x_shift = np.array(t_shift['w3-w4']).astype(float)
    names = t_shift['Name']

    y_else = np.array(t_else['w1-w2']).astype(float)
    x_else = np.array(t_else['w3-w4']).astype(float)
    else_names = t_else['Name']

    plt.clf()
    plt.cla()
    plt.close()
    plt.close()
    plt.close()
    plt.close()
    fig = plt.figure(81)
    shifted = plt.scatter(x_shift, y_shift, c='g')
    other = plt.scatter(x_else, y_else, c='m')
    plt.xlabel('W3-W4')
    plt.ylabel('W1-W2')
    plt.legend((shifted, other), ('Shifted', 'Other'))

    if label_points:
        for x in range(len(y_shift)):
            plt.annotate((names[x].split('.')[0])[:5], (x_shift[x], y_shift[x]), xytext=(0, 2),
                         textcoords='offset points', ha='right', va='bottom')
        for x in range(len(y_else)):
            plt.annotate((else_names[x].split('.')[0])[:5], (x_else[x], y_else[x]), xytext=(0, 2),
                         textcoords='offset points', ha='right', va='bottom')

    plt.savefig('wisecolors.png', overwrite=True, bbox_inches='tight', dpi=300)

    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()"""

    ######################################
    # w2-w3
    #

    Template_colors = True
    templates = sorted(
        glob.glob('/Users/graysonpetter/Desktop/Dartmouth/HIZEA/hizea-VLA-SFRs/Comprehensive_library/*.txt'))
    # wise bandpasses
    bandpass = sorted(glob.glob('/Users/graysonpetter/Desktop/Dartmouth/HIZEA/hizea-VLA-SFRs/bandpass/*.txt'))
    template_names = []
    # corrections to go from AB mags to WISE Vega mags. first entry is 0 so that W1 correction corresponds to 1st index
    wise_corr = [0, 2.699, 3.339, 5.174, 6.620]

    w_one_two_AB, w_two_three_AB, w_three_four_AB = [], [], []
    w_one_two, w_two_three, w_three_four, w_four = [], [], [], []

    # Simulate WISE galaxy colors for different types of templates
    if Template_colors:
        simulation = Templates.simulate_wise_fluxes_for_colors(0.6, templates, bandpass, False)
        lums_for_all = simulation[0]
        template_names = simulation[1]


        for set in lums_for_all:
            w_one_two.append((-2.5 * np.log10(set[0] / set[1])) - wise_corr[1] + wise_corr[2])
            w_two_three.append((-2.5 * np.log10(set[1] / set[2])) - wise_corr[2] + wise_corr[3])
            w_three_four.append(-2.5 * np.log10(set[2]/ set[3]) - wise_corr[3] + wise_corr[4])
            w_one_two_AB.append(-2.5 * np.log10(set[0] / set[1]))
            w_two_three_AB.append(-2.5 * np.log10(set[1] / set[2]))
            w_three_four_AB.append(-2.5 * np.log10(set[2]/ set[3]))
            w_four.append(set[3])

    t_fine = t[np.where((t['21 cm SFR'] < 1000))[0]]

    one_two, two_three, one_two_err, two_three_err = [], [], [], []
    names = t_fine['Name']
    ir_rad_ratio = np.array(np.array(t_fine[sfr_to_use]/t_fine['21 cm SFR']))
    ir_rad_detect = np.array(np.array(t_ok[sfr_to_use] / t_ok['21 cm SFR']))
    ir_rad_non = np.array(np.array(t_nondetect[sfr_to_use] / t_nondetect['21 cm SFR']))

    for name in names:
        wisecolrs = WISE.colors(name[:5])
        one_two.append(wisecolrs[0])
        one_two_err.append(wisecolrs[3])
        two_three.append(wisecolrs[2])
        two_three_err.append(wisecolrs[5])

    names_detect = t_ok['Name']
    names_non = t_nondetect['Name']

    onetwo_detect, twothree_detect = [], []

    for name in names_detect:
        wisecolrs = WISE.colors(name[:5])
        onetwo_detect.append(wisecolrs[0])
        twothree_detect.append(wisecolrs[2])

    onetwo_non, twothree_non = [], []

    for name in names_non:
        wisecolrs = WISE.colors(name[:5])
        onetwo_non.append(wisecolrs[0])
        twothree_non.append(wisecolrs[2])

    badcolors = WISE.colors('J0827')
    onetwo_bad = badcolors[0]
    two_three_bad = badcolors[2]
    print(one_two_err, two_three_err)

    avg_one_two_err = np.median(one_two_err)
    avg_two_three_err = np.median(two_three_err)

    fig = plt.figure(82)
    plt.style.use('default')
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"

    cmap = matcol.ListedColormap(['#CE03CE', '#D588D5', '#75D77D', '#5AC162', '#45AD4D', '#309A38', '#1E8526', '#13741B', '#0B6011', '#044D0A'])


    boundaries = [0, 0.5, 1., 1.5, 2., 2.5, 3, 3.5, 4., 4.5]
    norm = matcol.BoundaryNorm(boundaries, cmap.N, clip=True)
    colors_detect = plt.scatter(twothree_detect, onetwo_detect, c=ir_rad_detect, cmap=cmap, norm=norm, s=(marker_size*6))
    colors_non = plt.scatter(twothree_non, onetwo_non, c=ir_rad_non, cmap=cmap, norm=norm, s=(marker_size * 6))
    colors_bad = plt.scatter(two_three_bad, onetwo_bad, c='r', marker='x', s=marker_size*4)
    colors_non.set_facecolor('none')
    std_errbar = plt.errorbar(2.0, 0.7, yerr=avg_one_two_err, xerr=avg_two_three_err, ecolor='gray', elinewidth=0.5, capsize=2)
    cbar = plt.colorbar(colors_detect)
    cbar.set_label('$\mathrm{SFR_{IR}/SFR_{1.5\mathrm{GHz}}}$', size=12)

    agn_model = [w_two_three[:4], w_one_two[:4], template_names[:4]]
    comp_model = [w_two_three[4:8], w_one_two[4:8], template_names[4:8]]
    sfg_model = [w_two_three[8:11], w_one_two[8:11], template_names[8:11]]

    m1, = plt.plot([], [], c='magenta', marker='o', markersize=10,
                  fillstyle='left', linestyle='none', markeredgecolor='none')

    m2, = plt.plot([], [], c='g', marker='o', markersize=10,
                  fillstyle='right', linestyle='none', markeredgecolor='none')


    agn = plt.scatter(agn_model[0], agn_model[1], c='r', marker='+', alpha=0.5, s=30)
    comp = plt.scatter(comp_model[0], comp_model[1], c='#FF9E00', marker='D', alpha=0.5, s=15)
    sfg = plt.scatter(sfg_model[0], sfg_model[1], c='b', marker='*', s=40, alpha=0.5)
    plt.xlabel('W2-W3 (mag)', size=12)
    plt.ylabel('W1-W2 (mag)', size=12)
    plt.legend(((m2, m1), colors_bad, agn, comp, sfg), ('This sample', 'Radio-loud AGN', 'AGN Templates', 'Composite Templates', 'SFG Templates'))
    #ax = plt.gca()
    #legend = ax.get_legend()
    #legend.legendHandles[0].set_color(plt.cm.Reds(.8))

    """if label_points:
        for x in range(len(one_two)):
            plt.annotate((names[x].split('.')[0])[:5], (two_three[x], one_two[x]), xytext=(0, 2),
                         textcoords='offset points',
                         ha='right', va='bottom')

        #for x in range(len(w_two_three)):
            #plt.annotate(template_names[x], (w_two_three[x], w_one_two[x]), xytext=(0, 2), textcoords='offset points',
             #            ha='right', va='bottom')"""

    plt.savefig('wisecolors.pdf', overwrite=True, bbox_inches='tight', dpi=400)

    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close('all')

    """sofia_colors = True
    ######################################
    # sofia c-e vs w4-c


    sofia_bands = sorted(glob.glob('/Users/graysonpetter/Desktop/Dartmouth/HIZEA/hizea-VLA-SFRs/SOFIA/HAWC*.csv'))
    # templates.extend(sorted(glob.glob('/Users/graysonpetter/Desktop/Dartmouth/HIZEA/hizea-VLA-SFRs/swire_lib/*.sed')))
    examps = sorted(glob.glob('/Users/graysonpetter/Desktop/Dartmouth/HIZEA/hizea-VLA-SFRs/swire_lib/*.sed'))

    if sofia_colors:
        simulation_sof = simulate_colors(templates, sofia_bands, True)
        lums_for_all = simulation_sof[0]
        template_names = simulation_sof[1]

        sof_c_e, w_four_c = [], []
        for b in range(len(lums_for_all)):
            sof_c_e.append(-2.5 * np.log10(lums_for_all[b][0] / lums_for_all[b][1]))
            w_four_c.append(-2.5 * np.log10(w_four[b] / lums_for_all[b][0]))

        "simulation_ex_wise = simulate_colors(examps, bandpass, False)
        wise_exs = simulation_ex_wise[0]
        template_names.extend(simulation_ex_wise[1])


        simulation_ex_sof = simulate_colors(examps, sofia_bands, True)
        sof_exs = simulation_ex_sof[0]
        print wise_exs, sof_exs


        for b in range(len(sof_exs)):
            sof_c_e.append(-2.5 * np.log10(sof_exs[b][0] / sof_exs[b][1]))
            w_four_c.append(-2.5 * np.log10(wise_exs[b][3] / sof_exs[b][0]))"

        # HAWC C-E vs W4-C

        fig = plt.figure(37)

        sof_wise = plt.scatter(sof_c_e, w_four_c, c='k')
        plt.xlabel('HAWC C-E')
        plt.ylabel('W4-C')

        sofia_observed = pd.read_csv('/Users/graysonpetter/Desktop/Dartmouth/HIZEA/hizea-VLA-SFRs/SOFIA/SOFIA_observed.csv', header=None, engine='python')
        obs_y = np.array(sofia_observed.iloc[:, 0])
        obs_mag = np.array(sofia_observed.iloc[:, 1])
        hawc_c_e = obs_mag[len(obs_mag)-2]-obs_mag[len(obs_mag)-1]
        y_w_c = obs_mag[len(obs_mag)-3] - obs_mag[len(obs_mag)-2]

        obs_point = plt.scatter(hawc_c_e, y_w_c, c='g')


        if label_points:
            for x in range(len(template_names)):
                plt.annotate(template_names[x], (sof_c_e[x], w_four_c[x]), xytext=(0, 2),
                             textcoords='offset points',
                             ha='right', va='bottom')
            plt.annotate('J1613', (hawc_c_e, y_w_c), xytext=(0, 2),
                             textcoords='offset points',
                             ha='right', va='bottom')

        plt.savefig('sofia_colors.png', overwrite=True, bbox_inches='tight', dpi=300)

        plt.clf()
        plt.close()
        plt.clf()
        plt.cla()
        plt.close()

        "del sof_c_e[-1]
        del sof_c_e[-1]
        del template_names[-1]
        del template_names[-1]

        # HAWC C-E vs W2-W3

        fig = plt.figure(178)

        plt.scatter(sof_c_e, w_two_three_AB, c='k')
        plt.xlabel('HAWC C-E')
        plt.ylabel('W2-W3')


        if label_points:
            for x in range(len(template_names)):
                plt.annotate(template_names[x], (sof_c_e[x], w_two_three_AB[x]), xytext=(0, 2),
                             textcoords='offset points',
                             ha='right', va='bottom')

        plt.savefig('sofia_colors2.png', overwrite=True, bbox_inches='tight', dpi=300)

        plt.clf()
        plt.close()
        plt.clf()
        plt.cla()
        plt.close()

        # HAWC C-E vs W3-W4

        fig = plt.figure(6969)

        plt.scatter(sof_c_e, w_three_four_AB, c='k')
        plt.xlabel('HAWC C-E')
        plt.ylabel('W3-W4')

        w_three_four_obs = obs_mag[len(obs_mag)-4] - obs_mag[len(obs_mag)-3]
        plt.scatter(hawc_c_e, w_three_four_obs, c='g')

        if label_points:
            for x in range(len(template_names)):
                plt.annotate(template_names[x], (sof_c_e[x], w_three_four_AB[x]), xytext=(0, 2),
                             textcoords='offset points',
                             ha='right', va='bottom')

        plt.savefig('sofia_colors3.png', overwrite=True, bbox_inches='tight', dpi=300)

        plt.clf()
        plt.close()
        plt.clf()
        plt.cla()
        plt.close()"


    #################################
    # ir/radio vs w2-w3

    devi = np.array(t_shifto['IR SFR'] / t_shifto['21 cm SFR']).astype(float)
    col = np.array(t_shifto['w2-w3']).astype(float)
    namex = t_shifto['Name']

    devio = np.array(abs(t_elseo['IR SFR'] / t_elseo['21 cm SFR'])).astype(float)
    colo = np.array(t_elseo['w2-w3']).astype(float)
    else_namex = t_elseo['Name']

    fig = plt.figure(87)
    shiftedini = plt.scatter(col, devi, c='g')
    othelo = plt.scatter(colo, devio, c='m')
    plt.xlabel('W2-W3')
    plt.ylabel('IR SFR/Radio SFR')
    plt.legend((shiftedini, othelo), ('Shifted', 'Other'))

    if label_points:
        for x in range(len(devi)):
            plt.annotate((namex[x].split('.')[0])[:5], (col[x], devi[x]), xytext=(0, 2), textcoords='offset points',
                         ha='right', va='bottom')
        for x in range(len(devio)):
            plt.annotate((else_namex[x].split('.')[0])[:5], (colo[x], devio[x]), xytext=(0, 2),
                         textcoords='offset points',
                         ha='right', va='bottom')

    plt.savefig('wisecolors3.png', overwrite=True, bbox_inches='tight', dpi=300)

    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()

    ##############################
    # ir/radio vs z



    deviat = np.array(abs(t_shifto['IR SFR'] / t_shifto['21 cm SFR'])).astype(float)
    coli = np.array(t_shifto['Z']).astype(float)
    nameu = t_shifto['Name']

    deviato = np.array(abs(t_elseo['IR SFR'] / t_elseo['21 cm SFR'])).astype(float)
    colio = np.array(t_elseo['Z']).astype(float)
    else_nameu = t_elseo['Name']

    fig = plt.figure(88)
    linguini = plt.scatter(coli, deviat, c='g')
    oreo = plt.scatter(colio, deviato, c='m')
    plt.xlabel('Z')
    plt.ylabel('IR SFR/Radio SFR')
    plt.legend((linguini, oreo), ('Shifted', 'Other'))

    if label_points:
        for x in range(len(deviat)):
            plt.annotate((nameu[x].split('.')[0])[:5], (coli[x], deviat[x]), xytext=(0, 2), textcoords='offset points',
                         ha='right', va='bottom')
        for x in range(len(deviato)):
            plt.annotate((else_nameu[x].split('.')[0])[:5], (colio[x], deviato[x]), xytext=(0, 2),
                         textcoords='offset points',
                         ha='right', va='bottom')

    plt.savefig('wisecolors4.png', overwrite=True, bbox_inches='tight', dpi=300)

    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()"""


def sed():

    # separate galaxies which lie significantly to the right of 1-to-1 line from remainder of sample
    t_shift = t[np.where(t[sfr_to_use] > (1.5 * t['21 cm SFR']))[0]]
    t_else = t[np.where(t[sfr_to_use] < (1.5 * t['21 cm SFR']))[0]]

    # wavelengths of SDSS UGRIZ + WISE 1,2,3,4
    waves = [.3543, .4770, .6231, .7625, .9134, 3.3680, 4.6180, 12.0820, 22.1940]
    plt.clf()
    plt.cla()
    plt.close()
    plt.close()
    plt.close()
    plt.close()
    fig = plt.figure(274)

    # WISE corrections from VEGA mags to AB
    vega_to_ab = [2.699, 3.339, 5.174, 6.620]

    standard = 0
    magtrack = []
    colors = ['r', 'g', 'b', 'k', 'y', 'c', 'm', 'teal', 'orange', 'gold', 'pink', 'lime', 'dimgray', 'seagreen', 'olive', 'saddlebrown', 'lightgray', 'purple', 'tan', 'darkblue']

    # Normalize all
    for x in range(len(t['Name'])): #for x in range(len(t_shift['Name'])):

        wise = Table.read('unWISE/%s.fits' % (t['Name'][x][:5])) #wise = Table.read('unWISE/%s.fits' % (t_shift['Name'][x][:5]))
        SDSS = Table.read('SDSS/%s.fits' % (t['Name'][x][:5])) #SDSS = Table.read('SDSS/%s.fits' % (t_shift['Name'][x][:5]))

        # normalizing
        if x == 0:
            standard = wise['w1_mag'] #standard = SDSS['I']
        diff = float(standard) - float(wise['w1_mag']) #diff = float(standard) - float(SDSS['I'])

        if np.isnan(wise['w4_mag']):

            mags = ([SDSS['u'], SDSS['g'], SDSS['r'], SDSS['I'], SDSS['z'], (wise['w1_mag'] + vega_to_ab[0]),
                    (wise['w2_mag'] + vega_to_ab[1]), (wise['w3_mag'] + vega_to_ab[2]), np.nan])
        else:
            mags = ([SDSS['u'], SDSS['g'], SDSS['r'], SDSS['I'], SDSS['z'], (wise['w1_mag'] + vega_to_ab[0]),
                    (wise['w2_mag'] + vega_to_ab[1]), (wise['w3_mag'] + vega_to_ab[2]),
                    (wise['w4_mag'] + vega_to_ab[3])])

        mags = [y + diff for y in mags]
        magtrack.append(list(np.array(mags).astype(float).flatten()))
        plt.plot(waves[5:9], list(np.array(mags).astype(float).flatten())[5:9], c=colors[x], linewidth=1, label=t['Name'][x][:5])
    plt.legend(loc=2)

    """median_mag = []
    shift_std = []
    median_w_one = 0
    for p in range(0, 9):
        scoop = [item[p] for item in magtrack]
        scoop[:] = [b for b in scoop if b > 0]
        median_mag.append(np.median(np.array(scoop)))
        if p == 5:
            median_w_one = np.median(np.array(scoop))
        shift_std.append(np.std(np.array(scoop)))
    plt.plot(waves, median_mag, c='g', linewidth=.5)
    plt.errorbar(waves, median_mag, yerr=shift_std, fmt='none', elinewidth=1, color='g', capsize=2)

    magtrack, median_mag, shift_std = [], [], []

    for x in range(len(t_else['Name'])):

        wise = Table.read('unWISE/%s.fits' % (t_else['Name'][x][:5]))
        SDSS = Table.read('SDSS/%s.fits' % (t_else['Name'][x][:5]))
        diff = float(standard) - float(SDSS['I'])

        #scale = float(standard) / float(SDSS['I'])

        if np.isnan(wise['w4_mag']):
            mags = ([SDSS['u'], SDSS['g'], SDSS['r'], SDSS['I'], SDSS['z'], (wise['w1_mag'] + vega_to_ab[0]),
                    (wise['w2_mag'] + vega_to_ab[1]), (wise['w3_mag'] + vega_to_ab[2]), 0])
            mags = [y + diff for y in mags]
            magtrack.append(list(np.array(mags).astype(float).flatten()))
            # plt.plot(waves, mags, c='m', linewidth=.5)
        else:
            mags = ([SDSS['u'], SDSS['g'], SDSS['r'], SDSS['I'], SDSS['z'], (wise['w1_mag'] + vega_to_ab[0]),
                    (wise['w2_mag'] + vega_to_ab[1]), (wise['w3_mag'] + vega_to_ab[2]),
                    (wise['w4_mag'] + vega_to_ab[3])])
            mags = [y + diff for y in mags]
            magtrack.append(list(np.array(mags).astype(float).flatten()))
    # plt.plot(waves, mags, c='m', linewidth=.5)

    for p in range(0, 9):
        scoop = [item[p] for item in magtrack]
        scoop[:] = [b for b in scoop if b > 0]
        median_mag.append(np.median(np.array(scoop)))
        shift_std.append(np.std(np.array(scoop)))
    plt.plot(waves, median_mag, c='m', linewidth=.5)
    plt.errorbar(waves, median_mag, yerr=shift_std, fmt='none', elinewidth=0.5, color='m', capsize=2)"""

    """template_sed = True
    if template_sed:
        templates = sorted(
            glob.glob('/Users/graysonpetter/Desktop/Dartmouth/HIZEA/hizea-VLA-SFRs/Comprehensive_library/*.txt'))
        for tem in templates:
            red_spec = Templates.redshift_spectrum(0.6, tem, False)
            red_waves = red_spec[4]
            lumi = red_spec[1]
            cut = np.where(np.array(red_waves) < waves[8])[0]
            red_waves = red_waves[cut]
            lumi = lumi[cut]
            mags = []
            normalize_wave = (np.abs(np.array(red_waves) - waves[5])).argmin()
            difference = median_w_one - (-2.5*np.log10(lumi[normalize_wave]))
            for x in range(len(red_waves)):
                mags.append(-2.5*np.log10(lumi[x])+difference)

            if str(tem.split('/')[len(tem.split('/'))-1]).startswith('A'):
                temp_color='r'
            elif str(tem.split('/')[len(tem.split('/'))-1]).startswith('C'):
                temp_color = 'purple'
            else:
                temp_color = 'b'

            plt.plot(red_waves, mags, c=temp_color, linewidth=0.1)"""

    plt.gca().invert_yaxis()
    plt.xlabel('$\lambda \ (\mathrm{\mu m})$')
    plt.xscale('log')
    plt.ylabel('m')
    plt.savefig('SED.png', overwrite=True, dpi=400)
    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()


def radio_to_ir():
    tot_table = pd.read_csv('table.csv')
    idxs = np.where(np.array(tot_table['21 cm SFR']) > 1000.)[0]
    tot_table.drop(idxs, inplace=True)
    tot_table.reset_index(drop=True, inplace=True)
    idxs = np.where(np.array(tot_table['21 cm Flux']) < 0)[0]
    tot_table.drop(idxs, inplace=True)
    tot_table.reset_index(drop=True, inplace=True)

    gal_names = np.array(tot_table['Name'])
    radio_fluxes = np.array(tot_table['21 cm Flux'])

    vega_to_ab = [2.699, 3.339, 5.174, 6.620]

    rad_twelve, rad_z, ir_rad_ratio = [], [], []
    for x in range(len(gal_names)):
        wise = Table.read('unWISE/%s.fits' % (gal_names[x][:5]))
        SDSS = Table.read('SDSS/%s.fits' % (gal_names[x][:5]))
        w_ab = float(wise['w3_mag'] + vega_to_ab[2])
        z_band = float(SDSS['z'])
        w_flux = 3631*(10**(-0.4*w_ab))
        z_flux = 3631*(10**(-0.4*z_band))
        ir_rad_ratio.append(float(tot_table[sfr_to_use][x]/tot_table['21 cm SFR'][x]))
        rad_twelve.append(w_flux/z_flux)
        rad_z.append(radio_fluxes[x]/z_flux)

    rad_twelve = np.log10(np.array(rad_twelve))
    rad_z = np.log10(np.array(rad_z))

    plt.figure(151)
    plt.style.use('default')
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"
    cmap = matcol.ListedColormap(['#CE03CE', '#D242D2', '#D588D5', '#75D77D', '#5AC162', '#45AD4D', '#309A38', '#1E8526', '#13741B', '#0B6011', '#044D0A'])
    boundaries = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
    norm = matcol.BoundaryNorm(boundaries, cmap.N, clip=True)
    scatter = plt.scatter(rad_z, rad_twelve, c=ir_rad_ratio, cmap=cmap, norm=norm)
    cbar = plt.colorbar(scatter)
    cbar.set_label('$\mathrm{SFR_{IR}/SFR_{Radio}}$')
    plt.xlabel('$\mathrm{log(f_{1.5 GHz}/f_{Z})}$')
    plt.ylabel('$\mathrm{log(f_{W3}/f_{Z})}$')
    plt.savefig('radio_IR.png', overwrite=True, dpi=400)
    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()


def ratio_to_size():
    tot_table = pd.read_csv('table.csv')
    idxs = np.where(np.array(tot_table['21 cm SFR']) > 1000.)[0]
    tot_table.drop(idxs, inplace=True)
    tot_table.reset_index(drop=True, inplace=True)


    gal_names_detect = np.array(t_ok['Name'])
    zs_detect = np.array(t_ok['Z'])

    gal_names_non = np.array(t_nondetect['Name'])
    zs_non = np.array(t_nondetect['Z'])

    ir_rad_ratio_detect = np.array(t_ok[sfr_to_use]) / np.array(t_ok['21 cm SFR'])
    eff_rad_detect = np.array(t_ok['Re'])

    ir_rad_ratio_non = np.array(t_nondetect[sfr_to_use]) / np.array(t_nondetect['21 cm SFR'])
    eff_rad_non = np.array(t_nondetect['Re'])

    ratio_err = np.array(tot_table[sfr_to_use])/np.array(tot_table['21 cm SFR'])*np.sqrt(np.square(tot_table[err_to_use]/tot_table[sfr_to_use])+np.square(tot_table['21 cm SFR Error (stat.)']/tot_table['21 cm SFR']))
    avg_ratio_err = np.nanmean(ratio_err)

    fig = plt.figure(1419, (10, 5))
    fig.set_canvas(plt.gcf().canvas)
    plt.style.use('default')
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"

    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2, sharey=ax1)

    plt.setp(ax2.get_yticklabels(), visible=False)


    ir_sfr_detect = np.array(t_ok[sfr_to_use])
    ir_sfr_non = np.array(t_nondetect[sfr_to_use])


    half_detect = False
    compactness_detect, compact_err_detect = [], []
    location = '/Users/graysonpetter/Desktop/mac_copy/'
    for x in range(len(gal_names_detect)):
        data = fits.open(location+gal_names_detect[x]+'/'+gal_names_detect[x][:5]+'_HST.fits')[0].data
        cleaned = data[~np.isnan(data)]
        rms = biweight_midvariance(cleaned)

        background = Background2D(data, (500, 500))

        cosmo = FlatLambdaCDM(H0=70., Om0=0.3)
        pix_scale = 0.025   # arcsec/pix
        ang_scale = cosmo.kpc_proper_per_arcmin(zs_detect[x]).value/60.  # kpc/arcsec
        size_scale = ang_scale*pix_scale    # kpc/pixel
        ap_radii_phys = [0.5, 1.0, 5.]   # kpc
        ap_radii_pix = list(np.array(ap_radii_phys)/size_scale)     # pix
        apertures = [CircularAperture((1200., 1200.), r=r) for r in ap_radii_pix]
        tab = aperture_photometry(data - background.background, apertures, error=rms)

        if half_detect:
            compactness_detect.append(tab['aperture_sum_0']/tab['aperture_sum_2'])
            compact_err_detect.append(tab['aperture_sum_0']/tab['aperture_sum_2']*np.sqrt((tab['aperture_sum_err_0']/tab['aperture_sum_0'])**2+(tab['aperture_sum_err_2']/tab['aperture_sum_2'])**2))

        else:
            compactness_detect.append(tab['aperture_sum_1']/tab['aperture_sum_2'])
            #print(tab['aperture_sum_1'])
            #print(tab['aperture_sum_err_1'])
            #print('\n \n \n')
            compact_err_detect.append(tab['aperture_sum_1'] / tab['aperture_sum_2'] * np.sqrt(
                (tab['aperture_sum_err_1'] / tab['aperture_sum_1']) ** 2 + (
                tab['aperture_sum_err_2'] / tab['aperture_sum_2']) ** 2))

    half_non = False
    compactness_non, compact_err_non = [], []
    location = '/Users/graysonpetter/Desktop/mac_copy/'
    for x in range(len(gal_names_non)):
        data = fits.open(location + gal_names_non[x] + '/' + gal_names_non[x][:5] + '_HST.fits')[0].data
        cleaned = data[~np.isnan(data)]
        rms = biweight_midvariance(cleaned)

        background = Background2D(data, (500, 500))

        cosmo = FlatLambdaCDM(H0=70., Om0=0.3)
        pix_scale = 0.025  # arcsec/pix
        ang_scale = cosmo.kpc_proper_per_arcmin(zs_non[x]).value / 60.  # kpc/arcsec
        size_scale = ang_scale * pix_scale  # kpc/pixel
        ap_radii_phys = [0.5, 1.0, 5.]  # kpc
        ap_radii_pix = list(np.array(ap_radii_phys) / size_scale)  # pix
        apertures = [CircularAperture((1200., 1200.), r=r) for r in ap_radii_pix]
        tab = aperture_photometry(data - background.background, apertures, error=rms)

        if half_non:
            compactness_non.append(tab['aperture_sum_0'] / tab['aperture_sum_2'])
            compact_err_non.append(tab['aperture_sum_0'] / tab['aperture_sum_2'] * np.sqrt(
                (tab['aperture_sum_err_0'] / tab['aperture_sum_0']) ** 2 + (
                tab['aperture_sum_err_2'] / tab['aperture_sum_2']) ** 2))

        else:
            compactness_non.append(tab['aperture_sum_1'] / tab['aperture_sum_2'])
            # print(tab['aperture_sum_1'])
            # print(tab['aperture_sum_err_1'])
            # print('\n \n \n')
            compact_err_non.append(tab['aperture_sum_1'] / tab['aperture_sum_2'] * np.sqrt(
                (tab['aperture_sum_err_1'] / tab['aperture_sum_1']) ** 2 + (
                    tab['aperture_sum_err_2'] / tab['aperture_sum_2']) ** 2))

    average_sfrs = np.concatenate([ir_sfr_non, ir_sfr_detect], axis=0)
    min_, max_ = average_sfrs.min(), average_sfrs.max()

    ax2.annotate("", xy=(0.4, 4.1), xytext=(0.53, 4.1), arrowprops = dict(arrowstyle="<-"))
    ax2.text(0.51, 4.2, 'More Compact')
    #ax2.text('More Compact', 0.5, 3.2)
    ax2.set_xlabel(r'Compactness $\mathrm{\left(\frac{F_{1 kpc}}{F_{5 kpc}}\right)}$', size=12)
    avg_comp_err = np.average(np.concatenate((np.array(compact_err_detect), np.array(compact_err_non))))
    scatter2_detect = ax2.scatter(compactness_detect, ir_rad_ratio_detect, c=ir_sfr_detect, cmap='plasma', zorder=3)
    #plt.clim(min_, max_)
    scatter2_non = ax2.scatter(compactness_non, ir_rad_ratio_non, c=ir_sfr_non, cmap='plasma', marker='o', zorder=2)
    error2_non = ax2.errorbar(compactness_non, ir_rad_ratio_non, lolims=True, yerr=0.15, ecolor='k', ls='none', zorder=1)
    #plt.clim(min_, max_)
    ax2.errorbar(0.41, 2.9, yerr=avg_ratio_err, xerr=avg_comp_err, ecolor='gray', elinewidth=1, capsize=2)
    ax2.invert_xaxis()
    #plt.ylabel('$\mathrm{SFR_{IR}/SFR_{Radio}}$')
    #cbar = plt.colorbar(scatter)
    #cbar.set_label(r'$\mathrm{\overline{SFR}}$')

    groveslums = 10 ** np.array([42.119, 42.028, 41.936, 41.809, 41.75, 41.732])
    groveslumstwo = 10 ** np.array([42.138, 42.056, 41.968, 41.809, 41.893, 41.795])
    grovesrs = 200 * (10 ** (-1. / 2 * np.array([6.5, 6., 5.5, 5., 4.5, 4.])))

    xnew = np.linspace(grovesrs.min(), grovesrs.max(), 100)
    smoothfit = UnivariateSpline(grovesrs, groveslums, k=2, s=1)
    smoothfittwo = UnivariateSpline(grovesrs, groveslumstwo, k=2, s=1)
    cofs = smoothfit.get_coeffs()

    grovescofs = np.polyfit(grovesrs, groveslums, 2)
    grovescofstwo = np.polyfit(grovesrs, groveslumstwo, 2)

    def grovesfit(rs, a):
        return a * (grovescofs[0] * rs ** 2 + grovescofs[1] * rs + grovescofs[2])
    def grovesfittwo(rs, a):
        return a * (grovescofstwo[0] * rs ** 2 + grovescofstwo[1] * rs + grovescofstwo[2])

    norms, pcov = curve_fit(grovesfit, eff_rad_detect, ir_rad_ratio_detect)
    normstwo, pcovtwo = curve_fit(grovesfittwo, eff_rad_detect, ir_rad_ratio_detect)

    ax1.annotate("", xy=(2.1, 4.1), xytext=(1.3, 4.1), arrowprops=dict(arrowstyle="<-"))
    ax1.text(1.4, 4.2, 'More Compact')
    scatter1_detect = ax1.scatter(eff_rad_detect, ir_rad_ratio_detect, c=ir_sfr_detect, cmap='plasma', zorder=3)

    ax1.plot(xnew, norms[0] * smoothfit(xnew), c='k', ls='--')
    #ax1.plot(xnew, normstwo[0] * smoothfittwo(xnew), c='k', ls='-.')

    #plt.clim(min_, max_)
    scatter1_non = ax1.scatter(eff_rad_non, ir_rad_ratio_non, c=ir_sfr_non, cmap='plasma', marker='o', zorder=2)
    error1_non = ax1.errorbar(eff_rad_non, ir_rad_ratio_non, lolims=True, yerr=0.15, ecolor='k', ls='none', zorder=1)
    #plt.clim(min_, max_)
    ax1.set_xlabel('$\mathrm{R_{e}}$ (kpc)', size=12)

    ax1.scatter([0.29, 0.83], [2.42, 1.54], edgecolors='k', marker='D', s=60, label='', facecolor='white')

    ax1.text(1.5, 1.5, (r'$\frac{L_{\mathrm{mid-IR}}}{L_{\mathrm{IR}}}$(R)' + '\n(Groves+08)'))


    errorbar_avg = ax1.errorbar(1.7, 2.9, avg_ratio_err, 0, ecolor='gray', elinewidth=1., capsize=2.)
    # plt.legend(errorbar_avg, 'Average Uncertainty')

    fig.text(0.08, 0.5, '$\mathrm{SFR_{IR}/SFR_{1.5 \ \mathrm{GHz}}}$', ha="center", va="center", rotation=90, size=12)
    #cbar = plt.colorbar(scatter)
    #cbar.set_label(r'$\mathrm{\overline{SFR}}$')



    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes("right", size="5%", pad=0.05)

    divider2 = make_axes_locatable(ax2)
    cax2 = divider2.append_axes("right", size="5%", pad=0.05)

    # Create and remove the colorbar for the first subplot
    cbar1 = fig.colorbar(scatter1_detect, cax=cax1)
    fig.delaxes(fig.axes[2])

    # Create second colorbar
    cbar2 = fig.colorbar(scatter2_detect, cax=cax2)
    cbar2.set_label(r'SFR$_{\mathrm{IR}}$', size=12)

    plt.subplots_adjust(wspace=0.0)


    plt.savefig('ratio_v_compact.pdf', overwrite=True, dpi=400, bbox_inches='tight')

    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close('all')




    """plt.figure(1413)
    plt.style.use('default')
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"


    sfr_type = 'ir'
    if sfr_type == 'ir':
        x = np.array(tot_table[sfr_to_use]) / (np.pi * np.square(np.array(tot_table['Re'])))
        plt.xlabel(r'$\mathrm{\frac{\Sigma_{SFR_{IR}}}{M_{\odot} yr^{-1} kpc^{-2}}}$')

    else:
        x = np.array(average_sfr / (np.pi * np.square(np.array(tot_table['Re']))))
        plt.xlabel(r'$\mathrm{\frac{<\Sigma_{\overline{SFR}}>}{M_{\odot} yr^{-1} kpc^{-2}}}$')
    plt.xscale('log')

    scatter = plt.scatter(x, ir_rad_ratio, c=average_sfr, cmap='jet')
    plt.ylabel('$\mathrm{SFR_{IR}/SFR_{Radio}}$')
    cbar = plt.colorbar(scatter)
    cbar.set_label(r'$\mathrm{\overline{SFR}}$')

    errorbar_avg = plt.errorbar(1.5, 3, avg_ratio_err, 0, ecolor='gray', elinewidth=1., capsize=2.)
    plt.savefig('ratio_v_sigma.png', overwrite=True, dpi=400, bbox_inches='tight')

    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()

    plt.figure(1414)
    plt.style.use('default')
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"""""




def ratio_age():
    tot_table = pd.read_csv('table.csv')
    idxs = np.where(np.array(tot_table['21 cm SFR']) > 1000.)[0]
    tot_table.drop(idxs, inplace=True)
    tot_table.reset_index(drop=True, inplace=True)
    gal_names = np.array(tot_table['Name'])
    zs = np.array(tot_table['Z'])

    ages = np.array(tot_table['LW_Age'])
    ratio = np.array(tot_table[sfr_to_use]) / np.array(tot_table['21 cm SFR'])

    #binning(ages, ratio)



    ages_detect = np.array(t_ok['LW_Age'])
    ages_non = np.array(t_nondetect['LW_Age'])

    ir_rad_ratio_detect = np.array(t_ok[sfr_to_use]) / np.array(t_ok['21 cm SFR'])

    ir_rad_ratio_non = np.array(t_nondetect[sfr_to_use]) / np.array(t_nondetect['21 cm SFR'])

    ir_sfr_detect = np.array(t_ok[sfr_to_use])
    ir_sfr_non = np.array(t_nondetect[sfr_to_use])

    ratio_err = np.array(tot_table[sfr_to_use]) / np.array(tot_table['21 cm SFR']) * np.sqrt(
        np.square(tot_table[err_to_use] / tot_table[sfr_to_use]) + np.square(
            tot_table['21 cm SFR Error (stat.)'] / tot_table['21 cm SFR']))
    avg_ratio_err = np.nanmean(ratio_err)



    plt.figure(41)

    zs = np.concatenate([ir_sfr_non, ir_sfr_detect], axis=0)
    min_, max_ = zs.min(), zs.max()

    scatter_detect = plt.scatter(ages_detect, ir_rad_ratio_detect, c=ir_sfr_detect, cmap='plasma', zorder=3)
    plt.clim(min_, max_)
    scatter_non = plt.scatter(ages_non, ir_rad_ratio_non, c=ir_sfr_non, cmap='plasma', marker='o', zorder=2)
    plt.clim(min_, max_)
    scatter_non = plt.errorbar(ages_non, ir_rad_ratio_non, lolims=True, yerr=0.15, ecolor='k', ls='none', zorder=1)


    plt.ylabel('$\mathrm{SFR_{IR}/SFR_{1.5\mathrm{GHz}}}$', size=14)
    plt.xlabel('Mean Stellar Age (Myr)', size=14)
    plt.xscale('log')
    cbar = plt.colorbar()
    cbar.set_label(r'SFR$_{\mathrm{IR}}$', size=14)
    plt.errorbar(300, 2.5, yerr=avg_ratio_err, ecolor='gray', elinewidth=1., capsize=2.)
    plt.scatter([23, 89], [2.42, 1.02], edgecolors='k', marker='D', s=60, facecolor='white', label='Binned Median SFR Ratio')
    #plt.legend()

    plt.savefig('ratio_v_age.pdf', overwrite=True, dpi=400, bbox_inches='tight')



    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close('all')

def age_radius():
    tot_table = pd.read_csv('table.csv')
    idxs = np.where(np.array(tot_table['21 cm SFR']) > 1000.)[0]
    tot_table.drop(idxs, inplace=True)
    tot_table.reset_index(drop=True, inplace=True)
    gal_names = np.array(tot_table['Name'])
    zs = np.array(tot_table['Z'])
    ages = np.array(tot_table['LW_Age'])
    eff_rad = np.array(tot_table['Re'])
    plt.figure(55)
    plt.scatter(ages, eff_rad, c='k')
    plt.xlabel('Mean Stellar Age (Myr)')
    plt.ylabel('$R_{e}$')
    plt.savefig('R_v_age.png', overwrite=True, dpi=300, bbox_inches='tight')

def ratio_to_size_condon():
    tot_table = pd.read_csv('table.csv')
    idxs = np.where(np.array(tot_table['21 cm SFR']) > 1000.)[0]
    tot_table.drop(idxs, inplace=True)
    tot_table.reset_index(drop=True, inplace=True)
    gal_names = np.array(tot_table['Name'])
    zs = np.array(tot_table['Z'])
    spec_luminosity = np.array(tot_table['Luminosity'])*u.erg
    rad_luminosity = ((1.5e9*u.Hz)*spec_luminosity).to('W')

    ir_luminosity = (np.array(tot_table['IR Luminosity'])*u.solLum).to('W')




    ir_rad_ratio = np.log10(ir_luminosity/rad_luminosity)
    eff_rad = np.array(tot_table['Re'])

    #ratio_err = np.array(tot_table[sfr_to_use])/np.array(tot_table['21 cm SFR'])*np.sqrt(np.square(tot_table[err_to_use]/tot_table[sfr_to_use])+np.square(tot_table['21 cm SFR Error (stat.)']/tot_table['21 cm SFR']))
    #avg_ratio_err = np.average(ratio_err)

    plt.figure(1412)
    plt.style.use('default')
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"



    plt.xlabel('$\mathrm{R_{e}}$ (kpc)')
    scatter = plt.scatter(eff_rad, ir_rad_ratio, c='b')
    #errorbar_avg = plt.errorbar(1.5, 3, avg_ratio_err, 0, ecolor='gray', elinewidth=1., capsize=2.)
    #plt.legend(errorbar_avg, 'Average Uncertainty')

    plt.ylabel('$\mathrm{log(L_{IR}/L_{Radio})}$')

    IR_condon_table = fits.open('Condon91/Condon_1991_IR.fit')[1].data
    radio_condon_table = Table(fits.open('Condon91/Condon_1991_radio.fit')[1].data)


    def list_duplicates(seq):
        seen = set()
        seen_add = seen.add
        return [idx for idx, item in enumerate(seq) if item in seen or seen_add(item)]

    dupes = list_duplicates(radio_condon_table['Name'])
    radio_condon_table.remove_rows(dupes)



    condon_LFIR = ((10**(np.array(IR_condon_table['logLFIR'])))*u.solLum).to('W')
    distances = (np.array(IR_condon_table['Dist'])*u.Mpc).to('m')
    semi_major_size = ((((np.array(radio_condon_table['FWHMmaj']) + np.array(radio_condon_table['FWHMmin']))/2)*u.arcsec).to('rad')).value
    radio_radii = (distances*semi_major_size/2).to('kpc')
    radio_fluxes = np.array(radio_condon_table['S1_49GHz'])/1000.*u.Jy
    radio_lums = (radio_fluxes*(4*np.pi*distances**2)*(1.49e9*u.Hz)).to('W')
    condon_ratio = np.log10(condon_LFIR/radio_lums)

    cond = plt.scatter(radio_radii, condon_ratio, c='r', label='Condon+91')
    plt.xscale('log')
    plt.scatter(radio_radii[len(radio_radii)-1], condon_ratio[len(condon_ratio)-1], c='g', label='Mrk231')
    plt.legend()



    plt.savefig('ratio_v_size_condon.png', overwrite=True, dpi=400, bbox_inches='tight')

    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()

    def radio_FIR_corr(rad_L):
        return(np.log10(3.826e26*(10**((np.log10(rad_L)-8.66)/1.29))/(rad_L*1.5e9)))



    plt.figure(46)

    fig, ax = plt.subplots()
    plt.hist([condon_ratio, ir_rad_ratio], color=['r', 'b'], label=['Condon+91', 'This paper'], hatch='/', histtype='step', fill=False)
    plt.xlabel(r'$\mathrm{log(L_{IR}/\nu L_{1.5 GHz})}$')

    bars = ax.patches
    patterns = ['/', '']  # set hatch patterns in the correct order

    hatches = []  # list for hatches in the order of the bars
    for h in patterns:  # loop over patterns to create bar-ordered hatches
        for i in range(int(len(bars) / len(patterns))):
            hatches.append(h)
    for bar, hatch in zip(bars, hatches):  # loop over bars and hatches to set hatches in correct order
        bar.set_hatch(hatch)
    #bars[0].set_fill(True)
    #bars[0].set_alpha(0.2)

    plt.text(radio_FIR_corr(1e23)-0.07, 8, r'FIR-Radio Correlation at log(L$_{1.5 \mathrm{GHz}}$)=23', rotation=90, verticalalignment='center')

    plt.text(6.2, 12, 'Local Compact Starbursts \n (Condon+91)', color='darkred')
    plt.text(6.3, 5, 'Intermediate-z\nCompact Starbursts \n (This sample)', color='darkblue')

    plt.axvline(radio_FIR_corr(1e23), linestyle='--', c='k')
    plt.ylabel('Number')
    #plt.axvline(np.median(condon_ratio), linestyle='--', c='r', label='Median')
    #plt.axvline(np.median(ir_rad_ratio), linestyle='--', c='b', label='Median')

    plt.savefig('LIR_vs_Lrad_hist.pdf', overwrite=True, dpi=400, bbox_inches='tight')
    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()

def q_histogram():
    tot_table = pd.read_csv('table.csv')
    idxs = np.where(np.array(tot_table['21 cm SFR']) > 1000.)[0]
    tot_table.drop(idxs, inplace=True)
    tot_table.reset_index(drop=True, inplace=True)

    qs = np.array(tot_table['q'])
    #print(qs)
    lowlim = np.array(tot_table['detect'])
    kmf.fit(qs, lowlim)
    #print(kmf.median_)

    plt.figure(47)

    condon_qs = np.array([2.67, 2.58, 2.76, 2.51, 3.11, 2.1, 2.92, 2.54, 2.66, 2.24, 2.31, 2.38, 2.81, 2.63, 3.])+0.3
    print(condon_qs)


    fig, ax = plt.subplots()
    plt.hist(condon_qs, bins=5, color='r', label='Condon+91',
             histtype='step', fill=False)
    plt.xlabel('q$_{\mathrm{IR}}$', size=14)

    asurv_heights = [1.056, 2.111, 3.393, 2.644, 4.898, 4.898]
    starts = [2.3, 2.5, 2.7, 2.9, 3.1, 3.3]
    ends = np.array(starts) + 0.2

    plt.bar(starts, asurv_heights, width=0.2, edgecolor='b', fill=False, hatch='/', align='edge')

    plt.text(2.82, 7, 'Local Compact\n Starbursts \n (Condon+91)', color='darkred')
    plt.text(3.13, 5.9, 'Intermediate-z\nCompact Starbursts \n (This sample)', color='darkblue')
    #plt.text(2.28, 6, r'FIR-Radio Correlation (Condon+91)', rotation=90,
     #        verticalalignment='center')
    plt.text(2.59, 9.5, 'IR-Radio Correlation (Bell, 2003)', rotation=90)
    plt.text(2.74, 9, 'SFR$_{\mathrm{IR}}$ = SFR$_{1.4 \ \mathrm{GHz}}$', rotation=90)


    plt.ylim(0, 10)
    #plt.axvline(2.34, linestyle='--', c='k')
    plt.axvline(2.64, linestyle='--', c='k')
    plt.axvline(2.72, linestyle=':', c='k')
    plt.ylabel('Number', size=14)
    # plt.axvline(np.median(condon_ratio), linestyle='--', c='r', label='Median')
    # plt.axvline(np.median(ir_rad_ratio), linestyle='--', c='b', label='Median')

    plt.savefig('q_hist.pdf', overwrite=True, dpi=400, bbox_inches='tight')
    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()

def postageStamps():


    #####################################################################################
    # read in 2D Gaussian fit positions from text files to plot ellipses
    found_source = True
    # make plot taller than wide
    vertical = False
    #####################################################################################

    current_dir = os.getcwd()

    sorted_names = GetGalaxyList.return_galaxy_list()

    num_gals = len(sorted_names)

    if vertical:
        cols = 4
        rows = 5
    else:
        if num_gals / 5. >= 4.:
            cols = 5
        else:
            cols = 4
        rows = 5
    # Create figure
    fig = plt.figure(10, figsize=(cols * 6, rows * 6))

    plt.style.use('default')
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"

    spacing = 0.0

    x_size = (1. - spacing) / cols
    y_size = (1. - spacing) / rows
    x_start = spacing / (cols + 1.)
    y_start = 1. - y_size - spacing / (rows + 1.)
    x_iter = 0
    y_iter = 0

    for z in range(num_gals):
        os.chdir(sorted_names[z])

        if z % cols == 0 and z != 0:
            x_iter = 0
            y_iter = y_iter - (1 - y_start)

        x = x_start + x_iter
        y = y_start + y_iter

        x_iter = x_iter + x_size + x_start

        # if source was detected, write galaxy name in boldface
        with open('text/detect.txt', 'r') as f:
            truth = f.readline()
        if int(truth) == 1:
            weighting = 'heavy'
        else:
            weighting = 'normal'

        # open fits file to get beam stats, rms value
        imgname = '%s.cutout.pbcor.fits' % sorted_names[z]
        hdu = fits.open(imgname)
        bmaj = hdu[0].header['bmaj']
        bmin = hdu[0].header['bmin']
        angle = hdu[0].header['bpa']
        rms = np.std(hdu[0].data)

        f = aplpy.FITSFigure(imgname, figure=fig, subplot=[x, y, x_size, y_size], rasterize=True)

        # hst centroid to center the figure
        with open('text/center_HST.txt', 'r') as f_H:
            lines = f_H.readlines()
            HST_ra = float(lines[0])
            HST_dec = float(lines[1])

        # set positions of 2D gaussian ellipse
        if found_source:
            with open('text/center_radio_wcs.txt') as f_rad:
                lines = f_rad.readlines()
                ra = float(lines[0])
                dec = float(lines[1])
        else:
            with open('text/center_HST.txt', 'r') as f_center:
                lines = f_center.readlines()
                ra = float(lines[0])
                dec = float(lines[1])

        # get beam stats to set size of ellipse
        with open('text/width.txt', 'r') as f_width:
            lines = f_width.readlines()
            maj_width = float(lines[0])
            min_width = float(lines[1])
            PA = float(lines[2])

        # show 2D gaussian ellipse
        f.show_ellipses(ra, dec, min_width / 3600., maj_width / 3600., angle=PA, edgecolor='magenta', linewidth=3)
        # show image

        f.show_grayscale(vmin = -8e-5, vmax= 8e-5)
        # write galaxy name
        f.add_label(0.25, 0.9, ('%s' % sorted_names[z])[:5], relative=True, weight=weighting, color='orange', size=60)
        # recenter to HST centroid
        f.recenter(HST_ra, HST_dec, width=15. / 3600, height=15. / 3600)

        # put beam in upper right corner
        f.add_beam()
        f.beam.set_major(bmaj)
        f.beam.set_minor(bmin)
        f.beam.set_angle(angle)
        f.beam.show()
        f.beam.set_corner('bottom left')
        f.beam.set_color('cyan')
        f.beam.set_frame(True)

        """if z==0:
            f.add_scalebar(5./3600.)
            f.scalebar.show(5./3600.)
            f.scalebar.set_corner('bottom right')
            f.scalebar.set_label('5 arcseconds')
            f.scalebar.set_color('white')
            f.scalebar.set_alpha(1.)
            f.scalebar.set_font(size='xx-large')
            f.scalebar.set_linewidth(6.)"""

        # colorbar
        """f.add_colorbar()
        f.colorbar.show(log_format=False)
        f.colorbar.set_width(0.2)
        f.colorbar.set_pad(-0.7)
        f.colorbar.set_location('bottom')
        f.colorbar.set_ticks([-2 * rms, -1 * rms, 0, 1 * rms, 2 * rms])
        f.colorbar.set_frame_color('white')
        # f.colorbar.set_axis_label_text('Surface Brightness (Jy/beam)')
        # f.colorbar.set_axis_label_pad(9)
        # f.colorbar.set_axis_label_rotation(270)
        f.colorbar.set_font(size='small', weight='medium', stretch='extra-condensed', family='sans-serif',
                            style='normal',
                            variant='normal')
        f.colorbar.set_axis_label_font(size='x-small')"""

        # hide angle coordinates
        f.axis_labels.hide_y()
        f.axis_labels.hide_x()
        f.ticks.hide()
        f.tick_labels.hide()
        # f.axis_labels.set_ypad(0)
        # f.tick_labels.set_font(size='x-small', weight='medium', stretch='normal', family='sans-serif', style='normal',
        # variant='normal')

        os.chdir('..')
    cbar_ax = fig.add_axes([0.1, 0.15, .8, 0.03])
    #fig.colorbar(f, cax=cbar_ax)
    norm = mpl.colors.Normalize(vmin=-80, vmax=80)

    cb = mpl.colorbar.ColorbarBase(cbar_ax, cmap='gray', norm=norm, orientation='horizontal')
    cb.ax.tick_params(labelsize=40)
    cb.set_label('Surface Brightness ($\mu$Jy beam$^{-1}$)', size=40)

    os.chdir(current_dir)
    plt.savefig('Stamps.pdf', bbox_inches='tight')
    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()
    plt.close('all')

def HST_Stamps():

    vertical = False

    current_dir = os.getcwd()

    result = pd.read_csv('table.csv', engine='python')

    names = GetGalaxyList.return_galaxy_list()

    # Sort by RA
    stripped_names = []
    sorted_names = []

    for name in names:
        new_name = name.split('J')[1]
        stripped_names.append(new_name)

    sorted_nums = sorted(stripped_names)

    for thing in sorted_nums:
        sorted_names.append('J' + thing)

    #############
    rank = False
    if rank:
        rank_by_ratio = list(np.array(result['IR SFR']) / np.array(result['21 cm SFR']))
        # sort names by SFR ratio
        sorted_names = [x for _, x in sorted(zip(rank_by_ratio, names))]

    ##############

    num_gals = len(sorted_names)

    if vertical == True:
        cols = 4
        rows = 5
    else:
        if num_gals / 5. >= 4.:
            cols = 5
        else:
            cols = 4
        rows = 5

    # Create figure
    fig = plt.figure(11, figsize=(cols * 6, rows * 6))

    spacing = 0.0

    x_size = (1. - spacing) / cols
    y_size = (1. - spacing) / rows
    x_start = spacing / (cols + 1.)
    y_start = 1. - y_size - spacing / (rows + 1.)
    x_iter = 0
    y_iter = 0

    for z in range(len(sorted_names)):
        os.chdir(sorted_names[z])

        if z % cols == 0 and z != 0:
            x_iter = 0
            y_iter = y_iter - (1 - y_start)

        x = x_start + x_iter
        y = y_start + y_iter

        x_iter = x_iter + x_size + x_start

        with open('text/detect.txt', 'r') as f:
            truth = f.readline()

        if int(truth) == 1:
            weighting = 'bold'
        else:
            weighting = 'normal'

        HST_name = sorted_names[z][:5] + '_HST.fits'

        imgname = '%s.cutout.pbcor.fits' % sorted_names[z]

        with open('text/center_HST.txt', 'r') as f_center:
            lines = f_center.readlines()
            HST_ra = float(lines[0])
            HST_dec = float(lines[1])
        with open('text/stdev.txt', 'r') as f_rms:
            rms = float(f_rms.readline())

        f = aplpy.FITSFigure(HST_name, figure=fig, subplot=[x, y, x_size, y_size], rasterize=True, downsample=4)

        fitsimg = fits.open(HST_name)
        hstrms = np.std(fitsimg[0].data)
        print(hstrms)
        print(np.max(fitsimg[0].data))

        f.show_grayscale(vmin=0, vmax=30)
        f.show_contour(imgname, levels=[(2. * rms), (3. * rms), (6 * rms), (12 * rms)], alpha=1, linewidths=4.,
                       colors='cyan')
        f.recenter(HST_ra, HST_dec, width=15. / 3600, height=15. / 3600)

        # f.add_scalebar(5./3600.)
        # f.scalebar.show(5./3600.)
        # f.scalebar.set_corner('bottom right')
        # f.scalebar.set_label('5 arcseconds')

        f.axis_labels.hide_y()
        f.axis_labels.hide_x()
        f.ticks.hide()
        f.tick_labels.hide()

        f.add_label(0.25, 0.9, ('%s' % sorted_names[z])[:5], relative=True, weight=weighting, color='orange', size=60)
        f.set_theme('publication')

        os.chdir('..')

    cbar_ax = fig.add_axes([0.1, 0.15, .8, 0.03])
    #fig.colorbar(f, cax=cbar_ax)
    norm = mpl.colors.Normalize(vmin=0, vmax=50)

    cb = mpl.colorbar.ColorbarBase(cbar_ax, cmap='gray_r', norm=norm, orientation='horizontal')
    cb.ax.tick_params(labelsize=40)
    cb.set_label('Brightness (Counts)', size=40)

    os.chdir(current_dir)
    fig.subplots_adjust(right=0.8)
    plt.savefig('HST_Stamps.pdf', bbox_inches='tight', pad_inches=0, dpi=200)
    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close('all')


#plot_all_SFRs()
#sfr_v_speed()
#plot_lum_vs_z()
#wise_colors()
ratio_to_size()
ratio_age()
#q_histogram()
#postageStamps()
#HST_Stamps()

