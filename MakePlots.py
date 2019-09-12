import matplotlib.pyplot as plt
import matplotlib.colors as matcol
import seaborn as sns
import glob
import pandas as pd
import matplotlib.transforms as tf
from astropy.table import Table
import numpy as np
from scipy.optimize import curve_fit
import CalcSFRs

reload(CalcSFRs)
import Templates

reload(Templates)

####################################################################################
# parameters

# Keep this true
use_imfit = True
# Use plt.annotate to place the names of the galaxies next to points
label_points = True

marker_size = 9
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


# main plot for paper
# comparing SFRs derived from radio & IR, plus extra axes for comparing luminosities in radio & IR
def plot_all_SFRs():
    x_axis_lim = 1000
    y_axis_lim = 600

    # Making arrays of SFRs and uncertainties for non-AGN detections
    irsfrok = np.array(t_ok['IR SFR']).astype(float)
    irok_uncertainty = np.array(t_ok['IR SFR Err'])
    ##############################
    irok_syst = np.array(t_ok['IR Error syst.']) * irsfrok
    radiosfr_ok = np.array(t_ok['21 cm SFR']).astype(float)
    sfr_ok_uncertainty = np.array(t_ok['21 cm SFR Error (stat.)']).astype(float)
    names = t_ok['Name']
    labels = np.random.randint(0, 3, size=len(irsfrok))

    # for AGN (one deviant galaxy)
    flagIRSFR = np.array(t_bad['IR SFR'])
    flagIR_uncertainty = np.array(t_bad['IR SFR Err'])
    ##################
    flagIR_syst = np.array(t_bad['IR Error syst.']) * flagIRSFR
    flagRadioSFR = np.array(t_bad['21 cm SFR'])
    flag_SFR_uncertainty = np.array(t_bad['21 cm SFR Error (stat.)'])
    flagnames = t_bad['Name']

    # and for non-detections
    ir_non_detect = np.array(t_nondetect['IR SFR'])
    ir_non_unc = np.array(t_nondetect['IR SFR Err'])
    ##################
    ir_non_syst = np.array(t_nondetect['IR Error syst.']) * ir_non_detect
    radio_non = np.array(t_nondetect['21 cm SFR'])
    radio_non_unc = np.array(t_nondetect['21 cm SFR Error (stat.)'])
    non_names = t_nondetect['Name']

    # Get indices of sources which lie below proposed detection limit
    # detect_x_lims = (irsfrok < 30)
    # non_detect_x_lims = (ir_non_detect < 30)

    # Generate one-to-one line
    one_to_one = np.poly1d([1, 0])
    # Do the linear fits to non-AGN detections only
    fits = linear_fit(irsfrok, radiosfr_ok, irok_uncertainty, sfr_ok_uncertainty)

    # Generate subplots for the broken axis effect, ax2 for majority of data and ax for AGN
    fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 10), dpi=300, gridspec_kw={'height_ratios': [1, 4]})

    # Plot all data with different markers. For non-detections, make y-axis upper-limit arrows. For sources below
    # proposed detection limit, make x-axis upper limit arrows.
    ok_syst = ax2.errorbar(irsfrok, radiosfr_ok, yerr=0, xerr=irok_syst, fmt='none', ecolor='violet', elinewidth=1.,
                           capsize=3)
    ok = ax2.errorbar(irsfrok, radiosfr_ok, yerr=sfr_ok_uncertainty, xerr=irok_uncertainty, fmt='o', ecolor='k',
                      markeredgecolor='b', markerfacecolor='b', capsize=2, ms=marker_size, elinewidth=2.)
    flagged_syst = ax.errorbar(flagIRSFR, flagRadioSFR, yerr=0, xerr=flagIR_syst, fmt='none', ecolor='violet',
                               elinewidth=1., capsize=3)
    flagged = ax.errorbar(flagIRSFR, flagRadioSFR, yerr=flag_SFR_uncertainty, xerr=flagIR_uncertainty, fmt='o',
                          ecolor='k', c='r', capsize=2, marker='x', ms=marker_size, elinewidth=2.)
    non_detect_syst = ax2.errorbar(ir_non_detect, radio_non, yerr=0, xerr=ir_non_syst, fmt='none', elinewidth=1.,
                                   ecolor='violet', capsize=3)
    non_detect = ax2.errorbar(ir_non_detect, radio_non, yerr=2 * radio_non / 10, xerr=ir_non_unc, fmt='o', ecolor='k',
                              c='gold', capsize=2, uplims=True, marker='v', ms=marker_size, elinewidth=2.)

    # Plot the linear fits, the one-to-one line, and the detection limit dashed lines

    # fit_line = ax2.plot(np.linspace(0, x_axis_lim), fits[0](np.linspace(0, x_axis_lim)), 'g')
    one_to_one_line = ax2.plot(np.linspace(0, x_axis_lim), one_to_one(np.linspace(0, x_axis_lim)), 'k--')
    fixed_line = ax2.plot(np.linspace(0, x_axis_lim), np.poly1d([1. / 2.5, 0])(np.linspace(0, x_axis_lim)), 'c')
    ir_lim_line = ax.axvline(x=30, color='orange', ls='dashed')
    ax2.vlines(x=30, ymin=30, ymax=y_axis_lim, colors='orange', linestyles='dashed')
    ax2.hlines(y=30, xmin=30, xmax=x_axis_lim, colors='orange', linestyles='dashed')

    # copy of bottom axis to show luminosities
    tempax = ax2.twinx()
    ax3 = tempax.twiny()
    # copies of top axis to show luminosities
    ax4 = ax.twinx()
    ax5 = ax.twiny()

    # Titles
    # plt.suptitle('Star Formation Rate Comparison', y=0.91)
    fig.text(0.05, 0.5, '$\mathrm{SFR_{1.5 GHz} \ (M_{\odot} yr^{-1}})$', va='center', rotation='vertical', fontsize=18)
    fig.text(0.97, 0.5, '$\mathrm{L_{1.5 GHz}}$ (W Hz $^{-1}$)', va='center', rotation='vertical', fontsize=18)
    ax2.set_xlabel('$\mathrm{SFR_{IR}} \ (\mathrm{M}_{\odot} \mathrm{yr}^{-1})$', fontsize=18)
    ax5.set_xlabel('$\mathrm{L_{IR}}$ (W)', fontsize=18)

    # put equations of linear fits on plot
    # plt.annotate('%s' % fits[0], (45, 180))
    # plt.annotate('%s' % fits[1], (50, 80))

    # ratio between radio luminosity (erg/s/Hz) and SFR. Multiplying SFR bounds by this number translates them to lumiosity bounds (W/Hz)
    radio_conversion = 1.123517708e21

    low_sfr_lim = 5.

    # Log scales, axis limits
    ax.set_ylim(flagRadioSFR[0] - 10000, flagRadioSFR[0] + 10000)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax5.set_xscale('log')
    ax2.set_yscale('log')
    ax3.set_yscale('log')
    ax3.set_xscale('log')
    ax4.set_yscale('log')

    ax2.set_ylim(low_sfr_lim, y_axis_lim)
    ax3.set_ylim(low_sfr_lim * radio_conversion, y_axis_lim * radio_conversion)
    ax4.set_ylim((flagRadioSFR[0] - 10000) * radio_conversion, (flagRadioSFR[0] + 10000) * radio_conversion)

    # ratio between IR luminosity and SFR. 3.828E26 is a solar lumionsity in watts, and other factor is ratio between luminostiy in L_sol and SFR
    l_solar = 3.828 * (10 ** 26)
    # multiplying an IR SFR limit by this number converts to a luminosity in Watts
    ir_conversion = 6692265625.0 * l_solar

    ax.set_xlim(low_sfr_lim, x_axis_lim)
    ax2.set_xlim(low_sfr_lim, x_axis_lim)
    ax5.set_xlim(low_sfr_lim * ir_conversion, x_axis_lim * ir_conversion)
    ax3.set_xlim(low_sfr_lim * ir_conversion, x_axis_lim * ir_conversion)

    # star forming population defined by any galaxy with SFG = True
    SFG = t_survey[np.where(np.array(t_survey['SFG']) == 'T')[0]]
    # AGN population defined by any galaxy with EITHER Xray, MIR, or SED flagged AGN
    AGN = t_survey[np.where((np.array(t_survey['XrayAGN']) == 'T') | (np.array(t_survey['MIRAGN']) == 'T') | (
    np.array(t_survey['SEDAGN']) == 'T'))[0]]

    # plot contours of density of points on plot
    SFG_plot = sns.kdeplot(((10 ** SFG['logLTIRSF']) * (l_solar)), (10 ** SFG['logL21cm']), ax=ax3, n_levels=5,
                           cmap='Blues')
    AGN_plot = sns.kdeplot(((10 ** AGN['logLTIRSF']) * (l_solar)), (10 ** AGN['logL21cm']), ax=ax3, n_levels=5,
                           cmap='Greens')

    # Legend
    ax.legend((ok, flagged, non_detect, ok_syst, one_to_one_line[0], fixed_line[0], ir_lim_line),
              ('Detections', 'Radio AGN', '$3\sigma$ Upper Limits', 'Template variation',
               '$\mathrm{SFR}_{\mathrm{IR}} = \mathrm{SFR}_{\mathrm{1.5 GHz}}$',
               r'$\mathrm{SFR}_{\mathrm{IR}} = 2.5 \times \mathrm{SFR}_{\mathrm{1.5 GHz}}$',
               'Proposed Detection Limit'), prop={'size': 8})

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
    ax4.spines['bottom'].set_visible(False)
    ax5.spines['bottom'].set_visible(False)
    ax.tick_params(axis='x', which='both', bottom=False)
    ax4.tick_params(axis='x', which='both', bottom=False)
    ax5.tick_params(axis='x', which='both', bottom=False)
    ax2.spines['top'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.get_xaxis().set_visible(False)
    tempax.spines['top'].set_visible(False)

    ax3.axes.get_xaxis().set_visible(False)

    # ax.xaxis.tick_top()
    # ax.tick_params(labeltop='off')  # don't put tick labels at the top
    # ax2.xaxis.tick_bottom()

    # ax2.set_aspect('equal', adjustable='box')


    """# Trying to get the 2 subplots to line up
    bottom_params = np.array(ax2.get_position())
    x_left = bottom_params[0][0]
    x_right = bottom_params[1][0]

    top_params = np.array(ax.get_position())
    y_top = top_params[1][1]
    y_low = top_params[0][1]

    box = tf.Bbox(np.array([[x_left, y_low], [x_right, y_top]]))
    print(box)

    ax.set_position(box)"""

    if label_points:
        for x in range(len(irsfrok)):
            ax2.annotate((names[x].split('.')[0])[:5], (irsfrok[x], radiosfr_ok[x]), xytext=(0, 2),
                         textcoords='offset points',
                         ha='right', va='bottom')
        for x in range(len(flagIRSFR)):
            ax.annotate((flagnames[x].split('.')[0])[:5], (flagIRSFR[x], flagRadioSFR[x]), xytext=(0, 2),
                        textcoords='offset points',
                        ha='right', va='bottom')

        for x in range(len(ir_non_detect)):
            ax2.annotate((non_names[x].split('.')[0])[:5], (ir_non_detect[x], radio_non[x]), xytext=(0, 2),
                         textcoords='offset points',
                         ha='right', va='bottom')

    plt.savefig('SFR_all_plot.png', overwrite=True, bbox_inches='tight')
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
    # detections
    SFR_ok = np.log10(np.array(t_ok['21 cm SFR'])/(np.pi * np.square(np.array(t_ok['Re']))))
    SFR_uncertainty_ok = np.log10(np.array(t_ok['21 cm SFR Error (stat.)']))
    v_ok = np.array(t_ok['v_out'])
    idxs = np.where(v_ok > 0)[0]
    SFR_ok = SFR_ok[idxs]
    v_ok = v_ok[idxs]
    SFR_uncertainty_ok = SFR_uncertainty_ok[idxs]
    v_uncertainty_ok = 0
    names_ok = t_ok['Name']

    # non-detections
    SFR_non = np.log10(np.array(t_nondetect['21 cm SFR'])/(np.pi * np.square(np.array(t_nondetect['Re']))))
    SFR_uncertainty_non = np.log10(np.array(t_nondetect['21 cm SFR Error (stat.)']))
    v_non = np.array(t_nondetect['v_out'])
    idxs2 = np.where(v_non > 0)[0]
    SFR_non = SFR_non[idxs2]
    v_non = v_non[idxs2]
    v_uncertainty_non = 0
    names_non = t_nondetect['Name']

    fig = plt.figure(45)

    ok = plt.errorbar(SFR_ok, v_ok, xerr=SFR_uncertainty_ok, yerr=v_uncertainty_ok, fmt='o', ecolor='k', capsize=2,
                      c='b', ms=marker_size)

    non = plt.errorbar(SFR_non, v_non, xerr=SFR_non / 5, yerr=v_uncertainty_non, fmt='o', ecolor='k', capsize=2,
                       c='gold', xuplims=True, marker='<', ms=marker_size)

    fig.legend((ok, non), ('Detections', '$3\sigma$ Upper Limits'))

    plt.yscale('log')
    plt.ylabel('Outflow Velocity (km s$^{-1}$)', fontsize=18)
    plt.xlabel(r'log$\left(\mathrm{\frac{SFR_{1.5 GHz}}{M_{\odot} yr^{-1}}} \right)$', fontsize=18)

    if label_points:
        for x in range(len(SFR_ok)):
            plt.annotate((names_ok[x].split('.')[0])[:5], (v_ok[x], SFR_ok[x]), xytext=(0, 2),
                         textcoords='offset points',
                         ha='right', va='bottom')

        for x in range(len(SFR_non)):
            plt.annotate((names_non[x].split('.')[0])[:5], (v_non[x], SFR_non[x]), xytext=(0, 2),
                         textcoords='offset points',
                         ha='right', va='bottom')
    plt.savefig('lum_vs_speed.png', overwrite=True, bbox_inches='tight', dpi=300)
    plt.clf()
    plt.cla()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()


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

    fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 10), dpi=300, gridspec_kw={'height_ratios': [1, 4]})

    ok = ax2.errorbar(z_ok, lum_ok, yerr=lum_uncertainty_ok, xerr=z_uncertainty_ok,
                      fmt='o', ecolor='k', capsize=2, c='b', ms=marker_size)
    bad = ax.errorbar(z_bad, lum_bad, yerr=lum_uncertainty_bad, xerr=z_uncertainty_bad,
                      fmt='o', ecolor='k', capsize=2, c='r', marker='x', ms=marker_size)
    non = ax2.errorbar(z_non, lum_non, yerr=lum_non / 5, xerr=z_uncertainty_non,
                       fmt='o', ecolor='k', capsize=2, c='gold', uplims=True, marker='v', ms=marker_size)

    ax.legend((ok, bad, non), ('Detections', 'Radio AGN', '$3\sigma$ Upper Limits'))

    # plt.suptitle('Luminosity vs. Redshift', y=0.92)
    fig.text(0.05, 0.5, '$\mathrm{L_{1.5GHz}} \ (\mathrm{W \ Hz}^{-1}$)', va='center', rotation='vertical', fontsize=18)

    plt.yscale('log')
    plt.xlabel('$z$', fontsize=18)

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

    ax.yaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))

    if label_points:
        for x in range(len(lum_ok)):
            plt.annotate((names_ok[x].split('.')[0])[:5], (z_ok[x], lum_ok[x]), xytext=(0, 2),
                         textcoords='offset points',
                         ha='right', va='bottom')
        for x in range(len(lum_bad)):
            plt.annotate((names_bad[x].split('.')[0])[:5], (z_bad[x], lum_bad[x]), xytext=(0, 2),
                         textcoords='offset points',
                         ha='right', va='bottom')

        for x in range(len(lum_non)):
            plt.annotate((names_non[x].split('.')[0])[:5], (z_non[x], lum_non[x]), xytext=(0, 2),
                         textcoords='offset points',
                         ha='right', va='bottom')
    plt.savefig('lum_vs_z.png', overwrite=True, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()


def plot_relative():
    irsfrok = np.array(t_ok['IR SFR']).astype(float)
    irok_uncertainty = 0.2 * irsfrok
    y_ok = np.array(t_ok['21 cm SFR']).astype(float) / np.array(t_ok['IR SFR']).astype(float)
    y_ok_uncertainty = np.sqrt((np.array(t_ok['21 cm SFR Error (stat.)']).astype(float) / irsfrok) ** 2 +
                               (np.array(t_ok['21 cm SFR']).astype(float) / (irsfrok ** 2) * irok_uncertainty) ** 2)
    names = t_ok['Name']

    flagIRSFR = np.array(t_bad['IR SFR'])
    flagIR_uncertainty = 0.2 * flagIRSFR
    flag_y = np.array(t_bad['21 cm SFR']) / np.array(t_bad['IR SFR'])
    flag_y_uncertainty = np.sqrt((np.array(t_bad['21 cm SFR Error (stat.)']).astype(float) / flagIRSFR) ** 2 +
                                 (np.array(t_bad['21 cm SFR']).astype(float) / (
                                 flagIRSFR ** 2) * flagIR_uncertainty) ** 2)
    flagnames = t_bad['Name']

    ir_non_detect = np.array(t_nondetect['IR SFR']).astype(float)
    ir_non_unc = 0.2 * ir_non_detect
    radio_non = np.array(t_nondetect['21 cm SFR']).astype(float) / np.array(t_nondetect['IR SFR']).astype(float)
    radio_non_unc = np.array(t_nondetect['21 cm SFR Error (stat.)']) / np.array(t_nondetect['IR SFR'])
    non_names = t_nondetect['Name']

    plt.figure(6, figsize=(15, 12), dpi=300)
    # plt.scatter(irsfr, radiosfr, c='r')
    ok = plt.errorbar(irsfrok, y_ok, yerr=y_ok_uncertainty, xerr=irok_uncertainty, fmt='o', ecolor='k', c='b',
                      capsize=2)
    plt.yscale('log')
    plt.xscale('log')
    flagged = plt.errorbar(flagIRSFR, flag_y, yerr=flag_y_uncertainty, xerr=flagIR_uncertainty, fmt='o',
                           ecolor='k', c='r', capsize=2, marker='x')
    non_detect = plt.errorbar(ir_non_detect, radio_non, yerr=radio_non / 5, xerr=ir_non_unc, fmt='o', ecolor='k',
                              c='gold', capsize=2, uplims=True, marker='v')
    bar = plt.axhspan(0.5, 2.0, color='g', alpha=0.1)

    """plt.xlabel('IR SFR $(M_{\odot} yr^{-1})$')
    plt.ylabel('21cm SFR / IR SFR')
    plt.title('Relative Star Formation Comparison (All Sources)')
    plt.legend((ok, flagged, non_detect, bar), ('Detections', 'AGN', 'Non-Detection Upper Limits', 'Factor of 2 Tolerance'))

    for x in range(len(irsfrok)):
        plt.annotate(names[x].split('.')[0], (irsfrok[x], y_ok[x]), xytext=(0, 2), textcoords='offset points',
                     ha='right', va='bottom')
    for x in range(len(flagIRSFR)):
        plt.annotate(flagnames[x].split('.')[0], (flagIRSFR[x], flag_y[x]), xytext=(0, 2), textcoords='offset points',
                     ha='right', va='bottom')

    for x in range(len(ir_non_detect)):
        plt.annotate(non_names[x].split('.')[0], (ir_non_detect[x], radio_non[x]), xytext=(0, 2), textcoords='offset points',
                     ha='right', va='bottom')"""

    plt.savefig('SFR_relative_plot.png', overwrite=True)
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
    lums_ok = list(np.array(t_ok['Luminosity']) / (10000000.))
    lums_non = list(np.array(t_nondetect['Luminosity']) / (10000000.))
    plt.clf()
    plt.cla()
    plt.close()
    plt.close()
    plt.close()
    plt.close()
    colors = ['b', 'gold']
    fig = plt.hist([lums_ok, lums_non], bins=20, stacked=True, label='colors')
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
    plt.xlabel('$\mathrm{L_{1.5 GHz} \ (W \ Hz}^{-1}$)', fontsize=18)
    plt.ylabel('Number', fontsize=18)
    plt.savefig('lum_hist.png', overwrite=True, bbox_inches='tight', dpi=300)
    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()


def simulate_colors(tems, bands, csv):
    tot_mag_list, template_names = [], []
    # iterate through templates
    for tem in tems:
        # redshift template to z=0.6,
        red_spec = Templates.redshift_spectrum(0.6, tem, False)
        red_waves = np.array(red_spec[4])
        lumi = np.array(red_spec[1])

        normalized = []

        # iterate through WISE bands
        for y in range(len(bands)):
            if csv:
                band = pd.read_csv(bands[y], header=None, engine='python')
            else:
                band = pd.read_csv(bands[y], header=None, delim_whitespace=True, engine='python')
            bandwaves = np.array(band.iloc[:, 0])
            band_response = np.array(band.iloc[:, 1])

            # trim template to same wavelength range as WISE band
            cut = np.where((red_waves >= np.min(bandwaves)) & (red_waves <= np.max(bandwaves)))[0]
            trimmed_y = red_waves[cut]
            trimmed_L = lumi[cut]

            # interpolate template to band wavelengths, multiply by the response at that wavelength
            inter_lum = []
            for z in range(len(bandwaves)):
                inter_lum.append(band_response[z] * (np.interp(bandwaves[z], trimmed_y, trimmed_L)))

            # crude method
            """sum_lum = np.sum(np.array(inter_lum))
            sum_waves = np.sum(np.array(band_response))
            normalized.append(sum_lum/sum_waves)"""

            # integrate template multiplied by response function
            spectrum = [bandwaves, inter_lum]
            interped_again = Templates.interpolate_spec(spectrum, True)
            wise_lums = Templates.integrate_spectrum(interped_again[0], interped_again[1], interped_again[2])

            # integrate wise band
            band_spectrum = [bandwaves, band_response]
            interped_band = Templates.interpolate_spec(band_spectrum, True)
            integrated_band = Templates.integrate_spectrum(interped_band[0], interped_band[1], interped_band[2])
            # divide two
            normalized.append(wise_lums / integrated_band)

        tot_mag_list.append(normalized)

        template_names.append(tem.split('.txt')[0].split('/')[8])

    return tot_mag_list, template_names


def wise_colors():

    # w1-w2 vs w3-w4
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
    plt.close()

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
        simulation = simulate_colors(templates, bandpass, False)
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

    one_two = np.array(t_fine['w1-w2']).astype(float)
    two_three = np.array(t_fine['w2-w3']).astype(float)
    names = t_fine['Name']
    ir_rad_ratio = np.array(t_fine['IR SFR'])/np.array(t_fine['21 cm SFR'])

    fig = plt.figure(82)

    cmap = matcol.ListedColormap(['#CE03CE', '#D242D2', '#D588D5', '#75D77D', '#5AC162', '#45AD4D', '#309A38', '#1E8526', '#13741B', '#0B6011', '#044D0A'])
    boundaries = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
    norm = matcol.BoundaryNorm(boundaries, cmap.N, clip=True)
    colors = plt.scatter(two_three, one_two, c=ir_rad_ratio, cmap=cmap, norm=norm)
    cbar = plt.colorbar(colors)
    cbar.set_label('$\mathrm{SFR_{IR}/SFR_{Radio}}$')

    agn_model = [w_two_three[:4], w_one_two[:4], template_names[:4]]
    comp_model = [w_two_three[4:8], w_one_two[4:8], template_names[4:8]]
    sfg_model = [w_two_three[8:11], w_one_two[8:11], template_names[8:11]]

    agn = plt.scatter(agn_model[0], agn_model[1], c='r', marker='x', s=30)
    comp = plt.scatter(comp_model[0], comp_model[1], c='purple', marker='D', s=30)
    sfg = plt.scatter(sfg_model[0], sfg_model[1], c='b', marker='*', s=60)
    plt.xlabel('W2-W3 (mag)')
    plt.ylabel('W1-W2 (mag)')
    plt.legend((colors, agn, comp, sfg), ('unWISE Measured', 'AGN Templates', 'Composite Templates', 'SFG Templates'))

    """if label_points:
        for x in range(len(one_two)):
            plt.annotate((names[x].split('.')[0])[:5], (two_three[x], one_two[x]), xytext=(0, 2),
                         textcoords='offset points',
                         ha='right', va='bottom')

        #for x in range(len(w_two_three)):
            #plt.annotate(template_names[x], (w_two_three[x], w_one_two[x]), xytext=(0, 2), textcoords='offset points',
             #            ha='right', va='bottom')"""

    plt.savefig('wisecolors2.png', overwrite=True, bbox_inches='tight', dpi=300)

    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()

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
    t_shift = t[np.where(t['IR SFR'] > (1.5 * t['21 cm SFR']))[0]]
    t_else = t[np.where(t['IR SFR'] < (1.5 * t['21 cm SFR']))[0]]

    waves = [.3543, .4770, .6231, .7625, .9134, 3.3680, 4.6180, 12.0820, 22.1940]
    plt.clf()
    plt.cla()
    plt.close()
    plt.close()
    plt.close()
    plt.close()
    fig = plt.figure(274)

    vega_to_ab = [2.699, 3.339, 5.174, 6.620]

    standard = 0
    magtrack = []
    for x in range(len(t_shift['Name'])):

        wise = Table.read('unWISE/%s.fits' % (t_shift['Name'][x][:5]))
        SDSS = Table.read('SDSS/%s.fits' % (t_shift['Name'][x][:5]))
        if x == 0:
            standard = SDSS['I']
        #scale = float(standard)/float(SDSS['I'])
        diff = float(standard) - float(SDSS['I'])

        if np.isnan(wise['w4_mag']):

            mags = ([SDSS['u'], SDSS['g'], SDSS['r'], SDSS['I'], SDSS['z'], (wise['w1_mag'] + vega_to_ab[0]),
                    (wise['w2_mag'] + vega_to_ab[1]), (wise['w3_mag'] + vega_to_ab[2]), 0])
        else:
            mags = ([SDSS['u'], SDSS['g'], SDSS['r'], SDSS['I'], SDSS['z'], (wise['w1_mag'] + vega_to_ab[0]),
                    (wise['w2_mag'] + vega_to_ab[1]), (wise['w3_mag'] + vega_to_ab[2]),
                    (wise['w4_mag'] + vega_to_ab[3])])

        mags = [y + diff for y in mags]
        magtrack.append(list(np.array(mags).astype(float).flatten()))

    median_mag = []
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
    plt.errorbar(waves, median_mag, yerr=shift_std, fmt='none', elinewidth=0.5, color='m', capsize=2)

    template_sed = True
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

            plt.plot(red_waves, mags, c=temp_color, linewidth=0.1)




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
        ir_rad_ratio.append(float(tot_table['IR SFR'][x]/tot_table['21 cm SFR'][x]))
        rad_twelve.append(w_flux/z_flux)
        rad_z.append(radio_fluxes[x]/z_flux)

    rad_twelve = np.log10(np.array(rad_twelve))
    rad_z = np.log10(np.array(rad_z))

    plt.figure(151)
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

    gal_names = np.array(tot_table['Name'])

    ir_rad_ratio = np.array(tot_table['IR SFR'])/np.array(tot_table['21 cm SFR'])
    eff_rad = np.array(tot_table['Re'])


    plt.figure(1412)

    scatter = plt.scatter(eff_rad, ir_rad_ratio, c='k')

    plt.xlabel('$\mathrm{R_{e}}$ (kpc)')
    plt.ylabel('$\mathrm{SFR_{IR}/SFR_{radio}}$')
    plt.savefig('ratio_v_size.png', overwrite=True, dpi=400)
    plt.clf()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()





#plot_all_SFRs()
sfr_v_speed()
#plot_lum_vs_z()
#Lum_hist()
wise_colors()
sed()
radio_to_ir()
ratio_to_size()
