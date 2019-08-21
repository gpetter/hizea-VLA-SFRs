import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors
import matplotlib.transforms as tf
from astropy.table import Table
import numpy as np
from scipy.optimize import curve_fit
import CalcSFRs
reload(CalcSFRs)

####################################################################################
# parameters

# Keep this true
use_imfit = True
# Use plt.annotate to place the names of the galaxies next to points
label_points = True

marker_size=9
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



# survey data to plot on top of SFR comparison plot
surv_redshift = np.array(t_survey_raw['zbest'])
surv_rad_lum = np.array(t_survey_raw['logL21cm'])
surv_ir_lum = np.array(t_survey_raw['logLTIRSF'])

# multiplying columns, so we can eliminate data points which don't have all 3 (redshift, radio luminosity, and IR luminosity)
multed = np.multiply(np.multiply(surv_redshift, surv_rad_lum), surv_ir_lum)
t_survey_data = t_survey_raw[np.where(multed > 0)[0]]

# applying redshift cut similar to our sample
surv_z = np.array(t_survey_data['zbest'])
t_survey = t_survey_data[np.where((surv_z > .4) & (surv_z < 0.75))]



# Do two weighted linear fits to detections only. Fix the second fit to have a y-intercept of zero
def linear_fit(x_vals, y_vals, x_err, y_err):

    # Do first fit with just y errors
    tmp_fit = np.polyfit(x_vals, y_vals, 1, w=1 / y_err)
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
    radiosfr_ok = np.array(t_ok['21 cm SFR']).astype(float)
    sfr_ok_uncertainty = np.array(t_ok['21 cm SFR Error (stat.)']).astype(float)
    names = t_ok['Name']
    labels = np.random.randint(0, 3, size=len(irsfrok))

    # for AGN (one deviant galaxy)
    flagIRSFR = np.array(t_bad['IR SFR'])
    flagIR_uncertainty = np.array(t_bad['IR SFR Err'])
    flagRadioSFR = np.array(t_bad['21 cm SFR'])
    flag_SFR_uncertainty = np.array(t_bad['21 cm SFR Error (stat.)'])
    flagnames = t_bad['Name']

    # and for non-detections
    ir_non_detect = np.array(t_nondetect['IR SFR'])
    ir_non_unc = np.array(t_nondetect['IR SFR Err'])
    radio_non = np.array(t_nondetect['21 cm SFR'])
    radio_non_unc = np.array(t_nondetect['21 cm SFR Error (stat.)'])
    non_names = t_nondetect['Name']



    # Get indices of sources which lie below proposed detection limit
    #detect_x_lims = (irsfrok < 30)
    #non_detect_x_lims = (ir_non_detect < 30)

    # Generate one-to-one line
    one_to_one = np.poly1d([1, 0])
    # Do the linear fits to non-AGN detections only
    fits = linear_fit(irsfrok, radiosfr_ok, irok_uncertainty, sfr_ok_uncertainty)

    # Generate subplots for the broken axis effect, ax2 for majority of data and ax for AGN
    fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,10), dpi=300, gridspec_kw = {'height_ratios':[1, 4]})

    # Plot all data with different markers. For non-detections, make y-axis upper-limit arrows. For sources below
    # proposed detection limit, make x-axis upper limit arrows.
    ok = ax2.errorbar(irsfrok, radiosfr_ok, yerr=sfr_ok_uncertainty, xerr=irok_uncertainty, fmt='o', ecolor='k', markeredgecolor='b', markerfacecolor='b', capsize=2, ms=marker_size)
    flagged = ax.errorbar(flagIRSFR, flagRadioSFR, yerr=flag_SFR_uncertainty, xerr=flagIR_uncertainty, fmt='o',
                           ecolor='k', c='r', capsize=2, marker='x', ms=marker_size)
    non_detect = ax2.errorbar(ir_non_detect, radio_non, yerr=2 * radio_non / 10, xerr=ir_non_unc, fmt='o',
                              ecolor='k', c='gold', capsize=2, uplims=True, marker='v', ms=marker_size)




    # Plot the linear fits, the one-to-one line, and the detection limit dashed lines

    #fit_line = ax2.plot(np.linspace(0, x_axis_lim), fits[0](np.linspace(0, x_axis_lim)), 'g')
    one_to_one_line = ax2.plot(np.linspace(0, x_axis_lim), one_to_one(np.linspace(0, x_axis_lim)), 'k--')
    fixed_line = ax2.plot(np.linspace(0, x_axis_lim), np.poly1d([1./2.5, 0])(np.linspace(0, x_axis_lim)), 'c')
    ir_lim_line = ax.axvline(x=30, color='orange', ls='dashed')
    ax2.vlines(x=30, ymin=30, ymax=y_axis_lim, colors='orange', linestyles='dashed')
    ax2.hlines(y=30, xmin=30, xmax=x_axis_lim, colors='orange', linestyles='dashed')




    # copy of bottom axis to show luminosities
    ax3 = ax2.twinx().twiny()
    # copies of top axis to show luminosities
    ax4 = ax.twinx()
    ax5 = ax.twiny()



    # Titles
    #plt.suptitle('Star Formation Rate Comparison', y=0.91)
    fig.text(0.05, 0.5, '1.5 GHz SFR $(M_{\odot} yr^{-1})$', va='center', rotation='vertical', fontsize=18)
    fig.text(0.97, 0.5, '1.5 GHz Luminosity (W Hz $^{-1}$)', va='center', rotation='vertical', fontsize=18)
    ax2.set_xlabel('IR SFR $(M_{\odot} yr^{-1})$', fontsize=18)
    ax5.set_xlabel('IR Luminosity (W)', fontsize=18)

    

    # put equations of linear fits on plot
    #plt.annotate('%s' % fits[0], (45, 180))
    #plt.annotate('%s' % fits[1], (50, 80))

    #ratio between radio luminosity (erg/s/Hz) and SFR. Multiplying SFR bounds by this number translates them to lumiosity bounds (W/Hz)
    radio_conversion = 1.123517708e21

    low_sfr_lim = 1.




    # Log scales, axis limits
    ax.set_ylim(flagRadioSFR[0]-10000, flagRadioSFR[0]+10000)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax5.set_xscale('log')
    ax2.set_yscale('log')
    ax3.set_yscale('log')
    ax3.set_xscale('log')
    ax4.set_yscale('log')

    ax2.set_ylim(low_sfr_lim, y_axis_lim)
    ax3.set_ylim(low_sfr_lim*radio_conversion, y_axis_lim*radio_conversion)
    ax4.set_ylim((flagRadioSFR[0]-10000)*radio_conversion, (flagRadioSFR[0]+10000)*radio_conversion)

    # ratio between IR luminosity and SFR. 3.828E26 is a solar lumionsity in watts, and other factor is ratio between luminostiy in L_sol and SFR
    l_solar = 3.828*(10**26)
    # multiplying an IR SFR limit by this number converts to a luminosity in Watts
    ir_conversion = 6692265625.0*l_solar

    ax.set_xlim(low_sfr_lim, x_axis_lim)
    ax2.set_xlim(low_sfr_lim, x_axis_lim)
    ax5.set_xlim(low_sfr_lim*ir_conversion, x_axis_lim*ir_conversion)
    ax3.set_xlim(low_sfr_lim*ir_conversion, x_axis_lim*ir_conversion)





    # star forming population defined by any galaxy with SFG = True
    SFG = t_survey[np.where(np.array(t_survey['SFG'])=='T')[0]]
    # AGN population defined by any galaxy with EITHER Xray, MIR, or SED flagged AGN
    AGN = t_survey[np.where((np.array(t_survey['XrayAGN'])=='T') | (np.array(t_survey['MIRAGN'])=='T') | (np.array(t_survey['SEDAGN'])=='T'))[0]]

    # plot contours of density of points on plot
    SFG_plot = sns.kdeplot(((10**SFG['logLTIRSF'])*(l_solar)), (10**SFG['logL21cm']), ax=ax3, n_levels=5, cmap='Blues')
    AGN_plot = sns.kdeplot(((10**AGN['logLTIRSF'])*(l_solar)), (10**AGN['logL21cm']), ax=ax3, n_levels=5, cmap='Greens')


    # Legend
    ax.legend((ok, flagged, non_detect, one_to_one_line[0], fixed_line[0], ir_lim_line),
               ('Detections', 'AGN', 'Non-Detection Upper Limits', 'One to one',
                'Factor 2.5 Suppressed', 'Proposed Detection Limit'), prop={'size': 8})


    # Hack to make the diagonal hashes on broken axis
    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-3.5*d, 3.5*d), **kwargs)  # top-left diagonal
    ax.plot((1 - d, 1 + d), (-3.5*d, +3.5*d), **kwargs)  # top-right diagonal
    #plt.gca().set_aspect('equal', adjustable='box')

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
    ax2.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.get_xaxis().set_visible(False)
    ax3.axis('off')
    ax3.set_frame_on(False)
    ax3.axes.get_xaxis().set_visible(False)
    ax3.axes.get_yaxis().set_visible(False)

    #ax.xaxis.tick_top()
    #ax.tick_params(labeltop='off')  # don't put tick labels at the top
    #ax2.xaxis.tick_bottom()

    #ax2.set_aspect('equal', adjustable='box')


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
            ax2.annotate((names[x].split('.')[0])[:5], (irsfrok[x], radiosfr_ok[x]), xytext=(0, 2), textcoords='offset points',
                         ha='right', va='bottom')
        for x in range(len(flagIRSFR)):
            ax.annotate((flagnames[x].split('.')[0])[:5], (flagIRSFR[x], flagRadioSFR[x]), xytext=(0, 2), textcoords='offset points',
                         ha='right', va='bottom')

        for x in range(len(ir_non_detect)):
            ax2.annotate((non_names[x].split('.')[0])[:5], (ir_non_detect[x], radio_non[x]), xytext=(0, 2), textcoords='offset points',
                         ha='right', va='bottom')

    plt.savefig('SFR_all_plot.png', overwrite=True, bbox_inches='tight', pad_inches=0)
    plt.clf()
    plt.cla()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()
    plt.clf()
    plt.cla()
    plt.close()

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


def lum_v_speed():


    SFR_ok = np.log10(np.array(t_ok['21 cm SFR']))
    SFR_uncertainty_ok = np.log10(np.array(t_ok['21 cm SFR Error (stat.)']))
    v_ok = np.array(t_ok['v_out'])
    idxs = np.where(v_ok > 0)[0]
    SFR_ok = SFR_ok[idxs]
    v_ok = v_ok[idxs]
    SFR_uncertainty_ok = SFR_uncertainty_ok[idxs]
    v_uncertainty_ok = 0
    names_ok = t_ok['Name']



    SFR_non = np.log10(np.array(t_nondetect['21 cm SFR']))
    SFR_uncertainty_non = np.log10(np.array(t_nondetect['21 cm SFR Error (stat.)']))
    v_non = np.array(t_nondetect['v_out'])
    idxs2 = np.where(v_non > 0)[0]
    SFR_non = SFR_non[idxs2]
    v_non = v_non[idxs2]
    SFR_uncertainty_non = SFR_uncertainty_non[idxs2]
    v_uncertainty_non = 0
    names_non = t_nondetect['Name']

    fig = plt.figure(45)

    ok = plt.errorbar(SFR_ok, v_ok, xerr=SFR_uncertainty_ok, yerr=v_uncertainty_ok, fmt='o', ecolor='k', capsize=2, c='b', ms=marker_size)

    non = plt.errorbar(SFR_non, v_non, xerr=SFR_non/5, yerr=v_uncertainty_non, fmt='o', ecolor='k', capsize=2, c='gold', xuplims=True, marker='<', ms=marker_size)

    fig.legend((ok, non), ('Detections', 'Non-Detection Upper Limits'))

    #plt.suptitle('Luminosity vs. Wind Speed', y=0.92)
    #fig.text(0.5, 0.05, '1.5 GHz Luminosity (erg s$^{-1}$ Hz$^{-1}$)', va='center', rotation='horizontal', fontsize=18)

    plt.yscale('log')
    plt.ylabel('Outflow Velocity (km s$^{-1}$)', fontsize=18)
    plt.xlabel('log(1.5 GHz SFR) ($M_{\odot} yr^{-1}$)', fontsize=18)

   


    if label_points:
        for x in range(len(SFR_ok)):
            plt.annotate((names_ok[x].split('.')[0])[:5], (v_ok[x], SFR_ok[x]), xytext=(0, 2), textcoords='offset points',
                         ha='right', va='bottom')

        for x in range(len(SFR_non)):
            plt.annotate((names_non[x].split('.')[0])[:5], (v_non[x], SFR_non[x]), xytext=(0, 2), textcoords='offset points',
                         ha='right', va='bottom')
    plt.savefig('lum_vs_speed.png', overwrite=True, bbox_inches='tight')
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

    lum_ok = np.array(t_ok['Luminosity'])
    lum_uncertainty_ok = np.array(t_ok['Luminosity Error (stat.)'])
    z_ok = np.array(t_ok['Z'])
    z_uncertainty_ok = 0
    names_ok = t_ok['Name']

    lum_bad = np.array(t_bad['Luminosity'])
    lum_uncertainty_bad = np.array(t_bad['Luminosity Error (stat.)'])
    z_bad = np.array(t_bad['Z'])
    z_uncertainty_bad = 0
    names_bad = t_bad['Name']

    lum_non = np.array(t_nondetect['Luminosity'])
    lum_uncertainty_non = np.array(t_nondetect['Luminosity Error (stat.)'])
    z_non = np.array(t_nondetect['Z'])
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

    #plt.suptitle('Luminosity vs. Redshift', y=0.92)
    fig.text(0.05, 0.5, '1.5 GHz Luminosity (erg s$^{-1}$ Hz$^{-1}$)', va='center', rotation='vertical', fontsize=18)

    plt.yscale('log')
    plt.xlabel('Redshift', fontsize=18)

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
    y_ok = np.array(t_ok['21 cm SFR']).astype(float)/np.array(t_ok['IR SFR']).astype(float)
    y_ok_uncertainty = np.sqrt((np.array(t_ok['21 cm SFR Error (stat.)']).astype(float)/irsfrok)**2+
                               (np.array(t_ok['21 cm SFR']).astype(float)/(irsfrok**2)*irok_uncertainty)**2)
    names = t_ok['Name']


    flagIRSFR = np.array(t_bad['IR SFR'])
    flagIR_uncertainty = 0.2 * flagIRSFR
    flag_y = np.array(t_bad['21 cm SFR'])/np.array(t_bad['IR SFR'])
    flag_y_uncertainty = np.sqrt((np.array(t_bad['21 cm SFR Error (stat.)']).astype(float)/flagIRSFR)**2+
                                 (np.array(t_bad['21 cm SFR']).astype(float)/(flagIRSFR**2)*flagIR_uncertainty)**2)
    flagnames = t_bad['Name']

    ir_non_detect = np.array(t_nondetect['IR SFR']).astype(float)
    ir_non_unc = 0.2 * ir_non_detect
    radio_non = np.array(t_nondetect['21 cm SFR']).astype(float)/np.array(t_nondetect['IR SFR']).astype(float)
    radio_non_unc = np.array(t_nondetect['21 cm SFR Error (stat.)'])/np.array(t_nondetect['IR SFR'])
    non_names = t_nondetect['Name']

    plt.figure(6, figsize=(15, 12), dpi=300)
    # plt.scatter(irsfr, radiosfr, c='r')
    ok = plt.errorbar(irsfrok, y_ok, yerr=y_ok_uncertainty, xerr=irok_uncertainty, fmt='o', ecolor='k', c='b',
                      capsize=2)
    plt.yscale('log')
    plt.xscale('log')
    flagged = plt.errorbar(flagIRSFR, flag_y, yerr=flag_y_uncertainty, xerr=flagIR_uncertainty, fmt='o',
                           ecolor='k', c='r', capsize=2, marker='x')
    non_detect = plt.errorbar(ir_non_detect, radio_non, yerr=radio_non/5, xerr=ir_non_unc, fmt='o', ecolor='k',
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
	lums_ok = list(np.array(t_ok['Luminosity'])/(10000000.))
	lums_non = list(np.array(t_nondetect['Luminosity'])/(10000000.))
	plt.clf()
	plt.cla()
	plt.close()
	plt.close()
	plt.close()
	plt.close()
	colors = ['b', 'gold']
	fig = plt.hist([lums_ok, lums_non], bins=20, stacked=True, label='colors')
	plt.xlabel('1.5 GHz Luminosity (W Hz $^{-1}$)', fontsize=18)
	plt.ylabel('Counts', fontsize=18)
	plt.savefig('lum_hist.png', overwrite=True)
	plt.clf()
	plt.close()
	plt.clf()
	plt.cla()
	plt.close()

	


plot_all_SFRs()
lum_v_speed()
plot_lum_vs_z()
Lum_hist()
#plot_relative()
