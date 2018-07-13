import matplotlib.pyplot as plt
import matplotlib.colors
from astropy.table import Table
import numpy as np
from scipy.optimize import curve_fit

t = Table.read('table.csv')
t_detect = t[np.where(np.multiply(t['detect_pix'], t['detect_aper']))[0]]
correlation = np.multiply(t['detect_pix'], t['detect_aper'])
t_nondetect = t[np.where(correlation == 0)[0]]
t_ok = t_detect[np.where(t_detect['21 cm SFR'] < 1000)[0]]
t_bad = t_detect[np.where(t_detect['21 cm SFR'] > 1000)[0]]

label_points = True


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
    tmp_popt_cons, _ = curve_fit(func, x_vals, y_vals, bounds=([0, -np.inf], [1e-10, np.inf]), sigma=y_vals)

    # Use fitted function to calculate the total error as a result of x and y errors
    new_err_tot = np.sqrt((y_err) ** 2 + (tmp_popt_cons[1] * x_err) ** 2)

    # Use the total error to fit the data again
    popt_cons, _ = curve_fit(func, x_vals, y_vals, bounds=([0, -np.inf], [1e-10, np.inf]), sigma=new_err_tot)

    held_to_zero = np.poly1d([popt_cons[1], 0])
    print(held_to_zero)




    """reduced_chi_squared = residuals[0] / (float(len(x_vals) - 2))
    print(reduced_chi_squared)"""

    return [fit_fn, held_to_zero]


def plot_all_SFRs():

    irsfrok = np.array(t_ok['IR SFR']).astype(float)
    irok_uncertainty = 0.2*irsfrok
    radiosfr_ok = np.array(t_ok['21 cm SFR']).astype(float)
    sfr_ok_uncertainty = np.array(t_ok['21 cm SFR Error (stat.)']).astype(float)
    names = t_ok['Name']
    labels = np.random.randint(0, 3, size=len(irsfrok))

    flagIRSFR = np.array(t_bad['IR SFR'])
    flagIR_uncertainty = 0.2*flagIRSFR
    flagRadioSFR = np.array(t_bad['21 cm SFR'])
    flag_SFR_uncertainty = np.array(t_bad['21 cm SFR Error (stat.)'])
    flagnames = t_bad['Name']

    ir_non_detect = np.array(t_nondetect['IR SFR'])
    ir_non_unc = 0.2*ir_non_detect
    radio_non = np.array(t_nondetect['21 cm SFR'])
    radio_non_unc = np.array(t_nondetect['21 cm SFR Error (stat.)'])
    non_names = t_nondetect['Name']

    one_to_one = np.poly1d([1, 0])
    fits = linear_fit(irsfrok, radiosfr_ok, irok_uncertainty, sfr_ok_uncertainty)

    plt.figure(1, figsize=(15, 12), dpi=200)
    hax = plt.subplot(111)
    ok = plt.errorbar(irsfrok, radiosfr_ok, yerr=sfr_ok_uncertainty, xerr=irok_uncertainty, fmt='o', ecolor='k', c='b', capsize=2)

    plt.yscale('log')
    plt.xscale('log')

    fit_line = plt.plot(np.linspace(0, 600), fits[0](np.linspace(0, 600)), 'g')
    one_to_one_line = plt.plot(np.linspace(0, 600), one_to_one(np.linspace(0, 600)), 'k--')
    fixed_line = plt.plot(np.linspace(0, 600), fits[1](np.linspace(0, 600)), 'gold')

    flagged = plt.errorbar(flagIRSFR, flagRadioSFR, yerr=flag_SFR_uncertainty, xerr=flagIR_uncertainty, fmt='o', ecolor='k', c='r', capsize=2, marker='x')

    plt.xlabel('IR SFR $(M_{\odot} yr^{-1})$')
    plt.ylabel('21cm SFR $(M_{\odot} yr^{-1})$')
    plt.title('Star Formation Comparison (All Sources)')

    #hax2 = plt.subplot(122)

    non_detect = plt.errorbar(ir_non_detect, radio_non, yerr=10**.1, xerr=ir_non_unc, fmt='o',
                              ecolor='k', c='gold', capsize=2, uplims=True, marker='v')

    #hax.set_position([0, 0, 1, 1])
    #hax2.set_position([0, 0, 1, 1])

    #hax2.set_axis_off()

    plt.legend((ok, flagged, non_detect, fit_line[0], one_to_one_line[0], fixed_line[0]),
               ('Detections', 'Possible AGN', 'Non-Detections', 'Weighted fit', 'One to one', 'Y-int fixed 0'))

    if label_points:
        for x in range(len(irsfrok)):
            plt.annotate(names[x].split('.')[0], (irsfrok[x], radiosfr_ok[x]), xytext=(0, 2), textcoords='offset points',
                         ha='right', va='bottom')
        for x in range(len(flagIRSFR)):
            plt.annotate(flagnames[x].split('.')[0], (flagIRSFR[x], flagRadioSFR[x]), xytext=(0, 2), textcoords='offset points',
                         ha='right', va='bottom')

        for x in range(len(ir_non_detect)):
            plt.annotate(non_names[x].split('.')[0], (ir_non_detect[x], radio_non[x]), xytext=(0, 2), textcoords='offset points',
                         ha='right', va='bottom')

    plt.savefig('SFR_all_plot.png', overwrite=True)
    plt.clf()
    plt.close()


def plot_SFRs_detections():

    irsfr = np.array(t_ok['IR SFR'])
    ir_uncertainty = 0.2*irsfr
    radiosfr = np.array(t_ok['21 cm SFR'])
    sfr_uncertainty = np.array(t_ok['21 cm SFR Error (stat.)'])
    names = t_ok['Name']
    labels = np.random.randint(0, 3, size=len(irsfr))

    flagIRSFR = np.array(t_bad['IR SFR'])
    flagIR_uncertainty = 0.2*flagIRSFR
    flagRadioSFR = np.array(t_bad['21 cm SFR'])
    flag_SFR_uncertainty = np.array(t_bad['21 cm SFR Error (stat.)'])

    flagnames = t_bad['Name']

    fit = np.polyfit(irsfr, radiosfr, 1)
    fit_fn = np.poly1d(fit)
    print(fit_fn)


    def get_cmap(n, name='hsv'):
        '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
        RGB color; the keyword argument name must be a standard mpl colormap name.'''
        return plt.cm.get_cmap(name, n)


    cmap = get_cmap(len(irsfr))


    plt.figure(2)
    #plt.scatter(irsfr, radiosfr, c='r')
    good = plt.errorbar(irsfr, radiosfr, yerr=sfr_uncertainty, xerr=ir_uncertainty, fmt='o', ecolor='k', c='b', capsize=2)
    plt.plot(np.linspace(0, 300), fit_fn(np.linspace(0, 300)))
    plt.yscale('log')
    bad = plt.errorbar(flagIRSFR, flagRadioSFR, yerr=flag_SFR_uncertainty, xerr=flagIR_uncertainty, fmt='o', ecolor='k', c='r', capsize=2)
    plt.xlabel('IR SFR $(M_{\odot} yr^{-1})$')
    plt.ylabel('21cm SFR $(M_{\odot} yr^{-1})$')
    plt.title('Star Formation Comparison (Detected Sources)')
    plt.legend((good, bad), ('Detections', 'Possible AGN'))
    for x in range(len(irsfr)):
        plt.annotate(names[x].split('.')[0], (irsfr[x], radiosfr[x]), xytext=(0, 2), textcoords='offset points',
                     ha='right', va='bottom')
    for x in range(len(flagIRSFR)):
        plt.annotate(flagnames[x].split('.')[0], (flagIRSFR[x], flagRadioSFR[x]), xytext=(0, 2), textcoords='offset points',
                     ha='right', va='bottom')
    plt.savefig('SFR_detections.png', overwrite=True)
    plt.clf()
    plt.close()


def plot_good_SFRs():

    irsfr = np.array(t_ok['IR SFR']).astype(float)
    ir_uncertainty = 0.2*irsfr
    radiosfr = np.array(t_ok['21 cm SFR']).astype(float)
    sfr_uncertainty = np.array(t_ok['21 cm SFR Error (stat.)']).astype(float)
    names = t_ok['Name']
    labels = np.random.randint(0, 3, size=len(irsfr))

    one_to_one = np.poly1d([1, 0])
    fits = linear_fit(irsfr, radiosfr, ir_uncertainty, sfr_uncertainty)

    plt.figure(3)
    plt.clf()
    plt.errorbar(irsfr, radiosfr, yerr=sfr_uncertainty, xerr=ir_uncertainty, fmt='o', ecolor='k', capsize=2, c='b')

    fit_line = plt.plot(np.linspace(0, 600), fits[0](np.linspace(0, 600)), 'g')
    one_to_one_line = plt.plot(np.linspace(0,600), one_to_one(np.linspace(0, 600)), 'k--')
    fixed_line = plt.plot(np.linspace(0, 600), fits[1](np.linspace(0, 600)), 'gold')

    plt.legend((fit_line[0], one_to_one_line[0], fixed_line[0]), ('Weighted fit', 'One to one', 'Y-int fixed to 0'))

    plt.xlabel('IR SFR $(M_{\odot} yr^{-1})$')
    plt.ylabel('21cm SFR $(M_{\odot} yr^{-1})$')
    plt.xlim([0, 600])
    plt.ylim([0, 600])
    plt.title('Star Formation Comparison (Unflagged Sources)')

    if label_points:
        for x in range(len(irsfr)):
            plt.annotate(names[x].split('.')[0], (irsfr[x], radiosfr[x]), xytext=(45, -15), textcoords='offset points',
                         ha='right', va='bottom')

    plt.savefig('SFR_unflagged.png', overwrite=True)
    plt.clf()
    plt.close()


def plot_lum_vs_z():

    lum_ok = np.array(t_ok['Luminosity'])
    lum_uncertainty_ok = np.array(t_ok['Luminosity Error (stat.)'])
    z_ok = np.array(t_ok['Z'])
    z_uncertainty_ok = z_ok*.05
    names_ok = t_ok['Name']

    lum_bad = np.array(t_bad['Luminosity'])
    lum_uncertainty_bad = np.array(t_bad['Luminosity Error (stat.)'])
    z_bad = np.array(t_bad['Z'])
    z_uncertainty_bad = z_bad * .05
    names_bad = t_bad['Name']

    lum_non = np.array(t_nondetect['Luminosity'])
    lum_uncertainty_non = np.array(t_nondetect['Luminosity Error (stat.)'])
    z_non = np.array(t_nondetect['Z'])
    z_uncertainty_non = z_non * .05
    names_non = t_nondetect['Name']

    plt.figure(4, figsize=(15, 12), dpi=200)
    ok = plt.errorbar(z_ok, lum_ok, yerr=lum_uncertainty_ok, xerr=z_uncertainty_ok,
                      fmt='o', ecolor='k', capsize=2, c='b')
    bad = plt.errorbar(z_bad, lum_bad, yerr=lum_uncertainty_bad, xerr=z_uncertainty_bad,
                       fmt='o', ecolor='k', capsize=2, c='r', marker='x')
    non = plt.errorbar(z_non, lum_non, yerr=10**28*np.log10(lum_non), xerr=z_uncertainty_non,
                       fmt='o', ecolor='k', capsize=2, c='gold', uplims=True, marker='v')

    plt.legend((ok, bad, non), ('Detections', 'Possible AGN', 'Non-detections'))

    plt.yscale('log')
    plt.xlabel('Redshift')
    plt.ylabel('Luminosity (erg/s)')
    plt.title('Luminosity vs. Redshift')
    if label_points:
        for x in range(len(lum_ok)):
            plt.annotate(names_ok[x].split('.')[0], (z_ok[x], lum_ok[x]), xytext=(0, 2), textcoords='offset points',
                         ha='right', va='bottom')
        for x in range(len(lum_bad)):
            plt.annotate(names_bad[x].split('.')[0], (z_bad[x], lum_bad[x]), xytext=(0, 2), textcoords='offset points',
                         ha='right', va='bottom')

        for x in range(len(lum_non)):
            plt.annotate(names_non[x].split('.')[0], (z_non[x], lum_non[x]), xytext=(0, 2), textcoords='offset points',
                         ha='right', va='bottom')
    plt.savefig('lum_vs_z.png', overwrite=True)
    plt.clf()
    plt.close()


def plot_flux_vs_z():

    f = np.array(t['21 cm Flux'])
    f_uncertainty = np.array(t['21 cm Flux Error'])
    z = np.array(t['Z'])
    z_uncertainty = z*.05
    names = t['Name']

    plt.figure(5)
    # plt.scatter(irsfr, radiosfr, c='r')
    plt.errorbar(z, f, yerr=f_uncertainty, xerr=z_uncertainty, fmt='o', ecolor='k', capsize=2, c='b')
    plt.yscale('log')
    plt.xlabel('Redshift')
    plt.ylabel('Flux (Jy)')
    plt.title('Flux vs. Redshift')
    #for x in range(len(irsfr)):
        #plt.annotate(names[x].split('.')[0], (irsfr[x], radiosfr[x]), xytext=(0, 2), textcoords='offset points',
                     #ha='right', va='bottom')
    plt.savefig('flux_vs_z.png', overwrite=True)
    plt.clf()
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


    plt.figure(6, figsize=(15, 12), dpi=200)
    # plt.scatter(irsfr, radiosfr, c='r')
    ok = plt.errorbar(irsfrok, y_ok, yerr=y_ok_uncertainty, xerr=irok_uncertainty, fmt='o', ecolor='k', c='b',
                      capsize=2)
    plt.yscale('log')
    plt.xscale('log')
    flagged = plt.errorbar(flagIRSFR, flag_y, yerr=flag_y_uncertainty, xerr=flagIR_uncertainty, fmt='o',
                           ecolor='k', c='r', capsize=2, marker='x')
    non_detect = plt.errorbar(ir_non_detect, radio_non, yerr=np.abs(np.log10(radio_non)), xerr=ir_non_unc, fmt='o', ecolor='k',
                              c='gold', capsize=2, uplims=True, marker='v')
    bar = plt.axhspan(0.5, 2.0, color='g', alpha=0.1)

    plt.xlabel('IR SFR $(M_{\odot} yr^{-1})$')
    plt.ylabel('21cm SFR / IR SFR')
    plt.title('Relative Star Formation Comparison (All Sources)')
    plt.legend((ok, flagged, non_detect, bar), ('Detections', 'Possible AGN', 'Non-Detections', '0.5<y<2'))

    for x in range(len(irsfrok)):
        plt.annotate(names[x].split('.')[0], (irsfrok[x], y_ok[x]), xytext=(0, 2), textcoords='offset points',
                     ha='right', va='bottom')
    for x in range(len(flagIRSFR)):
        plt.annotate(flagnames[x].split('.')[0], (flagIRSFR[x], flag_y[x]), xytext=(0, 2), textcoords='offset points',
                     ha='right', va='bottom')

    for x in range(len(ir_non_detect)):
        plt.annotate(non_names[x].split('.')[0], (ir_non_detect[x], radio_non[x]), xytext=(0, 2), textcoords='offset points',
                     ha='right', va='bottom')

    plt.savefig('SFR_relative_plot.png', overwrite=True)
    plt.clf()
    plt.close()


def plot_all_reversed_axes():

    irsfrok = np.array(t_ok['IR SFR'])
    irok_uncertainty = 0.2 * irsfrok
    radiosfr_ok = np.array(t_ok['21 cm SFR'])
    sfr_ok_uncertainty = np.array(t_ok['21 cm SFR Error (stat.)'])
    names = t_ok['Name']
    labels = np.random.randint(0, 3, size=len(irsfrok))

    flagIRSFR = np.array(t_bad['IR SFR'])
    flagIR_uncertainty = 0.2 * flagIRSFR
    flagRadioSFR = np.array(t_bad['21 cm SFR'])
    flag_SFR_uncertainty = np.array(t_bad['21 cm SFR Error (stat.)'])
    flagnames = t_bad['Name']

    ir_non_detect = np.array(t_nondetect['IR SFR'])
    ir_non_unc = 0.2 * ir_non_detect
    radio_non = np.array(t_nondetect['21 cm SFR'])
    radio_non_unc = np.array(t_nondetect['21 cm SFR Error (stat.)'])
    non_names = t_nondetect['Name']

    fit = np.polyfit(radiosfr_ok, irsfrok, 1)
    fit_fn = np.poly1d(fit)

    plt.figure(7)
    # plt.scatter(irsfr, radiosfr, c='r')
    ok = plt.errorbar(radiosfr_ok, irsfrok, xerr=sfr_ok_uncertainty, yerr=irok_uncertainty, fmt='o', ecolor='k', c='b',
                      capsize=2)
    plt.plot(np.linspace(0, 600), fit_fn(np.linspace(0, 600)))
    plt.xscale('log')
    flagged = plt.errorbar(flagRadioSFR, flagIRSFR, xerr=flag_SFR_uncertainty, yerr=flagIR_uncertainty, fmt='o',
                           ecolor='k', c='r', capsize=2)
    non_detect = plt.errorbar(radio_non, ir_non_detect, xerr=radio_non_unc, yerr=ir_non_unc, fmt='o', ecolor='k',
                              c='y', capsize=2)
    plt.ylabel('IR SFR $(M_{\odot} yr^{-1})$')
    plt.xlabel('21cm SFR $(M_{\odot} yr^{-1})$')
    plt.title('Star Formation Comparison (All Sources)')
    plt.legend((ok, flagged, non_detect), ('Detections', 'Possible AGN', 'Non-Detections'))
    """for x in range(len(irsfrok)):
        plt.annotate(names[x].split('.')[0], (irsfrok[x], radiosfr_ok[x]), xytext=(0, 2), textcoords='offset points',
                     ha='right', va='bottom')
    for x in range(len(flagIRSFR)):
        plt.annotate(flagnames[x].split('.')[0], (flagIRSFR[x], flagRadioSFR[x]), xytext=(0, 2), textcoords='offset points',
                     ha='right', va='bottom')

    for x in range(len(ir_non_detect)):
        plt.annotate(non_names[x].split('.')[0], (ir_non_detect[x], radio_non[x]), xytext=(0, 2), textcoords='offset points',
                     ha='right', va='bottom')"""

    plt.savefig('SFR_all_reversed.png', overwrite=True)
    plt.clf()
    plt.close()


plot_all_SFRs()
#plot_SFRs_detections()
plot_good_SFRs()
plot_lum_vs_z()
#plot_flux_vs_z()
plot_relative()
#plot_all_reversed_axes()

