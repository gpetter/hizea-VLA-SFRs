x_axis_lim = 3000

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

    fig = plt.figure(1, figsize=(15, 12), dpi=300)
    ok = plt.errorbar(irsfrok, radiosfr_ok, yerr=sfr_ok_uncertainty, xerr=irok_uncertainty, fmt='o', ecolor='k', c='b', capsize=2)

    plt.yscale('log')
    plt.xscale('log')

    fit_line = plt.plot(np.linspace(0, x_axis_lim), fits[0](np.linspace(0, x_axis_lim)), 'g')
    one_to_one_line = plt.plot(np.linspace(0, x_axis_lim), one_to_one(np.linspace(0, x_axis_lim)), 'k--')
    fixed_line = plt.plot(np.linspace(0, x_axis_lim), fits[1](np.linspace(0, x_axis_lim)), 'c')
    ir_lim_line = plt.axvline(x=30, color='orange', ls='dashed')

    flagged = plt.errorbar(flagIRSFR, flagRadioSFR, yerr=flag_SFR_uncertainty, xerr=flagIR_uncertainty, fmt='o', ecolor='k', c='r', capsize=2, marker='x')

    plt.xlabel('IR SFR $(M_{\odot} yr^{-1})$')
    plt.ylabel('21cm SFR $(M_{\odot} yr^{-1})$')
    plt.title('Star Formation Comparison (All Sources)')


    non_detect = plt.errorbar(ir_non_detect, radio_non, yerr=2*radio_non/10, xerr=ir_non_unc, fmt='o',
                              ecolor='k', c='gold', capsize=2, uplims=True, marker='v')

    plt.legend((ok, flagged, non_detect, fit_line[0], one_to_one_line[0], fixed_line[0], ir_lim_line),
               ('Detections', 'AGN', 'Non-Detection Upper Limits', 'Weighted Linear Fit', 'One to one', 'Weighted Fit Fixed', 'IR Limit'))

    plt.annotate('%s' % fits[0], (45, 180))
    plt.annotate('%s' % fits[1], (50, 80))

    plt.ylim(10, 35000)
    plt.xlim(10, x_axis_lim)

    plt.gca().set_aspect('equal', adjustable='box')

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