import matplotlib.pyplot as plt
import matplotlib.colors
from astropy.table import Table
import numpy as np

t = Table.read('table.csv')
t_detect = t[np.where(np.multiply(t['detect_pix'], t['detect_aper']))[0]]
#t_nondetect = t[np.where(np.multiply(t['detect_pix'], t['detect_aper']))[0] == 0]
t_ok = t_detect[np.where(t_detect['21 cm SFR'] < 1000)[0]]

def plot_detections():

    t_bad = t_detect[np.where(t_detect['21 cm SFR'] > 1000)[0]]

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
    plt.errorbar(irsfr, radiosfr, yerr=sfr_uncertainty, xerr=ir_uncertainty, fmt='o', ecolor='k', c='b', capsize=2)
    plt.plot(np.linspace(0, 300), fit_fn(np.linspace(0, 300)))
    plt.errorbar(flagIRSFR, flagRadioSFR, yerr=flag_SFR_uncertainty, xerr=flagIR_uncertainty, fmt='o', ecolor='k', c='r', capsize=2)
    plt.xlabel('IR SFR $(M_{\odot} yr^{-1})$')
    plt.ylabel('21cm SFR $(M_{\odot} yr^{-1})$')
    plt.title('Detected Sources')
    for x in range(len(irsfr)):
        plt.annotate(names[x].split('.')[0], (irsfr[x], radiosfr[x]), xytext=(0, 2), textcoords='offset points',
                     ha='right', va='bottom')
    for x in range(len(flagIRSFR)):
        plt.annotate(flagnames[x].split('.')[0], (flagIRSFR[x], flagRadioSFR[x]), xytext=(0, 2), textcoords='offset points',
                     ha='right', va='bottom')
    plt.savefig('SFR_detect_plot.png', overwrite=True)
    plt.clf()
    plt.close()

def plot_good_data():
    irsfr = np.array(t_ok['IR SFR'])
    ir_uncertainty = 0.2*irsfr
    radiosfr = np.array(t_ok['21 cm SFR'])
    sfr_uncertainty = np.array(t_ok['21 cm SFR Error (stat.)'])
    names = t_ok['Name']
    labels = np.random.randint(0, 3, size=len(irsfr))

    fit = np.polyfit(irsfr, radiosfr, 1)
    fit_fn = np.poly1d(fit)
    print(fit_fn)

    def get_cmap(n, name='hsv'):
        '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
        RGB color; the keyword argument name must be a standard mpl colormap name.'''
        return plt.cm.get_cmap(name, n)

    cmap = get_cmap(len(irsfr))

    plt.figure(3)
    # plt.scatter(irsfr, radiosfr, c='r')
    plt.errorbar(irsfr, radiosfr, yerr=sfr_uncertainty, xerr=ir_uncertainty, fmt='o', ecolor='k', capsize=2, c='b')
    plt.plot(np.linspace(0, 300), fit_fn(np.linspace(0, 300)))
    plt.xlabel('IR SFR $(M_{\odot} yr^{-1})$')
    plt.ylabel('21cm SFR $(M_{\odot} yr^{-1})$')
    plt.title('Unflagged Sources')
    for x in range(len(irsfr)):
        plt.annotate(names[x].split('.')[0], (irsfr[x], radiosfr[x]), xytext=(0, 2), textcoords='offset points',
                     ha='right', va='bottom')
    plt.savefig('SFR_good_plot.png', overwrite=True)
    plt.clf()
    plt.close()

plot_detections()
plot_good_data()