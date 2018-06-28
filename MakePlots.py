import matplotlib.pyplot as plt
import matplotlib.colors
from astropy.table import Table
import numpy as np

t = Table.read('detected_table.csv')

irsfr = np.array(t['IR SFR'])
radiosfr = np.array(t['SFR'])
sfr_uncertainty = np.array(t['SFR_error'])
names = t['Name']
labels = np.random.randint(0, 3, size=len(irsfr))

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
plt.errorbar(irsfr, radiosfr, yerr=sfr_uncertainty, fmt='o', ecolor='k')
plt.plot(irsfr, fit_fn(irsfr))
plt.xlabel('IR SFR $(M_{\odot} yr^{-1})$')
plt.ylabel('21cm SFR $(M_{\odot} yr^{-1})$')
for x in range(len(irsfr)):
    plt.annotate(names[x].split('.')[0], (irsfr[x], radiosfr[x]), xytext=(0, 2), textcoords='offset points',
                 ha='right', va='bottom')
plt.savefig('SFRplot.png', overwrite=True)
plt.clf()
plt.close()