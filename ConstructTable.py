import os
import numpy as np
from astropy.table import Table
import CalcSFRs
CalcSFRs = reload(CalcSFRs)
import MeasureFluxes
MeasureFluxes = reload(MeasureFluxes)

t = Table.read('VLAsample.csv')
a = np.empty(len(t))
a[:] = np.nan
b = np.zeros(len(t))
t['data'], t['Flux'], t['Flux_error'], t['test_error'], t['Luminosity'], t['Luminosity_error'], t['SFR'], t['SFR_error'], \
t['detect_aper'], t['detect_pix'], t['rms'], t['Npixperbeam'], t['Nbeams'], t['MaxValue'], t['Max/noise'], \
t['Flux/error'] = b, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a

# Defining units of each column
t['RA (J2000)'].unit = 'deg'
t['Dec (J2000)'].unit = 'deg'
t['IR SFR'].unit = 'solMass/yr'
t['Flux'].unit = 'Jy'
t['Flux_error'].unit = 'Jy'
t['Luminosity'].unit = 'erg/s'
t['Luminosity_error'].unit = 'erg/s'
t['SFR'].unit = 'solMass/yr'
t['SFR_error'].unit = 'solMass/yr'
t['MaxValue'].unit = 'Jy/beam'

os.chdir('/users/gpetter/DATA')
names = os.listdir('data_v1')
os.chdir('data_v1')

cell_size = 0.2
aperture_size = 4.0
bkgd_subtract = False

for name in names:

    idx = np.where(t['Name'] == name)[0]
    t['data'][idx] = True

    #tmp = str(t['RA (J2000)'][idx][0])
    #ra = float(tmp[:-2])
    #print(ra)
    #tmp2 = str(t['Dec (J2000)'][idx][0])
    #dec = float(tmp2[:-2])
    #c = SkyCoord(ra, dec, frame='icrs', unit='deg')
    #print(c)

    flux_measured = MeasureFluxes.photometry(name)

    t['Flux'][idx] = flux_measured[0]
    t['Flux_error'][idx] = flux_measured[1]
    t['rms'][idx] = flux_measured[2]
    t['Npixperbeam'][idx] = flux_measured[3]
    t['test_error'][idx] = float(flux_measured[4])
    t['Nbeams'][idx] = flux_measured[5]
    t['MaxValue'][idx] = flux_measured[6]
    t['Max/noise'][idx] = flux_measured[7]
    t['Flux/error'][idx] = flux_measured[8]

    if flux_measured[0] > 3*flux_measured[2]:
        t['detect_aper'][idx] = True
    else: t['detect_aper'][idx] = False

    if flux_measured[6] > 3*flux_measured[2]:
        t['detect_pix'][idx] = True
    else: t['detect_pix'][idx] = False

    params_measured = CalcSFRs.calcparams(flux_measured[0], flux_measured[1], t['Z'][idx], 0.001)
    t['Luminosity'][idx] = params_measured[0]
    t['Luminosity_error'][idx] = params_measured[1]
    t['SFR'][idx] = params_measured[2]
    t['SFR_error'][idx] = params_measured[3]


t_data = t[np.where(t['data'])[0]]
print(t_data)


os.chdir('/users/gpetter/PycharmProjects/untitled')
t_data.write('table.csv', format='csv', overwrite=True)
