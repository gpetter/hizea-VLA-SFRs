import astropy
import os
import numpy as np
from photutils import CircularAperture
from photutils import EllipticalAperture
from photutils import EllipticalAnnulus
from photutils import CircularAnnulus
from photutils import aperture_photometry
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord  # High-level coordinates
from photutils import SkyCircularAperture
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
from astropy import wcs
import astropy.units as u
import CalcSFRs


t = Table.read('VLAsample.csv')
a = np.empty(len(t))
a[:] = np.nan
b = np.zeros(len(t))
t['data'], t['Flux'], t['Flux_error'], t['Luminosity'], t['Luminosity_error'], t['SFR'], t['SFR_error'], t['detect'], \
t['rms'], t['Npixperbeam'] = b, a, a, a, a, a, a, a, a, a

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

os.chdir('/users/gpetter/DATA')
names = os.listdir('data_v1')
os.chdir('data_v1')

cell_size = 0.2
aperture_size = 4.0
bkgd_subtract = False

def photometry(gal_name):
    os.chdir(gal_name)
    fits_name = '%s.cutout.pbcor.fits' % gal_name
    hdu = fits.open(fits_name)
    data = hdu[0].data[0][0]

    w = wcs.WCS(hdu[0].header)

    with open('stdev.txt', 'r') as f:
        std_dev = float(f.readline())

    ones = np.ones_like(data)
    error_image = std_dev*ones

    with open('beamarea.txt', 'r') as f:
        beamarea = float(f.readline())

    pix_per_beam = beamarea/(cell_size**2)

    positions = [(100, 100)]
    #radii = 4*u.arcsec
    aper_radius = aperture_size/cell_size
    beams_per_aper = np.pi*aperture_size**2/(beamarea)
    print(gal_name, np.sqrt(beams_per_aper)*std_dev)

    apertures = CircularAperture(positions, r=aper_radius)
    #pix_aper = apertures.to_pixel(w)

    if bkgd_subtract:
        global annuli
        annuli = CircularAnnulus(positions, r_in=20, r_out=30)
        apers = [apertures, annuli]
    else:     apers=[apertures]

    phot_table = aperture_photometry(data, apers, error=error_image)

    output = []

    output.append((phot_table['aperture_sum']/pix_per_beam)[0])
    output.append((phot_table['aperture_sum_err']/pix_per_beam)[0])
    output.append(std_dev)
    output.append(pix_per_beam)

    if bkgd_subtract:
        bkg_mean = phot_table['aperture_sum_1'] / annuli.area()
        bkg_sum = bkg_mean * apertures.area()
        final_sum = phot_table['aperture_sum_0'] - bkg_sum
        phot_table['residual_aperture_sum'] = final_sum
        print(phot_table['residual_aperture_sum'])

    os.chdir('..')
    return (output)


for name in names:

    idx = np.where(t['Name'] == name)[0]
    t['data'][idx] = True

    """tmp = str(t['RA (J2000)'][idx][0])
    ra = float(tmp[:-2])
    print(ra)
    tmp2 = str(t['Dec (J2000)'][idx][0])
    dec = float(tmp2[:-2])
    c = SkyCoord(ra, dec, frame='icrs', unit='deg')
    print(c)"""

    flux_measured = photometry(name)

    t['Flux'][idx] = flux_measured[0]
    t['Flux_error'][idx] = flux_measured[1]
    t['rms'][idx] = flux_measured[2]
    t['Npixperbeam'][idx] = flux_measured[3]

    if(flux_measured[0] > 3*flux_measured[1]):
        t['detect'][idx] = True
    else: t['detect'][idx] = False

    params_measured = CalcSFRs.calc_params(flux_measured[0], flux_measured[1], t['Z'][idx], 0.001)
    t['Luminosity'][idx] = params_measured[0]
    t['Luminosity_error'][idx] = params_measured[1]
    t['SFR'][idx] = params_measured[2]
    t['SFR_error'][idx] = params_measured[3]

print(t)


os.chdir('/users/gpetter/PycharmProjects/untitled')
t.write('table', format='ascii', overwrite=True)