import os
import numpy as np
from random import randint
import matplotlib.pyplot as plt
from photutils import CircularAperture
from photutils import EllipticalAperture
from photutils import EllipticalAnnulus
from photutils import CircularAnnulus
from photutils import aperture_photometry
from astropy.io import fits
from astropy.coordinates import SkyCoord  # High-level coordinates
from photutils import SkyCircularAperture
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
from astropy import wcs
import astropy.units as u
import CalcSFRs
CalcSFRs = reload(CalcSFRs)


#############################################################################
# parameters
cell_size = 0.2  # arcsec/pixel
aperture_size = 4  # arcsec (radius)
bkgd_subtract = False
make_growth_curves = False
#############################################################################


def find_aperture_size(position, img):

    ap_rad_list = np.arange(0.1, 20, 0.1)
    rad_pix = ap_rad_list/cell_size
    f = []
    for x in range(len(ap_rad_list)):
        ap = CircularAperture(position, r=rad_pix[x])
        phot = aperture_photometry(img, ap)
        f.append(phot['aperture_sum'][0])
    maximum = max(f)
    idxs = np.where(f>0.9*maximum)[0]
    plt.figure(randint(0, 10000))
    plt.plot(ap_rad_list, f, 'ro')
    plt.xlabel('Aperture radius (")')
    plt.ylabel('Flux Density * pixels per beam')
    plt.savefig('curveofgrowth.png', overwrite=True)
    plt.close()
    if len(idxs) > 0:
        return(rad_pix[min(idxs)])
    else: return(10)


def photometry(gal_name):

    os.chdir(gal_name)
    fits_name = '%s.cutout.pbcor.fits' % gal_name
    hdu = fits.open(fits_name)
    data = hdu[0].data[0][0]

    #w = wcs.WCS(hdu[0].header)

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
    beams_per_aper = np.pi*aperture_size**2/beamarea
    testsig = np.sqrt(beams_per_aper)*std_dev/2

    if make_growth_curves:
        aper_radius = find_aperture_size(positions, data)

    apertures = CircularAperture(positions, r=aper_radius)
    mask = apertures.to_mask(method='center')[0].to_image((200, 200))

    max_val_in_aperture = max(map(max, data[np.where(mask)[0]]))

    x_max = np.where(data == max_val_in_aperture)[0][0]
    y_max = np.where(data == max_val_in_aperture)[1][0]
    print(x_max, y_max)

    #pix_aper = apertures.to_pixel(w)

    if bkgd_subtract:
        global annuli
        annuli = CircularAnnulus(positions, r_in=20, r_out=30)
        apers = [apertures, annuli]
    else:
        apers = [apertures]

    phot_table = aperture_photometry(data, apers, error=error_image)

    output = []

    output.append((phot_table['aperture_sum']/pix_per_beam)[0])
    output.append((phot_table['aperture_sum_err']/pix_per_beam)[0])
    output.append(std_dev)
    output.append(pix_per_beam)
    output.append(testsig)
    output.append(beams_per_aper)
    output.append(max_val_in_aperture)
    output.append(max_val_in_aperture/std_dev)
    output.append(((phot_table['aperture_sum']/pix_per_beam)[0])/testsig)

    if bkgd_subtract:
        bkg_mean = phot_table['aperture_sum_1'] / annuli.area()
        bkg_sum = bkg_mean * apertures.area()
        final_sum = phot_table['aperture_sum_0'] - bkg_sum
        phot_table['residual_aperture_sum'] = final_sum
        print(phot_table['residual_aperture_sum'])

    os.chdir('..')
    return (output)

