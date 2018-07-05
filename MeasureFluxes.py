#
# Author: Grayson Petter
# Date: 6/28/18
# Description: Given a path to a galaxy directory, this code will do aperture photometry on the central source
# in the image, and return a flux, error on flux, as well as other relevant values.
#
#

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
reload(CalcSFRs)


#############################################################################
# parameters
cell_size = 0.2  # arcsec/pixel
aperture_size = 3  # arcsec (radius)
bkgd_subtract = False
make_growth_curves = False
source_find = True
#############################################################################


# If enabled, this function will test many aperture sizes on a given galaxy and determine the optimal one.
# In the process, a curve of growth plot is generated and saved in each galaxy's directory.
# The optimal aperture size is then returned so the photometry() function can use it.
def find_aperture_size(position, img):

    # list of aperture radii in arcsec
    ap_rad_list = np.arange(0.1, 20, 0.1)
    rad_pix = ap_rad_list/cell_size
    f = []

    # Do photometry with each aperture
    for x in range(len(ap_rad_list)):
        ap = CircularAperture(position, r=rad_pix[x])
        phot = aperture_photometry(img, ap)
        f.append(phot['aperture_sum'][0])

    # Make curve of growth plot
    plt.figure(randint(0, 10000))
    plt.plot(ap_rad_list, f, 'ro')
    plt.xlabel('Aperture radius (")')
    plt.ylabel('Flux Density * pixels per beam')
    plt.savefig('curveofgrowth.png', overwrite=True)
    plt.clf()
    plt.close()

    # Select aperture at which it first reaches 90% of maximum (can tweak this).
    maximum = max(f)
    idxs = np.where(f > 0.9 * maximum)[0]
    if len(idxs) > 0:
        return rad_pix[min(idxs)]
    else:
        return 10


# Do photometry with a circular aperture
# Arguments: path to directory of pertinent galaxy
# Returns: A flux, flux error, image RMS, number of pixels per beam, number of beams in aperture,
# the maximum pixel in aperture, and the relative errors: flux/error, and max/rms.
def photometry(gal_name):

    # open fits file as 2D array
    os.chdir(gal_name)
    fits_name = '%s.cutout.pbcor.fits' % gal_name
    hdu = fits.open(fits_name)
    data = hdu[0].data[0][0]

    # open text file containing the scaled MAD and beam area saved by generate_cleans.statistics()
    with open('stdev.txt', 'r') as f:
        std_dev = float(f.readline())  # Jy/beam
    with open('beamarea.txt', 'r') as f:
        beamarea = float(f.readline())  # arcsec^2

    # create error image to pass to astropy photometry to calculate flux errors
    ones = np.ones_like(data)
    error_image = std_dev*ones

    # Flux density = Sum of pixels in aperture / number of pixels per beam
    # Pixels per beam = angular area of beam / angular area of one pixel
    pix_per_beam = beamarea/(cell_size**2)

    # aperture radius in pixels
    aper_radius = aperture_size/cell_size

    # number of beams in the aperture
    beams_per_aper = np.pi*aperture_size**2/beamarea

    # flux error calculation
    flux_error = np.sqrt(beams_per_aper)*std_dev/2

    # center of image
    positions = [(100, 100)]

    # If you don't know proper aperture size, let the function find the optimal one, create aperture object
    if make_growth_curves:
        aper_radius = find_aperture_size(positions, data)
    apertures = CircularAperture(positions, r=5)

    # make mask image where pixels outside aperture are zero, then find the maximum value in the aperture
    mask = apertures.to_mask(method='center')[0].to_image((201, 201))
    masked_data = np.multiply(data, mask)
    max_val_in_aperture = max(map(max, masked_data))

    # get coordinates of maximum point
    y_max = np.where(data == max_val_in_aperture)[0][0]
    x_max = np.where(data == max_val_in_aperture)[1][0]
    #print(x_max, y_max)
    #print(x_max, y_max)

    if source_find:
        positions = [(x_max, y_max)]
    apertures = CircularAperture(positions, r=aper_radius)

    # Do background subtraction with an annulus if desired
    if bkgd_subtract:
        global annuli
        annuli = CircularAnnulus(positions, r_in=20, r_out=30)
        apers = [apertures, annuli]
    else:
        apers = [apertures]

    # Do the photometry, background subtraction if desired
    photo_table = aperture_photometry(data, apers, error=error_image)

    if bkgd_subtract:
        bkg_mean = photo_table['aperture_sum_1'] / annuli.area()
        bkg_sum = bkg_mean * apertures.area()
        final_sum = photo_table['aperture_sum_0'] - bkg_sum
        photo_table['residual_aperture_sum'] = final_sum
        print(photo_table['residual_aperture_sum'])

    # return all relevant values
    output = []

    flux = (photo_table['aperture_sum']/pix_per_beam)[0]


    output.append(flux)
    output.append(flux_error)
    output.append(std_dev)
    output.append(pix_per_beam)
    output.append(beams_per_aper)
    output.append(max_val_in_aperture)
    output.append(max_val_in_aperture/std_dev)
    output.append(flux/flux_error)

    os.chdir('..')
    return output

