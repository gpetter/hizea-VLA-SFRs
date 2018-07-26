import os
from astropy.io import fits
from astropy import wcs

current_dir = os.getcwd()

##########################################################################################
# parameters
data_path = '/users/gpetter/DATA/data_v1'
fix_width = True
##########################################################################################

# Go to directory, get list of galaxies
names = os.listdir(data_path)
os.chdir(data_path)

paths_to_files, paths_to_dirs = [], []


def new_imfit():

    for name in names:
        os.chdir(name)

        # Get sky coordinates of centroid of HST image
        hdu_hst = fits.open('%s_HST.fits' % name[:5])
        w_hst = wcs.WCS(hdu_hst[0])
        trans_hst = w_hst.all_pix2world(1199., 1199., 0)
        ra, dec = trans_hst[0], trans_hst[1]

        # Save centroid RA, Dec to file
        with open('center_HST.txt', 'w') as f_center:
            f_center.write('%s\n' % ra)
            f_center.write('%s\n' % dec)

        # get pixel coordinates corresponding to HST centroid
        hdu = fits.open('%s.cutout.pbcor.fits' % name)
        w = wcs.WCS(hdu[0])
        trans = w.all_world2pix(ra, dec, 1, 1, 0)
        x_pix, y_pix = trans[0], trans[1]
        print(x_pix, y_pix)

        # Get beam size
        bmaj = hdu[0].header['bmaj']
        bmin = hdu[0].header['bmin']

        # Retrieve max and rms to feed as estimates to imfit
        with open('max.txt', 'r') as f_max:
            maximum = f_max.readline()
        with open('stdev.txt', 'r') as f_rms:
            rms = f_rms.readline()

        if fix_width:
            fix_str = 'xyab'
            summary_str = 'fixed_summary.log'
        else:
            fix_str = 'xy'
            summary_str = 'summary.log'

        with open('estimates.txt', 'w') as f_est:
            f_est.write('%s, %s, %s, %sdeg, %sdeg, 0deg, %s' % (maximum, x_pix, y_pix, bmaj, bmin, fix_str))

        with open('run_imfit.py', 'w') as f:
            paths_to_dirs.append(os.getcwd())
            paths_to_files.append(os.path.realpath(f.name))

            f.write("""imfit(imagename='%s.cutout.pbcor', box='85,85,115,115', estimates = 'estimates.txt', 
            logfile = 'imfit.log', append=False, residual='test_residual', 
            model='test_model', dooff=False, rms='%sJy/beam', summary='summary.log')""" % (name, rms))

        os.chdir('..')


new_imfit()

os.chdir('/users/gpetter/DATA')
with open('imfitrun', 'w') as f:
    for x in range(len(paths_to_files)):
        f.write("""cd %s; xvfb-run -d casa -r 5.3.0-143 --nogui -c %s\n""" % (paths_to_dirs[x], paths_to_files[x]))

os.chdir(current_dir)