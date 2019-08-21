import os
from astropy.io import fits
from astropy import wcs
import GetGalaxyList
reload(GetGalaxyList)

current_dir = os.getcwd()

##########################################################################################
# parameters
fix_width = True
find_source = True
allow_resolved = False
if allow_resolved:
    resolved_list = ['J082733.87', 'J082638.41', 'J110702.87', 'J121955.77', 'J122949.83', 'J134136.79', 'J211824.06',
            'J161332.52']
else:
    resolved_list = []
##########################################################################################

# Go to directory, get list of galaxies
names = GetGalaxyList.return_galaxy_list()
paths_to_files, paths_to_dirs = [], []


def new_imfit():

    for name in names:
        os.chdir(name)
	
	
	if find_source:
		with open('text/center_radio.txt', 'r') as f_center:
			lines = f_center.readlines()
			x_pix = float(lines[0])
			y_pix = float(lines[1])
		hdu_vla = fits.open('%s.cutout.pbcor.fits' % name)
		w = wcs.WCS(hdu_vla[0])
		trans = w.all_pix2world(x_pix, y_pix, 1, 1, 0)
		with open('text/center_radio_wcs.txt', 'w') as f:
			f.write('%s\n' % trans[0])
			f.write('%s\n' % trans[1])
	else:
		# Get sky coordinates of centroid of HST image
		hdu_hst = fits.open('%s_HST.fits' % name[:5])
		w_hst = wcs.WCS(hdu_hst[0])
		trans_hst = w_hst.all_pix2world(1199., 1199., 0)
		ra, dec = trans_hst[0], trans_hst[1]

		# Save centroid RA, Dec to file
		with open('text/center_HST.txt', 'w') as f_center:
		    f_center.write('%s\n' % ra)
		    f_center.write('%s\n' % dec)

		# get pixel coordinates corresponding to HST centroid
		hdu = fits.open('%s.cutout.pbcor.fits' % name)
		w = wcs.WCS(hdu[0])
		trans = w.all_world2pix(ra, dec, 1, 1, 0)
		x_pix, y_pix = trans[0], trans[1]

	hdu = fits.open('%s.cutout.pbcor.fits' % name)

        # Get beam size
        bmaj = hdu[0].header['bmaj']
        bmin = hdu[0].header['bmin']
	pa = hdu[0].header['bpa']

        # Retrieve max and rms to feed as estimates to imfit
        with open('text/max.txt', 'r') as f_max:
            maximum = f_max.readline()
        with open('text/stdev.txt', 'r') as f_rms:
            rms = f_rms.readline()

        if fix_width and (name not in resolved_list):
            fix_str = 'xyabp'
            summary_str = 'fixed_summary.log'
        else:
            fix_str = 'xy'
            summary_str = 'summary.log'

        with open('text/estimates.txt', 'w') as f_est:
            f_est.write('%s, %s, %s, %sdeg, %sdeg, %sdeg, %s' % (maximum, x_pix, y_pix, bmaj, bmin, pa, fix_str))

        with open('run_imfit.py', 'w') as f:
            paths_to_dirs.append(os.getcwd())
            paths_to_files.append(os.path.realpath(f.name))

            f.write("""imfit(imagename='%s.cutout.pbcor', box='85,85,115,115', estimates = 'text/estimates.txt', 
            logfile = 'imfit.log', append=False, residual='test_residual', 
            model='test_model', dooff=False, rms='%sJy/beam', summary='summary.log')""" % (name, rms))

        os.chdir('..')


new_imfit()

os.chdir(current_dir)
with open('imfitrun', 'w') as f:
    for x in range(len(paths_to_files)):
        f.write("""cd %s; xvfb-run -d casa -r 5.3.0-143 --nogui -c %s\n""" % (paths_to_dirs[x], paths_to_files[x]))

st = os.stat('imfitrun')
os.chmod('imfitrun', st.st_mode | 0111)
