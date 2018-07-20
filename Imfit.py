import os
from astropy.io import fits
from astropy import wcs

# go to DATA directory
os.chdir('/users/gpetter/DATA')

# get list of directories within 'data_v1'
# each directory corresponds to a galaxy
names = os.listdir('data_v1')
os.chdir('data_v1')

paths_to_files, paths_to_dirs = [], []

fix_width = True

def new_imfit():


    for x in range(len(names)):
        os.chdir(names[x])

        hdu_hst = fits.open('%s_HST.fits' % names[x][:5])
        w_hst = wcs.WCS(hdu_hst[0])
        trans_hst = w_hst.all_pix2world(1199., 1199., 0)
        ra, dec = trans_hst[0], trans_hst[1]

        hdu = fits.open('%s.cutout.pbcor.fits' % names[x])
        w = wcs.WCS(hdu[0])
        trans = w.all_world2pix(ra, dec, 1, 1, 0)
        xpix, ypix = trans[0], trans[1]

        bmaj = hdu[0].header['bmaj']
        bmin = hdu[0].header['bmin']



        with open('max.txt', 'r') as f_max:
            maximum = f_max.readline()
        with open('stdev.txt', 'r') as f_rms:
            rms = f_rms.readline()

        if fix_width:
            fix_str = 'xyab'
        else:
            fix_str = 'xy'

        with open('estimates.txt', 'w') as f_est:
            f_est.write('%s, %s, %s, %sdeg, %sdeg, 0deg, %s' % (maximum, xpix, ypix, bmaj, bmin, fix_str))



        with open('run_imfit.py', 'w') as f:
            paths_to_dirs.append(os.getcwd())
            paths_to_files.append(os.path.realpath(f.name))

            f.write("""imfit(imagename='%s.cutout.pbcor', box='85,85,115,115', estimates = 'estimates.txt', 
            logfile = 'imfit.log', append=False, residual='test_residual', 
            model='test_model', dooff=False, rms='%sJy/beam', summary='summary.log')""" % (names[x], rms))

        os.chdir('..')


new_imfit()

os.chdir('/users/gpetter/DATA')
with open('imfitrun', 'w') as f:
    for x in range(len(paths_to_files)):
        f.write("""cd %s; xvfb-run -d casa -r 5.3.0-143 --nogui -c %s\n""" % (paths_to_dirs[x], paths_to_files[x]))

os.chdir('/users/gpetter/PycharmProjects/untitled')