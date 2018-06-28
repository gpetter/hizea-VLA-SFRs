import os

vises = []

# go to DATA directory
os.chdir('/users/gpetter/DATA')

# get list of directories within 'data_v1'
# each directory corresponds to a galaxy
names = os.listdir('data_v1')
os.chdir('data_v1')

# for each galaxy, find all .ms files and append to visibilities list
for name in names:
    all_dirs = os.listdir(name)
    only_ms = [y for y in all_dirs if '.ms' in y]
    vises.append(only_ms)

paths_to_files, paths_to_dirs = [], []
dirty_has_run = False
clean_has_run = False
pb_has_run = False
cutout_has_run = False


############################################################################
# parameters
clean_sigma = 3
img_size = 12000
cutout_size = 200
cut_frame = [img_size/2-cutout_size/2, img_size/2+cutout_size/2]
############################################################################


# Creates a python script in each galaxy's directory
# This script will construct a dirty image, or a very lightly cleaned image (niter~=100)
# Measures the MAD of that image, scales it by 1.4826, and multiplies it by some number to set the threshold
# In the process, a PSF and residual is saved so it does not have to be recalculated when the images are cleaned
def make_dirty_images():

    global dirty_has_run

    for x in range(len(names)):

        os.chdir(names[x])

        with open('run_tclean_%s.py' %(names[x].split('.')[0]), 'w') as f:

            # adding paths to simplify creating pipeline script later
            paths_to_dirs.append(os.getcwd())
            paths_to_files.append(os.path.realpath(f.name))

            # write out dirty tclean command
            f.write("""tclean(vis=%s, imagename='%s', field='0', datacolumn='data',
                verbose=True, gridder='wproject', wprojplanes=128, pblimit=-1, robust=0.5, imsize=[%s], 
                cell='0.2arcsec', specmode='mfs', deconvolver='mtmfs', nterms=2, scales=[0,11,28], 
                interactive=False, niter=50,
                weighting='briggs', usemask='auto-multithresh', stokes='I', threshold='0.0Jy', calcpsf=True,
                calcres=True, savemodel='modelcolumn', restart=False) \n \n""" % (vises[x], names[x], img_size))

            # will call imstat to measure the MAD of each image, scaled by number*1.4826*MAD
            # saves threshold value to text file for CleanImages() script's access
            f.write("""stats=imstat('%s.image.tt0')\n""" %(names[x]))
            f.write("""thresh=%s*1.4826*stats['medabsdevmed'][0]\n""" % clean_sigma)
            f.write("""print(thresh)\n""")
            f.write("""with open('threshold.txt', 'w') as f:\n""")
            f.write("""\tf.write('%s' %(thresh))\n \n""")

        os.chdir('..')
    dirty_has_run = True


# Generates a script which will run tclean with the threshold determined in the dirty_image() step
def clean_images():

    global dirty_has_run, clean_has_run

    if dirty_has_run:
        edit_mode = 'a'
    else:
        edit_mode = 'w'

    for x in range(len(names)):

        os.chdir(names[x])

        with open('run_tclean_%s.py' %(names[x].split('.')[0]), edit_mode) as f:
            if not dirty_has_run:
                paths_to_dirs.append(os.getcwd())
                paths_to_files.append(os.path.realpath(f.name))

            # retrieve previously saved threshold value
            f.write("""with open('threshold.txt', 'r') as f:\n""")
            f.write("""\tglobal threshold\n""")
            f.write("""\tlines=f.readlines()\n""")
            f.write("""\tthreshold=float(lines[0])\n""")
            f.write("""print(threshold)\n \n""")

            # run tclean with said threshold, high niter value
            f.write(("""tclean(vis=%s, imagename='%s', field='0', datacolumn='data',
                   verbose=True, gridder='wproject', wprojplanes=128, pblimit=-1, robust=0.5, imsize=[%s], 
                   cell='0.2arcsec', specmode='mfs', deconvolver='mtmfs', nterms=2, scales=[0,11,28], """
                    % (vises[x], names[x], img_size)) + """ interactive=False, niter=20000, weighting='briggs',
                   usemask='auto-multithresh', sidelobethreshold = 4.0, stokes='I', threshold='%sJy' %(threshold), 
                   savemodel='modelcolumn', calcres=False, calcpsf=False, restart=True) \n \n""")

        os.chdir('..')
    clean_has_run = True


def pb_cor():

    global clean_has_run, pb_has_run
    if clean_has_run:
        editmode='a'
    else: editmode='w'
    for x in range(len(names)):
        os.chdir(names[x])
        with open('run_tclean_%s.py' %(names[x].split('.')[0]), editmode) as f:
            if not clean_has_run:
                paths_to_dirs.append(os.getcwd())
                paths_to_files.append(os.path.realpath(f.name))

            f.write("""widebandpbcor(vis='%s', imagename='%s',""" %(vises[x][0], names[x]) + """ nterms=2,
                   threshold='%sJy' %(threshold), action='pbcor', field='0', spwlist=[0,7,15], chanlist=[0,0,0],
                   weightlist=[1,1,1]) \n \n""")

        os.chdir('..')
    pb_has_run = True


def cutout():

    global pb_has_run, cutout_has_run

    if pb_has_run:
        edit_mode='a'
    else: edit_mode='w'
    for x in range(len(names)):
        os.chdir(names[x])

        with open('run_tclean_%s.py' %(names[x].split('.')[0]), edit_mode) as f:
            if not pb_has_run:
                paths_to_dirs.append(os.getcwd())
                paths_to_files.append(os.path.realpath(f.name))

            # make cutout image
            f.write("""imsubimage(imagename='%s.pbcor.image.tt0', outfile='%s.cutout.pbcor', overwrite=True,
                    region='box[[%spix, %spix], [%spix, %spix]]')\n \n""" % (names[x], names[x], cut_frame[0],
                                                                             cut_frame[0], cut_frame[1], cut_frame[1]))
            f.write("""imsubimage(imagename='%s.image.tt0', outfile='%s.cutout', overwrite=True,
                    region='box[[%spix, %spix], [%spix, %spix]]')\n \n""" % (names[x], names[x], cut_frame[0],
                                                                             cut_frame[0], cut_frame[1], cut_frame[1]))

            f.write("""exportfits(imagename='%s.cutout.pbcor', fitsimage='%s.cutout.pbcor.fits', overwrite=True)\n"""
                    % (names[x], names[x]))
            f.write("""exportfits(imagename='%s.cutout', fitsimage='%s.cutout.fits', overwrite=True)\n"""
                    % (names[x], names[x]))
            f.write("""exportfits(imagename='%s.pbcor.image.tt0', fitsimage='%s.pbcor.fits', overwrite=True)\n"""
                    % (names[x], names[x]))
            f.write("""exportfits(imagename='%s.image.tt0', fitsimage='%s.fits', overwrite=True)\n \n"""
                    % (names[x], names[x]))

        os.chdir('..')
    cutout_has_run = True


def statistics():
    global cutout_has_run

    if cutout_has_run:
        edit_mode = 'a'
    else:
        edit_mode = 'w'
    for x in range(len(names)):
        os.chdir(names[x])

        with open('run_tclean_%s.py' % (names[x].split('.')[0]), edit_mode) as f:
            if not cutout_has_run:
                paths_to_dirs.append(os.getcwd())
                paths_to_files.append(os.path.realpath(f.name))

            f.write("""stats=imstat('%s.residual.tt0')\n""" % (names[x]))
            f.write("""stdev=1.4826*stats['medabsdevmed'][0]\n""")
            f.write("""with open('stdev.txt', 'w') as f:\n""")
            f.write("""\tf.write('%s' %(stdev))\n \n""")
            f.write("""majoraxis = imhead(imagename='%s.image.tt0', mode='get', hdkey='bmaj')['value']\n"""
                    % (names[x]))
            f.write("""minoraxis = imhead(imagename='%s.image.tt0', mode='get', hdkey='bmin')['value']\n \n"""
                    % (names[x]))
            f.write("""beamarea = np.pi*majoraxis*minoraxis/(4*np.log(2))\n""")
            f.write("""with open('beamarea.txt', 'w') as f:\n""")
            f.write("""\tf.write('%s' %(beamarea))\n \n""")

        os.chdir('..')


# can make dirty images, then later clean
# or can run both one after another
clean_images()
pb_cor()
cutout()
statistics()

# generates the pipeline script
os.chdir('/users/gpetter/DATA')
with open('pipelinerun', 'w') as f:
    for x in range(len(paths_to_files)):
        f.write("""cd %s; xvfb-run -d casa -r 5.3.0-143 --nogui -c %s\n""" % (paths_to_dirs[x], paths_to_files[x]))




