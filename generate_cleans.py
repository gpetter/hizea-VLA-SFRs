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
    alldirs = os.listdir(name)
    onlyms = [y for y in alldirs if '.ms' in y]
    vises.append(onlyms)

pathstofiles, pathstodirs = [], []
dirtyhasrun = False
cleanhasrun = False
pbhasrun = False
cutouthasrun = False

# parameters
############################################################################
cleansigma = 3
imgsize = 12000
cutoutsize = 200
cutframe = [imgsize/2-cutoutsize/2, imgsize/2+cutoutsize/2]
############################################################################

# Creates a python script in each galaxy's directory
# This script will construct a dirty image, or a very lightly cleaned image (niter~=100)
# Measures the MAD of that image, scales it by 1.4826, and multiplies it by some number to set the threshold
# In the process, a PSF and residual is saved so it does not have to be recalculated when the images are cleaned

# Generate a script in each directory to construct a dirty image, then measure that image's variance
def makedirtyimages():
    global dirtyhasrun
    for x in range(len(names)):
        os.chdir(names[x])
        with open('run_tclean_%s.py' %(names[x].split('.')[0]), 'w') as f:

            # adding paths to simplify creating pipeline script later
            pathstodirs.append(os.getcwd())
            pathstofiles.append(os.path.realpath(f.name))

            # write out dirty tclean command
            f.write("""tclean(vis=%s, imagename='%s', field='0', datacolumn='data',
                verbose=True, gridder='wproject', wprojplanes=128, pblimit=-1, robust=0.5, imsize=[%s], cell='0.2arcsec',
                specmode='mfs', deconvolver='mtmfs', nterms=2, scales=[0,11,28], interactive=False, niter=50,
                weighting='briggs', usemask='auto-multithresh', stokes='I', threshold='0.0Jy', calcpsf=True,
                calcres=True, savemodel='modelcolumn', restart=False) \n \n""" % (vises[x], names[x], imgsize))

            # will call imstat to measure the MAD of each image, scaled by number*1.4826*MAD
            # saves threshold value to text file for CleanImages() script's access
            f.write("""stats=imstat('%s.image.tt0')\n""" %(names[x]))
            f.write("""thresh=%s*1.4826*stats['medabsdevmed'][0]\n""" % (cleansigma))
            f.write("""print(thresh)\n""")
            f.write("""with open('threshold.txt', 'w') as f:\n""")
            f.write("""\tf.write('%s' %(thresh))\n \n""")

        os.chdir('..')
    dirtyhasrun=True


def cleanimages():

    global dirtyhasrun, cleanhasrun
    if dirtyhasrun:
        editmode='a'
    else: editmode='w'
    for x in range(len(names)):
        os.chdir(names[x])

        with open('run_tclean_%s.py' %(names[x].split('.')[0]), editmode) as f:
            if not dirtyhasrun:
                pathstodirs.append(os.getcwd())
                pathstofiles.append(os.path.realpath(f.name))

            # retrieve previously saved threshold value
            f.write("""with open('threshold.txt', 'r') as f:\n""")
            f.write("""\tglobal threshold\n""")
            f.write("""\tlines=f.readlines()\n""")
            f.write("""\tthreshold=float(lines[0])\n""")
            f.write("""print(threshold)\n \n""")

            # run tclean with said threshold, high niter value
            f.write(("""tclean(vis=%s, imagename='%s', field='0', datacolumn='data',
                   verbose=True, gridder='wproject', wprojplanes=128, pblimit=-1, robust=0.5, imsize=[%s], cell='0.2arcsec', 
                   specmode='mfs', deconvolver='mtmfs', nterms=2, scales=[0,11,28], """ % (vises[x], names[x], imgsize))
                    + """ interactive=False, niter=20000, weighting='briggs',
                   usemask='auto-multithresh', sidelobethreshold = 4.0, stokes='I', threshold='%sJy' %(threshold), savemodel='modelcolumn', calcres=False, calcpsf=False, 
                   restart=True) \n \n""")


        os.chdir('..')
    cleanhasrun = True


def pbcor():

    global cleanhasrun, pbhasrun
    if cleanhasrun:
        editmode='a'
    else: editmode='w'
    for x in range(len(names)):
        os.chdir(names[x])
        with open('run_tclean_%s.py' %(names[x].split('.')[0]), editmode) as f:
            if not cleanhasrun:
                pathstodirs.append(os.getcwd())
                pathstofiles.append(os.path.realpath(f.name))

            f.write("""widebandpbcor(vis='%s', imagename='%s',""" %(vises[x][0], names[x]) + """ nterms=2,
                   threshold='%sJy' %(threshold), action='pbcor', field='0', spwlist=[0,7,15], chanlist=[0,0,0],
                   weightlist=[1,1,1]) \n \n""")

        os.chdir('..')
    pbhasrun = True


def cutout():

    global pbhasrun, cutouthasrun

    if pbhasrun:
        editmode='a'
    else: editmode='w'
    for x in range(len(names)):
        os.chdir(names[x])

        with open('run_tclean_%s.py' %(names[x].split('.')[0]), editmode) as f:
            if not pbhasrun:
                pathstodirs.append(os.getcwd())
                pathstofiles.append(os.path.realpath(f.name))

            # make cutout image
            f.write("""imsubimage(imagename='%s.pbcor.image.tt0', outfile='%s.cutout.pbcor', overwrite=True,
                    region='box[[%spix, %spix], [%spix, %spix]]')\n \n""" % (names[x], names[x], cutframe[0],
                                                                             cutframe[0], cutframe[1], cutframe[1]))
            f.write("""imsubimage(imagename='%s.image.tt0', outfile='%s.cutout', overwrite=True,
                    region='box[[%spix, %spix], [%spix, %spix]]')\n \n""" % (names[x], names[x], cutframe[0],
                                                                             cutframe[0], cutframe[1], cutframe[1]))

            f.write("""exportfits(imagename='%s.cutout.pbcor', fitsimage='%s.cutout.pbcor.fits', overwrite=True)\n"""
                    % (names[x], names[x]))
            f.write("""exportfits(imagename='%s.cutout', fitsimage='%s.cutout.fits', overwrite=True)\n"""
                    % (names[x], names[x]))
            f.write("""exportfits(imagename='%s.pbcor.image.tt0', fitsimage='%s.pbcor.fits', overwrite=True)\n"""
                    % (names[x], names[x]))
            f.write("""exportfits(imagename='%s.image.tt0', fitsimage='%s.fits', overwrite=True)\n \n"""
                    % (names[x], names[x]))

        os.chdir('..')
    cutouthasrun = True


def statistics():
    global cutouthasrun

    editmode = ''

    if cutouthasrun:
        editmode = 'a'
    else:
        editmode = 'w'
    for x in range(len(names)):
        os.chdir(names[x])

        with open('run_tclean_%s.py' % (names[x].split('.')[0]), editmode) as f:
            if not cutouthasrun:
                pathstodirs.append(os.getcwd())
                pathstofiles.append(os.path.realpath(f.name))

            f.write("""stats=imstat('%s.residual.tt0')\n""" % (names[x]))
            f.write("""stdev=1.4826*stats['medabsdevmed'][0]\n""")
            f.write("""with open('stdev.txt', 'w') as f:\n""")
            f.write("""\tf.write('%s' %(stdev))\n \n""")
            f.write("""majoraxis = imhead(imagename='%s.image.tt0', mode='get', hdkey='bmaj')['value']\n""" % (names[x]))
            f.write("""minoraxis = imhead(imagename='%s.image.tt0', mode='get', hdkey='bmin')['value']\n \n""" % (names[x]))
            f.write("""beamarea = np.pi*majoraxis*minoraxis/(4*np.log(2))\n""")
            f.write("""with open('beamarea.txt', 'w') as f:\n""")
            f.write("""\tf.write('%s' %(beamarea))\n \n""")

        os.chdir('..')
# can make dirty images, then later clean
# or can run both one after another
makedirtyimages()

# generates the pipeline script
os.chdir('/users/gpetter/DATA')
with open('pipelinerun', 'w') as f:
    for x in range(len(pathstofiles)):
        f.write("""cd %s; xvfb-run -d casa -r 5.3.0-143 --nogui -c %s\n""" % (pathstodirs[x], pathstofiles[x]))




