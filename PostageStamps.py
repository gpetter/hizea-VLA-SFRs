import aplpy
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import GetGalaxyList
reload(GetGalaxyList)

current_dir = os.getcwd()

sorted_names = GetGalaxyList.return_galaxy_list()



num_gals = len(sorted_names)
if num_gals/5. >= 4.:
    cols = 5
else:
    cols = 4
rows = 5

# Create figure
fig = plt.figure(10, figsize=(cols*6, rows*6), dpi=200)

spacing = 0.35

x_size = (1. - spacing)/cols
y_size = (1. - spacing)/rows
x_start = spacing/(cols+1.)
y_start = 1. - y_size - spacing/(rows+1.)
x_iter = 0
y_iter = 0


for z in range(num_gals):
    os.chdir(sorted_names[z])

    if z % cols == 0 and z != 0:
        x_iter = 0
        y_iter = y_iter - (1-y_start)

    x = x_start + x_iter
    y = y_start + y_iter

    x_iter = x_iter + x_size + x_start

    with open('text/detect.txt', 'r') as f:
        truth = f.readline()

    if int(truth) == 1:
        weighting = 'bold'
    else:
        weighting = 'normal'

    imgname = '%s.cutout.pbcor.fits' % sorted_names[z]
    hdu = fits.open(imgname)
    bmaj = hdu[0].header['bmaj']
    bmin = hdu[0].header['bmin']
    angle = hdu[0].header['bpa']

    f = aplpy.FITSFigure(imgname, figure=fig, subplot=[x, y, x_size, y_size])

    #center = f.pixel2world(100, 100)
    #f.show_circles([center[0]], [center[1]], [4.0 / 3600.0], edgecolor='magenta')
    with open('text/center_HST.txt', 'r') as f_center:
        lines = f_center.readlines()
        HST_ra = float(lines[0])
        HST_dec = float(lines[1])
    with open('text/width.txt', 'r') as f_width:
        lines = f_width.readlines()
        maj_width = float(lines[0])
        min_width = float(lines[1])
        PA = float(lines[2])

    f.show_ellipses(HST_ra, HST_dec, min_width/3600., maj_width/3600., angle=PA, edgecolor='magenta')
    #f.show_markers(HST_ra, HST_dec, marker='x', facecolor='c')

    f.show_colorscale()

    f.set_title('%s' % sorted_names[z], weight=weighting)

    f.add_beam()
    f.beam.set_major(bmaj)
    f.beam.set_minor(bmin)
    f.beam.set_angle(angle)
    f.beam.show()
    f.beam.set_corner('bottom left')
    f.beam.set_color('white')

    f.add_colorbar()
    f.colorbar.show(log_format=False)
    f.colorbar.set_width(0.1)
    f.colorbar.set_pad(0.03)
    f.colorbar.set_axis_label_text('Surface Brightness (Jy/beam)')
    f.colorbar.set_axis_label_pad(9)
    f.colorbar.set_axis_label_rotation(270)
    f.colorbar.set_font(size='xx-small', weight='medium', stretch='extra-condensed', family='sans-serif', style='normal',
                        variant='normal')
    f.colorbar.set_axis_label_font(size='x-small')

    #f.axis_labels.hide_y()
    f.axis_labels.set_ypad(0)
    f.tick_labels.set_font(size='x-small', weight='medium', stretch='normal', family='sans-serif', style='normal',
                        variant='normal')

    os.chdir('..')

os.chdir(current_dir)
plt.savefig('Stamps.png')
plt.clf()
plt.close()