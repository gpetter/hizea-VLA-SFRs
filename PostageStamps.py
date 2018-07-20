import aplpy
import matplotlib.pyplot as plt
import os
from astropy.io import fits

# Go to directory, get list of galaxies
os.chdir('/users/gpetter/DATA')
names = os.listdir('data_v1')
os.chdir('data_v1')

# Sort by RA
stripped_names = []
sorted_names = []

for name in names:
    new_name = name.split('J')[1]
    stripped_names.append(new_name)

sorted_nums = sorted(stripped_names)

for thing in sorted_nums:
    sorted_names.append('J' + thing)

# Create figure
fig = plt.figure(10, figsize=(23, 30), dpi=200)

x_size = .16
y_size = .12
x_start = .072
y_start = .8133333
x_iter = 0
y_iter = 0

for z in range(len(sorted_names)):
    os.chdir(sorted_names[z])

    if z % 4 == 0 and z != 0:
        x_iter = 0
        y_iter = y_iter - (1-y_start)

    x = x_start + x_iter
    y = y_start + y_iter

    x_iter = x_iter + x_size + x_start

    with open('detect.txt', 'r') as f:
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
    center = f.pixel2world(100, 100)
    f.show_circles([center[0]], [center[1]], [4.0 / 3600.0], edgecolor='magenta')
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

os.chdir('/users/gpetter/PycharmProjects/untitled')
plt.savefig('Stamps.png')
plt.clf()
plt.close()