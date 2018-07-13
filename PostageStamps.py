import aplpy
import matplotlib.pyplot as plt
import os
from astropy.io import fits

# Go to directory, get list of galaxies
os.chdir('/users/gpetter/DATA')
names = os.listdir('data_v1')
os.chdir('data_v1')

fig = plt.figure(10, figsize=(23, 30), dpi=100)

x_size = .16
y_size = .12
x_start = .072
y_start = .8133333
x_iter = 0
y_iter = 0

for z in range(len(names)):
    os.chdir(names[z])

    if z % 4 == 0 and z != 0:
        x_iter = 0
        y_iter = y_iter - (1-y_start)

    x = x_start + x_iter
    y = y_start + y_iter

    x_iter = x_iter + x_size + x_start


    imgname = '%s.cutout.pbcor.fits' % names[z]
    hdu = fits.open(imgname)
    bmaj = hdu[0].header['bmaj']
    bmin = hdu[0].header['bmin']
    angle = hdu[0].header['bpa']



    f = aplpy.FITSFigure(imgname, figure=fig, subplot=[x, y, x_size, y_size])
    center = f.pixel2world(100, 100)
    f.show_circles([center[0]], [center[1]], [4.0 / 3600.0], edgecolor='magenta')
    f.show_colorscale()
    f.set_title('%s' % names[z])
    f.add_beam()
    f.beam.set_major(bmaj)
    f.beam.set_minor(bmin)
    f.beam.set_angle(angle)
    f.beam.show()
    f.beam.set_corner('bottom left')
    f.beam.set_color('white')










    os.chdir('..')

os.chdir('/users/gpetter/PycharmProjects/untitled')
plt.savefig('Stamps.png')
plt.clf()
plt.close()