import aplpy
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import GetGalaxyList
reload(GetGalaxyList)

current_dir = os.getcwd()

names = GetGalaxyList.return_galaxy_list(1)

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

    HST_name = sorted_names[z][:5] + '_HST.fits'

    imgname = '%s.cutout.pbcor.fits' % sorted_names[z]

    f = aplpy.FITSFigure(HST_name, figure=fig, subplot=[x, y, x_size, y_size])

    f.show_grayscale()
    f.show_contour(imgname, alpha=.5)

    f.set_title('%s' % sorted_names[z], weight=weighting)

    os.chdir('..')

os.chdir(current_dir)
plt.savefig('HST_Stamps.png')
plt.clf()
plt.close()