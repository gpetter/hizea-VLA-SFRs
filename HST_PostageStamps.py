import aplpy
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import GetGalaxyList
reload(GetGalaxyList)


current_dir = os.getcwd()

names = GetGalaxyList.return_galaxy_list()

# Sort by RA
stripped_names = []
sorted_names = []

for name in names:
    new_name = name.split('J')[1]
    stripped_names.append(new_name)

sorted_nums = sorted(stripped_names)

for thing in sorted_nums:
    sorted_names.append('J' + thing)

num_gals = len(sorted_names)
if num_gals/5. >= 4.:
    cols = 5
else:
    cols = 4
rows = 5

# Create figure
fig = plt.figure(10, figsize=(cols*6, rows*6), dpi=200)

spacing = 0.0

x_size = (1. - spacing)/cols
y_size = (1. - spacing)/rows
x_start = spacing/(cols+1.)
y_start = 1. - y_size - spacing/(rows+1.)
x_iter = 0
y_iter = 0

for z in range(len(sorted_names)):
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

    HST_name = sorted_names[z][:5] + '_HST.fits'

    imgname = '%s.cutout.pbcor.fits' % sorted_names[z]

    with open('text/center_HST.txt', 'r') as f_center:
        lines = f_center.readlines()
        HST_ra = float(lines[0])
        HST_dec = float(lines[1])
    with open('text/stdev.txt', 'r') as f_rms:
	rms = float(f_rms.readline())

    f = aplpy.FITSFigure(HST_name, figure=fig, subplot=[x, y, x_size, y_size])

    f.show_grayscale()
    f.show_contour(imgname, levels=[(2.*rms), (3.*rms), (6*rms), (12*rms)], alpha=1, linewidths=4., colors='cyan')
    f.recenter(HST_ra, HST_dec, width=15./3600, height=15./3600)
 
    #f.add_scalebar(5./3600.)
    #f.scalebar.show(5./3600.)
    #f.scalebar.set_corner('bottom right')
    #f.scalebar.set_label('5 arcseconds')

    f.axis_labels.hide_y()
    f.axis_labels.hide_x()
    f.ticks.hide()
    f.tick_labels.hide()

    f.add_label(0.2, 0.9, ('%s' % sorted_names[z])[:5], relative=True, weight=weighting, color='orange', size=50)
    f.set_theme('publication')

    os.chdir('..')

os.chdir(current_dir)
plt.savefig('HST_Stamps.png', bbox_inches='tight', pad_inches=0)
plt.clf()
plt.close()
plt.clf()
plt.cla()
plt.close()
