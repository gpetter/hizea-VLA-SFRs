# Author: Grayson Petter
import numpy as np


def sig_figs(lum_flag, figs, value, uncertainty):
    output = []

    if lum_flag:
        sig_fig_uncertainty = ('%.' + '%sg' % figs) % float(uncertainty)
        thing2 = int(sig_fig_uncertainty.split('+')[1])

        thing3 = '%s' % value
        thing4 = int(thing3.split('+')[1])

        difference = thing4 - thing2
        idx = thing4 - difference - 1

        output.append(round(value, -idx))

    elif float(uncertainty) > 1.0:
        sig_fig_uncertainty = '%f' % float(('%.' + '%sg' % figs) % float(uncertainty))
        arrayed = map(int, sig_fig_uncertainty.split('.')[0])

        reversed_arr = np.array(arrayed[::-1])

        idx = min(np.where(reversed_arr > 0.0)[0])
        output.append(round(float(value), -idx))
    else:
        sig_fig_uncertainty = '%f' % float(('%.' + '%sg' % figs) % float(uncertainty))
        arrayed = np.array(map(float, sig_fig_uncertainty.split('.')[1]))

        idx = min(np.where(arrayed > 0.0)[0]) + figs
        output.append(round(float(value), idx))

    output.append(float(sig_fig_uncertainty))
    return output
