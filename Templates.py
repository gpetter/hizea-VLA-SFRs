#
# Author: Grayson Petter

import numpy as np
from astropy import units as u
from astropy import constants as const
import math
import pickle
from scipy.io import readsav
from scipy.optimize import curve_fit
import WISE
from astropy.cosmology import FlatLambdaCDM
import pandas as pd
import glob

# speed of light
c = 299792458. # m/s

# path to project
projpath = '/Users/graysonpetter/Desktop/Dartmouth/HIZEA/hizea-VLA-SFRs/'

# set cosmology for calculating luminosity distance
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

# Read in all templates except for AGN (SFGs and Composite)
templates = glob.glob(projpath + 'Comprehensive_library/SFG*.txt')
templates.extend(
    glob.glob(projpath + 'Comprehensive_library/Comp*.txt'))

char_e = readsav('/Users/graysonpetter/Downloads/chary_elbaz_codes/chary_elbaz.save')
ce_temps = np.transpose(char_e['nuLnuinLsun'])
ce_waves = char_e['lambda']
ce_tots = char_e['Lir']
wise_bandpasses_3_4 = sorted(glob.glob(projpath + 'bandpass/*.txt'))[2:4]


# take in a template and redshift every wavelength
def redshift_spectrum(z, template, trim, table):
    if table:
        t = pd.read_csv(template, delim_whitespace=True, engine='python', header=None)
        #t.columns = ['lambda_e', 'L', 'dL']

        wavelengths = np.array(t.iloc[:, 0])
        Lums = np.array(t.iloc[:, 1])
    else:
        wavelengths = ce_waves
        Lums = template


    if trim:
        # cut template down to 8-1000 microns
        spec_range = np.where((wavelengths >= 8.) & (wavelengths <= 1000.))[0]
        wavelengths = wavelengths[spec_range]
        Lums = Lums[spec_range]

    shifted_len = np.array(wavelengths)*(1+z)

    # get luminosity at 12 & 22 micron
    twelve_mu = (np.abs(shifted_len - 12)).argmin()
    twenty_two_mu = (np.abs(shifted_len - 22)).argmin()

    return wavelengths, Lums, Lums[twelve_mu], Lums[twenty_two_mu], shifted_len


def interpolate_spec(shifted_spec, model):
    # convert wavelengths in microns to frequencies in Hz
    nus = (10 ** 6) * c / (shifted_spec[0])

    # reverse lists so frequencies go from low to high for simplicity
    reversed_nus = np.flipud(nus).flatten()
    reversed_lums = np.flipud(shifted_spec[1])

    # calculate constant interval to interpolate on
    if model:
        step = reversed_nus[1] - reversed_nus[0]
        dx = round(step, -(len(str(int(step))) - 1))
    else:
        dx = 10000000000

    # find smallest factor of dx Hz greater than the smallest frequency in the list
    start = (reversed_nus[0] + int(dx)) - (reversed_nus[0] % int(dx))

    # number of dx Hz intervals in entire template
    span = reversed_nus[len(reversed_nus) - 1] - reversed_nus[0]
    chunks = int(math.floor(span / dx))

    new_nus, new_lums = [], []
    current_nu = start

    # linearly interpolate to frequencies in dx Hz steps
    for x in range(chunks):
        new_nus.append(current_nu)
        new_lums.append(np.interp(current_nu, reversed_nus, reversed_lums))
        current_nu += dx

    return new_nus, new_lums, dx


# integrate spectrum using trapezoid method (rectangle below
def integrate_spectrum(freqs, Ls, dx):
    tot_ir = 0.
    for x in range(len(freqs) - 1):
        if Ls[x + 1] > Ls[x]:
            rect = dx * Ls[x]
            tri = (dx * (Ls[x + 1] - Ls[x])) / 2.

        else:
            rect = dx * Ls[x + 1]
            tri = (dx * (Ls[x] - Ls[x + 1])) / 2.

        tot_ir += (rect + tri)

    return tot_ir


def test_templates(zz):
    ratio_list = []
    with open(projpath + 'integrations/kirk.txt', 'rb') as fb:
        total_ir = np.array(pickle.load(fb))

    for x in range(len(templates)):
        shifted_spectrum = redshift_spectrum(zz, templates[x], True, True)
        twenty_two_lum = shifted_spectrum[3]
        ratio_list.append(total_ir[x] / twenty_two_lum)

    # ratio_list = np.log10(np.array(ratio_list))
    ratio_list = (np.array(ratio_list))
    averg = np.mean(ratio_list)
    stdev = np.std(ratio_list)
    # range_ratios = np.ptp(ratio_list)

    return stdev / averg


def simulate_wise_fluxes_for_colors(z, tems, bands, csv):
    tot_mag_list, template_names = [], []
    # iterate through templates
    for tem in tems:
        # redshift template
        red_spec = redshift_spectrum(z, tem, False, True)
        red_waves = np.array(red_spec[4])
        lumi = np.array(red_spec[1])

        normalized = []

        # iterate through WISE bands
        for y in range(len(bands)):
            if csv:
                band = pd.read_csv(bands[y], header=None, engine='python')
            else:
                band = pd.read_csv(bands[y], header=None, delim_whitespace=True, engine='python')
            bandwaves = np.array(band.iloc[:, 0])
            band_response = np.array(band.iloc[:, 1])


            # trim template to same wavelength range as WISE band
            cut = np.where((red_waves >= np.min(bandwaves)) & (red_waves <= np.max(bandwaves)))[0]
            trimmed_y = red_waves[cut]
            trimmed_L = lumi[cut]

            # interpolate template to band wavelengths, multiply by the response at that wavelength
            inter_lum = []
            for j in range(len(bandwaves)):
                inter_lum.append(band_response[j] * (np.interp(bandwaves[j], trimmed_y, trimmed_L)))

            # crude method
            """sum_lum = np.sum(np.array(inter_lum))
            sum_waves = np.sum(np.array(band_response))
            normalized.append(sum_lum/sum_waves)"""

            # integrate template multiplied by response function
            spectrum = [bandwaves, inter_lum]
            interped_again = interpolate_spec(spectrum, True)
            wise_lums = integrate_spectrum(interped_again[0], interped_again[1], interped_again[2])

            # integrate wise band
            band_spectrum = [bandwaves, band_response]
            interped_band = interpolate_spec(band_spectrum, True)
            integrated_band = integrate_spectrum(interped_band[0], interped_band[1], interped_band[2])
            # divide two
            normalized.append(wise_lums / integrated_band)

        tot_mag_list.append(normalized)

        template_names.append(tem.split('.txt')[0].split('/')[8])

    return tot_mag_list, template_names

def simulate_wise_fluxes(z, tem, bands, csv):
    tot_mag_list = []

    # redshift template
    red_spec = redshift_spectrum(z, tem, False, True)
    red_waves = np.array(red_spec[4])
    lumi = np.array(red_spec[1])

    normalized = []

    # iterate through WISE bands
    for y in range(len(bands)):
        if csv:
            band = pd.read_csv(bands[y], header=None, engine='python')
        else:
            band = pd.read_csv(bands[y], header=None, delim_whitespace=True, engine='python')
        bandwaves = np.array(band.iloc[:, 0])
        band_response = np.array(band.iloc[:, 1])
        band_apple = np.multiply(bandwaves, band_response)

        # trim template to same wavelength range as WISE band
        cut = np.where((red_waves >= np.min(bandwaves)) & (red_waves <= np.max(bandwaves)))[0]
        trimmed_y = red_waves[cut]
        trimmed_L = lumi[cut]

        # interpolate template to band wavelengths, multiply by the response at that wavelength
        inter_lum = []
        for i in range(len(bandwaves)):
            inter_lum.append(band_apple[i]  * (np.interp(bandwaves[i], trimmed_y, trimmed_L)))

        # crude method
        """sum_lum = np.sum(np.array(inter_lum))
        sum_waves = np.sum(np.array(band_response))
        normalized.append(sum_lum/sum_waves)"""

        # integrate template multiplied by response function
        spectrum = [bandwaves, inter_lum]
        interped_again = interpolate_spec(spectrum, True)
        wise_lums = integrate_spectrum(interped_again[0], interped_again[1], interped_again[2])

        # integrate wise band
        band_spectrum = [bandwaves, band_apple]
        interped_band = interpolate_spec(band_spectrum, True)
        integrated_band = integrate_spectrum(interped_band[0], interped_band[1], interped_band[2])
        # divide two
        normalized.append(wise_lums / integrated_band)
    return normalized



def chary_elbaz(zz):
    ratios = []
    for i in range(len(ce_tots)):
        shifted = redshift_spectrum(zz, ce_temps[i], True, False)
        twenty_two_lum = shifted[3]
        ratios.append(ce_tots[i]/twenty_two_lum)

    ratio_list = (np.array(ratios))
    averg = np.mean(ratio_list)
    stdev = np.std(ratio_list)
    # range_ratios = np.ptp(ratio_list)
    return stdev / averg

def writetotals():
    totlist = []
    for x in range(len(templates)):
        shifted_spectrum = redshift_spectrum(0, templates[x], True, True)
        interped_spectrum = interpolate_spec(shifted_spectrum, False)
        total_ir = integrate_spectrum(interped_spectrum[0], interped_spectrum[1], interped_spectrum[2])
        totlist.append(total_ir)
    with open(projpath + 'integrations/kirk.txt', 'wb') as fb:
        pickle.dump(totlist, fb)
writetotals()


def Kennicut1998(L_IR, L_ir_err):
    L_IR = L_IR.to('erg/s').value
    L_ir_err = L_ir_err.to('erg/s').value
    SFR = 4.5e-44*L_IR
    SFR_err = L_ir_err/L_IR*SFR
    return SFR, SFR_err

def murphyIRSFR(L_IR, L_ir_err):
    L_IR = L_IR.to('erg/s').value
    L_ir_err = L_ir_err.to('erg/s').value
    SFR = 3.88e-44 * L_IR
    SFR_err = L_ir_err / L_IR * SFR
    return SFR, SFR_err

def test_SFRs(z, name, table, tems=templates):

    d = cosmo.luminosity_distance(z)
    fluxes = WISE.mag_to_flux(name)
    w3_flux = fluxes[0] * u.Jy
    w3_flux_err = fluxes[2] * u.Jy
    w_four_good = False
    simulate_flux = True
    w3_lum = (w3_flux*4*np.pi*d**2).to('W/Hz')
    w3_lum_err = ((4*np.pi*d**2)*w3_flux_err).to('W/Hz')
    band_nine_fluxes = []

    # if there's data for W4
    if not np.isnan(fluxes[1]):
        w4_flux = fluxes[1] * u.Jy
        w4_flux_err = fluxes[3] * u.Jy
        w_four_good = True
        w4_lum = (w4_flux*4*np.pi*d**2).to('W/Hz')
        w4_lum_err = ((4*np.pi*d**2)*w4_flux_err).to('W/Hz')
    # which templates to use (kirk = kirkpatrick 2015, chary = chary & elbaz, both= both of them)
    key = 'kirk'

    SFRs, SFR_errs = [], []
    if key=='both':
        with open(projpath + 'integrations/kirk.txt', 'rb') as fb:
            total_ir = np.array(pickle.load(fb))
        for i, tem in enumerate(tems):

            tem_lum = redshift_spectrum(z, tem, True, table)

            l_ratio = w3_lum.value/tem_lum[2]
            #shifted_lums = np.array(tem_lum[1])*l_ratio

            #shifted_tem = [tem_lum[0], shifted_lums]
            """if i==0:
                plt.figure(0)
                plt.plot(tem_lum_w3[0], tem_lum_w3[1])
                plt.scatter(12, lum.value, c='k')
                plt.plot(tem_lum_w3[0], shifted_lums)
                plt.xlim(8, 15)
                plt.show()
                plt.close()
                plt.clf()
                plt.cla()"""

            #interped_spectrum = interpolate_spec(shifted_tem, True)
            #total_ir = integrate_spectrum(interped_spectrum[0], interped_spectrum[1], interped_spectrum[2])*u.W

            SFR = Kennicut1998(total_ir[i]*l_ratio*u.W)
            SFRs.append(SFR)

        for i, tem in enumerate(ce_temps):
            shifted = redshift_spectrum(z, tem, True, False)
            freq = (const.c/((12*u.micron))).to('Hz')
            twelve_micron = shifted[2]*u.solLum/freq
            l_ratio = w3_lum/(twelve_micron.to('W/Hz'))

            SFR = Kennicut1998((ce_tots[i]*u.solLum * l_ratio).to('W'))
            SFRs.append(SFR)

    elif key == 'chary':
        for i, tem in enumerate(ce_temps):
            shifted = redshift_spectrum(z, tem, True, False)
            freq = (const.c/((12*u.micron))).to('Hz')
            twelve_micron = shifted[2]*u.solLum/freq
            l_ratio = w3_lum/(twelve_micron.to('W/Hz'))

            SFR = Kennicut1998((ce_tots[i]*u.solLum * l_ratio).to('W'))
            SFRs.append(SFR)

    elif key == 'kirk':
        with open(projpath + 'integrations/kirk.txt', 'rb') as fb:
            total_ir = np.array(pickle.load(fb))
        for i, tem in enumerate(tems):
            tem_lum = redshift_spectrum(z, tem, False, table)

            if w_four_good:
                wavelengthsaye = tem_lum[4]
                maxcutoff = min((np.where(wavelengthsaye > 25.0)[0]))
                mincutoff = max(np.where(wavelengthsaye < 7.)[0])
                lums = tem_lum[1][mincutoff:maxcutoff]
                lams = np.array(tem_lum[4])[mincutoff:maxcutoff]

                fit = np.flipud(np.polyfit(lams, lums, 12))

                lambdas = np.array([float(12.082), float(22.194)])
                measured_lums = np.array([float(w3_lum.value), float(w4_lum.value)])
                measured_lum_errs = np.array([float(w3_lum_err.value), float(w4_lum_err.value)])


                def func(x, a):
                    return a * sum((q * x ** j for j, q in enumerate(fit)))
                if simulate_flux:
                    simulated = np.array(simulate_wise_fluxes(z, tem, wise_bandpasses_3_4, False))

                    #def func(a):
                    #    return a*simulated

                    l_ratio = (measured_lums[0]*simulated[0]/(measured_lum_errs[0])**2 + measured_lums[1]*simulated[1]/(measured_lum_errs[1])**2)/((simulated[0]/measured_lum_errs[0])**2 + (simulated[1]/measured_lum_errs[1])**2)
                    normalization_percent_err = 0
                else:
                    popt, pcov = curve_fit(func, lambdas, measured_lums, sigma=measured_lum_errs)
                    l_ratio = float(popt)
                    normalization_percent_err = np.sqrt(float(pcov)) / float(popt)

            else:
                l_ratio = float(w3_lum.value/tem_lum[2])
                """shifted_lums = np.array(tem_lum[1])*l_ratio
                shifted_tem = [tem_lum[0], shifted_lums]
                interped_spectrum = interpolate_spec(shifted_tem, True)
                total_ir = integrate_spectrum(interped_spectrum[0], interped_spectrum[1], interped_spectrum[2])*u.W"""
                normalization_percent_err = float(w3_lum_err/w3_lum)
            L_ir_tot = total_ir[i]*l_ratio*u.W
            L_ir_tot_err = normalization_percent_err*L_ir_tot
            #print(normalization_percent_err)
            SFR = murphyIRSFR(L_ir_tot, L_ir_tot_err)
            SFRs.append(SFR[0])
            """idxboob = np.abs((np.array(tem_lum[4])-450.)).argmin()
            band_nine_lum = ((tem_lum[1])[idxboob]*l_ratio)*u.J
            band_nine_fluxes.append(((band_nine_lum/(4*np.pi*d**2))[0].to('Jy')).value)

            SFR_errs.append(SFR[1])
        print(name, min(band_nine_fluxes), max(band_nine_fluxes))"""

    SFR_uncertainty = np.sqrt(np.sum(np.square(SFR_errs)))/len(SFR_errs)


    return np.average(SFRs), np.std(SFRs), SFR_uncertainty


