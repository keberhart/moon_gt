#!/usr/bin/env python
#
#   lets calculate the atmospheric attenuation from each of our
#   antennas.
#
#   1/2019 - Kyle Eberhart - of cours this wouldn't work without the itur lib
#
#----------------------------------------------------------------

import sys

arguments = len(sys.argv) - 1
if arguments is not 5:
    print "\n"
    print "This command requires Temperature(C), Pressure(hPa), Humidity(%), Frequency(GHz) and Elevation angle(Degrees)."
    print "For example 'atmo_atten.py T P H F E'"
    sys.exit()

temp = float(sys.argv[1])
press = float(sys.argv[2])
humid = float(sys.argv[3])
freq = float(sys.argv[4])
elev = float(sys.argv[5])

import itur

def calc_values(lat, lon, diam, alt):
    '''
    Calculate and print the atten for a given site using the command
        line arguments.
    '''
    T = temp * itur.u.deg_C
    P = press * itur.u.hPa
    H = humid
    f = freq * itur.u.GHz
    el = elev

    hs = alt * itur.u.km
    D = diam * itur.u.m

    p = 0.1

    # calculated atmospheric parameters
    rho_p = itur.surface_water_vapour_density(lat, lon, p, hs)
#    rho_sa = itur.models.itu835.water_vapour_density(lat, hs)
#    T_sa = itur.models.itu835.temperature(lat, hs)
#    V = itur.models.itu836.total_water_vapour_content(lat, lon, p, hs)

    # rain and cloud prediction parameters
#    R_prob = itur.models.itu618.rain_attenuation_probability(lat, lon, el, hs)
#    R_pct_prob = itur.models.itu837.rainfall_probability(lat, lon)
#    R001 = itur.models.itu837.rainfall_rate(lat, lon, p)
#    h_0 = itur.models.itu839.isoterm_0(lat, lon)
#    h_rain = itur.models.itu839.rain_height(lat, lon)
#    L_red = itur.models.itu840.columnar_content_reduced_liquid(lat, lon, p)
#    A_w = itur.models.itu676.zenit_water_vapour_attenuation(lat, lon, p , f, h=hs)

    # compute attenuation values
    A_g = itur.gaseous_attenuation_slant_path(f, el, rho_p, P, T)
    A_r = itur.rain_attenuation(lat, lon, f, el, hs=hs, p=p)
    A_c = itur.cloud_attenuation(lat, lon, el, f, p)
    A_s = itur.scintillation_attenuation(lat, lon, f, el, p, D)
    A_t = itur.atmospheric_attenuation_slant_path(lat, lon, f, el, p, D)

    print("\n")
    print("- Rain attenuation          ", A_r, " dB")
    print("- Gaseous attenuation       ", A_g, " dB")
    print("- Cloud attenuation         ", A_c, " dB")
    print("- Scintillation attenuation ", A_s, " dB")
    print("- Total attenuation         ", A_t, " dB")


# 13M-1
calc_values(64.974, -147.50989, 13, .367)
