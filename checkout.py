#!/usr/bin/env python

from figure_of_merit import MoonGoT as mgt

moon = mgt()
moon.set_dtg()
moon.phi()
moon.theta_moon()
moon.freq_input("2200MHz")
moon.y_input("2.2db")
moon.calc_wavelength()
moon.calc_k1()
moon.calc_k2()
moon.calc_moon_flux()

print moon
