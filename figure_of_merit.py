#!/usr/bin/env python
#
#   Can we calculate the flux density of the moon and an antenna G/T based off
#   that?
#
#   1/2019 - Kyle Eberhart
#
#   Based on - Darko Sekuljica - "Using the Moon as a calibrated noise source
#   to measure the G/T figure-of-merit of an X-band satellite receiving station
#   with a large antenna 200...400 wavelengths in diameter"
#   http://antena.fe.uni-lj.si/literatura/Razno/Diplome/Sekuljica/Master%20Thesis%20-%20popravki%2017-01-2017.pdf
#
#   William C. Daywitt - "An Error Analysis for the Use of Presently Avaliable
#   Lunar Radio Flux Data in Broadbeam Antenna-System Measurements"
#
#   This project made use of Skyfield, http://rhodesmill.org/skyfield/
#
#   This project made use of Astropy, http://www.astropy.org a
#   community-developed core Python package for Astronomy astropy:2013,
#   astropy:2018
#
#------------------------------------------------------------------

import math
from skyfield.api import load
from skyfield.api import Topos
from skyfield.almanac import phase_angle
from astropy import units as u
from datetime import timedelta
import itur

class MoonGoT:
    '''Wrap up the Moon Flux and G/T in a class'''

    def __init__(self, test=None):
        # some constants and globals to use
        self.press = None
        self.temp = None
        self.humid = None
        self.k1 = None
        # Boltzmanns constant (m^2 kg s^-2 K^-1)
        self.k_bolt = 1.38064852e-23
        # Plancks constant (m^2 kg s^-1)
        self.h_planck = 6.62607004e-34
        # the speed of light in (m s^-1)
        self.c = 299792458
        self.moon_dia_km = 3476

        # Configure and pre-calculate a bit
        self.ts = load.timescale()

        self.e = load('de421.bsp')
        self.earth, self.luna = self.e['earth'], self.e['moon']

        if test is None:
            # The 13M antennas are close enough together to make no real difference.
            # We will use 13M-2 as the location of reference.

#-------------------------------------------------------------
#   Change this value to the correct date and time of measurement
            self.t = self.ts.utc(2019, 1, 27, 14, 10)
#
#------------------------------------------------------------
            #self.t = self.ts.now()
            self.diam = 13.0
            self.lat = 64.97381
            self.lon = -147.50575
            self.alt = 385.0
            self.observer = self.earth + Topos(latitude_degrees=self.lat,
                                                longitude_degrees=self.lon,
                                                elevation_m=self.alt)
            self.phi()
            self.theta_moon()
        else:
            # use the values from the paper to check stuff
            self.t = self.ts.utc(2016, 10, 13, 19, 23)
            self.diam = 11.28
            self.frequency = 8177500000
            self.y = 2.226
            self.k2est = 5.71
            self.lat = 40.6664
            self.lon = 16.6043
            self.alt = 100.0
            self.observer = self.earth + Topos(latitude_degrees=self.lat,
                                                longitude_degrees=self.lon,
                                                elevation_m=self.alt)

            self.phi()
            self.theta_moon()
            self.calc_g_over_t()


    def phi(self, time=None):
        '''Calculate the lunar phase angle in degrees'''
        # uses the skyfield library which seems to calculate the phase
        # angle of the moon weird. Its 180 out from the regular convention that
        # I can find.

        if time is None:
            time = self.t

        pa_ = phase_angle(self.e, "moon", time)
        pa = 180 - (pa_.to(u.deg)).value

        return pa

    def theta_moon(self, time=None):
        # Moons angular diameter in degrees
        if time is None:
            time = self.t

        astrometric = self.observer.at(time).observe(self.luna)
        el, az, d = astrometric.apparent().altaz()
        self.elevation = el.to(u.deg)
        self.azimuth = az.to(u.deg)
        self.distance = d.to(u.km)

        arcsec = ((206265 * self.moon_dia_km*u.km) / self.distance)*u.arcsec
        self.theta = (arcsec.to(u.deg)).value

        return self.theta

    def calc_t_moon(self, deg, f=None):
        # The average brightness temperature of the Moon in (K)
        # f is frequency in Hz
        # deg is phi, the lunar phase angle in degrees 
        #   - if phase angle is degreasing use 360 - phi
        #
        # from equation [4.14 - 4.16]
        #
        if f is None:
            f = self.frequency

        self.deg = deg

        five_min = timedelta(minutes=5)

        now_five = self.ts.utc(self.t.utc_datetime() + five_min)

        p1 = self.phi(now_five)

        self.t0 = 207.7 + (24.43/(f*10**-9))

        self.t1overt0 = 0.004212 * math.pow((f*10**-9),1.224)

        self.psi = 43.83 / (1 + 0.0109 * (f*10**-9))

        if deg > p1:
            self.t_moon = self.t0 * ( 1 - self.t1overt0 * math.cos(math.radians((360-deg) - self.psi)))
        #    print(360-deg, "phi -lunar phase angle in degrees, descending value calc.")
            return self.t_moon
        else:
            self.t_moon = self.t0 * ( 1 - self.t1overt0 * math.cos(math.radians(deg - self.psi)))
         #   print(deg, "phi- lunar phase angle in degrees, ascending value calc.")
            return self.t_moon

    def calc_moon_flux(self, f=None):
        # this is all from a paper by Darko Sekuljica, "Using the Moon as a calibrated
        # noise source to measure the G/T ...." it goes on. At any rate equation [4.13]
        #
        #   s_moon = (2*k_bolt*pi^3)/(4*180^2*c^2) * f^2 * T_moon * Omega_moon
        #
        # Flux density of the moon in (W m^-2 Hz-1)
        #
        if f is None:
            f = self.frequency

        # The first term? is all constants...
        term_1 = (2*self.k_bolt*math.pow(math.pi, 3))/(4*math.pow(180,2)*math.pow(self.c,2))

        term_2 = math.pow(f, 2)

        term_3 = self.calc_t_moon(self.phi())

        term_4 = math.pow(self.theta, 2)

        self.s_moon = term_1 * term_2 * term_3 * term_4

        return self.s_moon

    def calc_wavelength(self, f=None):
        # convert frequency to wavelength
        if f is None:
            f = self.frequency
        self.wavelength = float(self.c)/float(f)
        return self.wavelength

    def calc_hpbw(self):
        '''calculate the Half Power Beam Width of our antenna.'''
        # Our 13M antennas have the S-Band feed at the prime focus and the
        # X-band feed is double-shaped, with two reflectors.
        #
        # also I have learned that this is not necissarily the best formula for
        # calculaing HPBW, but we ned to know the edge taper of the dish for
        # better results.
        if self.frequency > 3*10**9:
            self.hpbw = 59.0 * self.wavelength / self.diam
        else:
            self.hpbw = 67.5 * self.wavelength / self.diam
        return self.hpbw

    def calc_k2(self, diam=None):
        # Calculate the correction factor for angular source size and antennas HPBW
        # This seems to work for X-Band frequencies but it is throwing an error
        # when the source is smaller than the beamwidth. Such as at S-Band frequencies.
        if diam is None:
            diam = self.diam

        moonbw = self.theta
        self.calc_hpbw()

        if self.hpbw/moonbw > 1:
            # When the source is narrower than the HPBW from the DOC study
#            x2 = 0.6441*math.pow(moonbw/self.hpbw, 2)
#            self.k2est = (1 - math.pow(math.e, -x2))/x2

            # from our sun worksheet not sure of its origin, I wouldn't expect
            # it to work in this case...
#            self.k2est = math.pow(1 + .18 * math.pow(self.theta/self.hpbw, 2), 2)

            # from our cas-a worksheet, from Datron we think
            x2 = (self.hpbw*60)
            self.k2est = (1-0.041025/x2
                    +8.00685/math.pow(x2, 2)
                    -10.673775/math.pow(x2, 3)
                    +41.662351/math.pow(x2, 4))
        else:
            # when the source is wider than the HPBW
            common_bit = math.log(2) * math.pow(moonbw/self.hpbw, 2)
            self.k2est = common_bit / 1 - (1 / math.pow(math.e, common_bit))

        return self.k2est

    def lin_to_db(self, val):
        '''turn a amp factor to a dB value.'''
        result = 10*math.log10(val)
        return result

    def db_to_lin(self, val):
        '''turn a dB to amp factor.'''
        result = math.pow(10, (val/10))
        return result

    def calc_g_over_t(self, y=None, f=None, k1=None, k2=None):
        '''calculate the G/T using this fun pile of a class.'''
        if y is None:
            y = self.y
        else:
            self.y = y
        if f is None:
            f = self.frequency
        else:
            self.frequency = f

        self.calc_wavelength()

        if k1 is None:
            if self.temp is None:
                k1 = 1.03
            else:
                k1 = self.calc_k1()
        else:
            self.k1 = k1
        if k2 is None:
            k2 = self.calc_k2()
        else:
            k2 = self.k2est

        S = self.calc_moon_flux()

        self.gt_lin = ((8*math.pi*self.k_bolt*(y-1))/(math.pow(self.wavelength,2)*S)*k1*k2)
        try:
            self.gt_db = 10*math.log10(self.gt_lin)
        except ValueError,e:
            self.gt_db = 0
            print "Somethings wrong with a log by less than 0."
            print str(e)

    def calc_k1(self):
        '''Calculate the k1 atmospheric attenuation value.'''
        # this uses the cool itur equations, but it is really slow...
        T = self.temp * itur.u.deg_C
        P = self.press * itur.u.hPa
        H = self.humid
        f = (self.frequency*10**-9) * itur.u.GHz
        el = self.elevation
        hs = (self.alt*10**-3) * itur.u.km
        D = self.diam * itur.u.m
        p = 0.1

        # calculated atmospheric parameters
        rho_p = itur.surface_water_vapour_density(self.lat, self.lon, p, hs)

        # compute attenuation values
#        A_g = itur.gaseous_attenuation_slant_path(f, el, rho_p, P, T)
#        A_r = itur.rain_attenuation(lat, lon, f, el, hs=hs, p=p)
#        A_c = itur.cloud_attenuation(lat, lon, el, f, p)
#        A_s = itur.scintillation_attenuation(lat, lon, f, el, p, D)
        A_t = itur.atmospheric_attenuation_slant_path(self.lat, self.lon, f, el,
                p, D, hs=hs, rho=rho_p, T=T, H=H, P=P)

#        print("\n")
#        print("- Rain attenuation          ", A_r, " dB")
#        print("- Gaseous attenuation       ", A_g, " dB")
#        print("- Cloud attenuation         ", A_c, " dB")
#        print("- Scintillation attenuation ", A_s, " dB")
#        print("- Total attenuation         ", A_t, " dB")

        print A_t
        print A_t.value

        self.k1 = self.db_to_lin(A_t.value)
        return self.k1

    def __str__(self):
        '''Lets print what we have done.'''
        compose = [
                "{:16} {:>11.4f} {:<13}".format('Distance', self.distance, ''),
                "{:16} {:>11.4f} {:<13}".format('Wavelength', self.wavelength, 'm'),
                "{:16} {:>11.4f} {:<13}".format('HPBW', self.hpbw, 'deg'),
                "{:16} {:>11.4f} {:<13}".format('T0', self.t0, 'K'),
                "{:16} {:>11.4f} {:<13}".format('T1/T0', self.t1overt0, ""),
                "{:16} {:>11.4f} {:<13} (apparent temperature of the moon)".format('Tmoon', self.t_moon, "K"),
                "{:16} {:>10g} {:<13} (Moon flux density)".format('Smoon', self.s_moon, "W m^-2 Hz^-1"),
                "{:16} {:>11.4f} {:<13} (lunar phase lag)".format('psi', self.psi, "deg"),
                "{:16} {:>11.4f} {:<13} (lunar phase angle)".format('phi', self.deg, "deg"),
                "{:16} {:>11.4f} {:<13} (lunar angular diameter)".format('Theta_moon', self.theta, "deg"),
                "{:16} {:>11.4f} {:<13} (atmospheric correction factor [linear])".format('K1', self.k1, ""),
                "{:16} {:>11.4f} {:<13} (extended source size correction factor [linear])".format('K2', self.k2est, ""),
                "\n Measurement Data",
                "{:16} {:>11.4f} {:<13}".format('Frequency', self.frequency, "Hz"),
                "{:16} {:>11.4f} {:<13}".format('Y-factor', self.y, ""),
                "{:16} {:>11.4f} {:<13}".format('Elevation', self.elevation, ''),
                "{:16} {:>11.4f} {:<13}".format('Azimuth', self.azimuth, ''),
                "{:16} {:>11.4f} {:<13} ([linear])".format('G/T', self.gt_lin, ""),
                "{:16} {:>11.4f} {:<13} ".format('G/T', self.gt_db, "dB/K"),
                "\n",

                  ]

        output = "\n".join(compose)
        return output

    def freq_input(self, freq):
        '''ingest the frequency fed in by the user. Convert whatever to Hz.'''
        if 'MHZ' in freq.upper():
            f_ = (freq.upper()).split('MHZ')
            f = float(f_[0])*10**6
            self.frequency = f
        elif 'GHZ' in freq.upper():
            f_ = (freq.upper()).split('GHZ')
            f = float(f_[0])*10**9
            self.frequency = f
        elif 'HZ' in freq.upper():
            f_ = (freq.upper()).split('HZ')
            f = float(f_[0])
            self.frequency = f
        elif freq.find('.') is 1:
            f = float(freq)*10**9
            self.frequency = f
            print "Labels are nice. I think you entered a frequency in GHz."
        elif freq.find('.') is 4:
            f = float(freq)*10**6
            self.frequency = f
            print "Labels are nice. I think you entered a frequency in MHz."
        elif freq.isdigit():
            f = float(freq)
            self.frequency = f
            print "Labels are nice. I think you entered a frequency in Hz."

        if self.frequency > 10*10**9:
            print str(self.frequency*10**-9)+" GHz is out of range. Try between 10 and 1 GHz."
            sys.exit()
        elif self.frequency < 1*10**9:
            print str(self.frequency*10**-9)+" GHZ is out of range. Try between 10 and 1 GHz."
            sys.exit()

        print str(f*10**-9)+" GHz is what I'll use."

    def y_input(self, y):
        '''Ingest the y-factor for later use. convert whatever to linear units[W].'''
        if 'DB' in y.upper():
            y_ = (y.upper()).split('DB')
            self.y = self.db_to_lin(float(y_[0]))
        elif y.isdigit():
            self.y = float(y)
            print "No Y-factor labels, it must have been a linear value."

    def temp_input(self, temp):
        '''store the temperature.'''
        self.temp = float(temp)

    def press_input(self, press):
        '''Store the pressure.'''
        self.press = float(press)

    def humid_input(self, humid):
        '''Store the humidity'''
        self.humid = float(humid)


if __name__ == "__main__":

    import sys

    arguments = len(sys.argv) - 1
    if arguments is 0:
        print "\n"
        print "This command requires a Y factor[W] and a Frequency [Hz] Temp[C] Pressure[hPA] Humidity[%] "
        print "For example 'figure_of_merit.py Y.YY FFFFFFFFFF TT PPPP HH'"
        sys.exit()
    if arguments is 1:
        merit = MoonGoT(test='yes')
        print merit
    if arguments is 2:
        y_fac = sys.argv[1]
        freq = sys.argv[2]
        merit = MoonGoT()
        merit.freq_input(freq)
        merit.y_input(y_fac)
        merit.calc_g_over_t()
        print merit
    if arguments is 5:
        merit = MoonGoT()
        print "\n"
        merit.y_input(sys.argv[1])
        merit.freq_input(sys.argv[2])
        merit.temp_input(sys.argv[3])
        merit.press_input(sys.argv[4])
        merit.humid_input(sys.argv[5])
        print "\n"

        merit.calc_g_over_t()

        print merit

#        got = MoonGoT()
#        got.calc_g_over_t(2.266, 8177500000)
#        print got

#        iot = MoonGoT(test='yes')
#        print iot

#g_over_t(2.266, 8177500000, 11.28)
