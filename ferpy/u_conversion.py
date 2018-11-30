#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 14:20:10 2017

@author: frico
"""

#==============================================================================
#                    
#                   Conversion btw units
#
#==============================================================================

from math import pi
import astropy.units as u
import astropy.constants.si as _si
import numpy as np

import astropy.cosmology
from astropy.cosmology import FlatLambdaCDM

# Cosmology
cosmo = FlatLambdaCDM(H0=69.6, Om0=0.286)



#__all__ = ["wm2_to_jykms", "lum_to_flux", "lsun_to_watt", "lsun_to_ergs", "Jybeam_to_T"]

# Luminostiy to Flux
def lum_to_flux(luminosity, z_source, z=True):
    """
    From Luminosity in Lsun to flux in W/m^2
    F = L/(4pid^2)
    if z is false, z_source is dist in Mpc
    """
    if z==True:
        dist_lum = cosmo.luminosity_distance(z_source)
        dist_lum = dist_lum.to(u.m) # dist lum in meters
    else:
        dist_lum = u.Quantity(z_source, u.Mpc)
        dist_lum = dist_lum.to(u.m)
    lum = u.Quantity(luminosity, u.L_sun)
    lum = lsun_to_watt(lum)
    flux = lum/(4.0*pi*dist_lum**2)
    return flux # Flux in W/m^2

def flux_to_lum(flux, z_source):
    """
    From flux in W/m^2 to Luminosity in Lsun
    F = L/(4pid^2)
    """
    dist_lum = cosmo.luminosity_distance(z_source)
    dist_lum = dist_lum.to(u.m) # dist lum in meters
    flux = u.Quantity(flux, u.W / u.m**2)
    lum = (4.0*pi*dist_lum**2) * flux
    lum = lum.to(u.L_sun)
    return lum # Luminosity in Lsun

def fluxwcm2_to_lum(flux, dist_lum):
    """
    From flux in W/cm^2 to Luminosity in Lsun
    F = L/(4pid^2)
    """
    flux = u.Quantity(flux, u.W / u.cm**2)
    flux = flux.to(u.W/u.m**2)
    dist_lum = dist_lum.to(u.m) # dist lum in meters
    lum = (4.0*pi*dist_lum**2) * flux
    lum = lum.to(u.L_sun)
    return lum # Luminosity in Lsun

# Stefan-Boltzmann Law
def stef_boltz(r,T):
    """
    Stefan-Boltzmann law
    r must be given in meters
    T in Kelvins
    returns luminostiy in W
    """
    s_si = _si.sigma_sb
    T = u.Quantity(T, u.K)
    r = u.Quantity(r, u.m)
    L = 4*pi*(r**2)*s_si*T**4
    return L

# Stefan-Boltzmann Law Error
def stef_boltz_error(r, r_err, T, T_err):
    """
    Stefan-Boltzmann law
    r must be given in meters
    T in Kelvins
    returns luminostiy in W
    """
    s_si = _si.sigma_sb.value
    T = u.Quantity(T, u.K).value
    r = u.Quantity(r, u.m).value
    L_err = np.sqrt(((4*4*pi*(r**2)*s_si*(T**3)*T_err)**2)+((4*pi*2*(r)*s_si*(T**4)*r_err)**2))
    return L_err

# Planck Function
def planck(freq, T):
    """ Maaaaaaaaaal
    Planck Function
    T in kelvins
    freq in Hz
    """
    h_si = _si.h
    c_si = _si.c
    k_si = _si.k_B
    T = u.Quantity(T, u.K)
    freq = u.Quantity(freq, u.Hz)
    B = 2.*h_si*freq**3./(c_si**2. * np.exp(h_si*freq/(k_si*T)) - (1.0 *u.m**2 /u.s**2))
    return B
    
# Linear Size
def lin_size(D, r):
    """
    Linear size from angular size
    D is distance in Mpc
    r is radius in arcsec
    returns the linear size in meters
    """
    D = u.Quantity(D, u.Mpc)
    d = D.to(u.m)
    L = r*d/206265.0
    return L

# Angular Size
def ang_size(D, L):
    """
    Angular size from linear size
    D is distance in Mpc
     linear size in meters
    returns the angular size in arcsec
    """
    D = u.Quantity(D, u.Mpc)
    L = u.Quantity(L, u.m)
    d = D.to(u.m)
    r = (u.arcsec)*L*206265.0/d
    return r

# Flux conversions
def wm2_to_jykms(flux, freq_wave_obs, *wavelength):
    """
    From W/m^2 to Jy km/s
    Intensity (W/m^2)
    freq_wave_obs (GHz or mm)
    c_si (km/s)
    wavelength True --> fobs is lambda_obs
    """
    c_si = _si.c.to(u.km / u.s) # km / s
    
    flux = u.Quantity(flux, u.W/ u.m**2)
    
    freq_obs = u.Quantity(freq_wave_obs, u.GHz)
    fobs = freq_obs.to(u.Hz, equivalencies=u.spectral())
    
    if wavelength:
        wave_obs = u.Quantity(freq_wave_obs, u.mm)
        fobs = wave_obs.to(u.Hz, equivalencies=u.spectral())
    
    dflux = flux * c_si / fobs #  W * km / s * Hz * m^2
    dflux_kms = dflux.to(u.Jy * u.km / u.s) # Jy * km / s
    return dflux_kms # Flux in Jy km/s
    
# Flux conversions
def wm2_to_jy(flux, fwhm, *wavelength):
    """
    From W/m^2 to Jy km/s
    Intensity (W/m^2)
    freq_wave_obs (GHz or mm)
    c_si (km/s)
    wavelength True --> fobs is lambda_obs
    """
    
    flux = u.Quantity(flux, u.W/ u.m**2)
    
    freq_obs = u.Quantity(fwhm, u.GHz)
    fobs = freq_obs.to(u.Hz, equivalencies=u.spectral())
    
    if wavelength:
        wave_obs = u.Quantity(fwhm, u.mm)
        fobs = wave_obs.to(u.Hz, equivalencies=u.spectral())
    
    dflux = flux / fobs #  W / Hz * m^2
    dflux_k = dflux.to(u.Jy) # Jy
    return dflux_k # Flux in Jy
    
def jy_to_wm2(flux, fwhm, *wavelength):
    """
    From W/m^2 to Jy km/s
    Intensity (W/m^2)
    freq_wave_obs (GHz or mm)
    c_si (km/s)
    wavelength True --> fobs is lambda_obs
    """
    
    flux = u.Quantity(flux, u.Jy)
    
    freq_obs = u.Quantity(fwhm, u.GHz)
    fobs = freq_obs.to(u.Hz, equivalencies=u.spectral())
    
    if wavelength:
        wave_obs = u.Quantity(fwhm, u.mm)
        fobs = wave_obs.to(u.Hz, equivalencies=u.spectral())
    
    dflux = flux * fobs #  W / Hz * m^2
    dflux_k = dflux.to(u.W/ u.m**2) # Jy
    return dflux_k # Flux in Jy

def jykms_to_wm2(flux, freq_wave_obs, *wavelength):
    """
    From Jy km/s to W/m^2
    Intensity (W/m^2)
    freq_wave_obs (GHz or mm)
    c_si (km/s)
    wavelength True --> fobs is lambda_obs
    """
    c_si = _si.c.to(u.km / u.s) # km / s
    
    flux = u.Quantity(flux, u.Jy * u.km / u.s)
    
    freq_obs = u.Quantity(freq_wave_obs, u.GHz)
    fobs = freq_obs.to(u.Hz, equivalencies=u.spectral())
    
    if wavelength:
        wave_obs = u.Quantity(freq_wave_obs, u.mm)
        fobs = wave_obs.to(u.Hz, equivalencies=u.spectral())
    
    dflux = flux * fobs / c_si #  Jy * Hz
    dflux_w = dflux.to(u.W / u.m**2) # W/m^2
    return dflux_w # Flux in Jy km/s

def jykms_to_wm2_w(flux, wave_obs):
    """
    From Jy km/s to W/m^2
    Intensity (W/m^2)
    freq_wave_obs (GHz or mm)
    c_si (km/s)
    wavelength True --> fobs is lambda_obs
    """
    c_si = _si.c.to(u.km / u.s) # km / s
    
    flux = u.Quantity(flux, u.Jy * u.km / u.s)
    
    wave_obs = u.Quantity(wave_obs, u.um)
    fobs = wave_obs.to(u.Hz, equivalencies=u.spectral())
    
    dflux = flux * fobs / c_si #  Jy * Hz
    dflux_w = dflux.to(u.W / u.m**2) # W/m^2
    return dflux_w # Flux in Jy km/s

def wm_to_ergscm(flux):
    """
    From W/m^2 to erg/s/cm^2
    """
    flux = u.Quantity(flux, u.W/u.m**2)
    dflux = flux.to(u.erg*u.s**(-1)*u.cm**(-2))
    return dflux # Flux in erg/s/cm^2

def wcm2_to_ergscm2(flux):
    """
    From W/cm^2 to erg/s/cm^2
    """
    flux = u.Quantity(flux, u.W/u.cm**2)
    dflux = flux.to(u.erg*u.s**(-1)*u.cm**(-2))
    return dflux # Flux in erg/s/cm^2

def ergscm_to_wm(flux):
    """
    From erg/s/cm^2 to W/m^2
    """
    flux = u.Quantity(flux, u.erg*u.s**(-1)*u.cm**(-2))
    dflux = flux.to(u.W/u.m**2)
    return dflux # Flux in W/m^2

def ergscm_to_jy(flux):
    """
    From erg/s/cm^2 to jy
    """
    flux = u.Quantity(flux, u.erg*u.s**(-1)*u.cm**(-2))
    dflux = flux.to(u.Jy)
    return dflux # Flux in Jy


def ergscmum_to_mjy(flux, lambda_microns):
    """
    From erg/s/cm^2/um to mJy
    """
    #flux = u.Quantity(flux, u.erg*u.s**(-1)*u.cm**(-2)*(u.um))
    #lambda_microns = u.Quantity(lambda_microns*(u.um))
    dflux = flux*(lambda_microns**2)*1e26/3e14
    return dflux # Flux in mJy

def jy_to_ergscm(flux):
    """
    From erg/s/cm^2 to jy
    """
    flux = u.Quantity(flux, u.Jy)
    dflux = flux.to(u.erg*u.s**(-1)*u.cm**(-2))
    return dflux # Flux in erg/s/cm^2

# Luminosities conversions
def lsun_to_watt(luminosity):
    """
    From luminostiy in Lsun to W
    """
    lum = u.Quantity(luminosity, u.L_sun)
    return lum.to(u.W)

def watt_to_lsun(luminosity):
    """
    From luminostiy in W to Lsun
    """
    lum = u.Quantity(luminosity, u.W)
    return lum.to(u.L_sun)

def lsun_to_ergs(luminosity):
    """
    From luminostiy in Lsun to erg/s
    """
    lum = u.Quantity(luminosity, u.L_sun)
    return lum.to(u.erg / u.s)

def ergs_to_lsun(luminosity):
    """
    From luminostiy in erg/s to Lsun
    """
    lum = u.Quantity(luminosity, u.erg / u.s)
    return lum.to(u.L_sun)

# Flux to luminosity
def wm2_to_lsun(flux, dist_lum):
    """
    From luminostiy in wm2 to Lsun
    distance luminosity requiered in Mpc
    """
    d_lum = u.Quantity(dist_lum, u.Mpc)
    f = u.Quantity(flux, u.W / u.m**2)
    lum = 4.*np.pi*(d_lum.to(u.m)**2.)*f
    return lum.to(u.L_sun)

# Jy/beam to Brightness Temperature
def Jybeam_to_T(Jybeam, freq, Bmin, Bmaj):
    # https://science.nrao.edu/facilities/vla/proposing/TBconv
    """
    Flux density in Jy/Beam yo Brightness temperature in Kelvin
    freq in GHz
    Bmin and Bmaj in arcseconds
    """
    T = 1.222E3 * (Jybeam * 1000.) / ((freq**2.) * Bmin *Bmaj)
    return T



