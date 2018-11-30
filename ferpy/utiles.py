#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 12:41:35 2017

@author: frico
"""


#==============================================================================
#                    
#                   Utils
#
#==============================================================================

import numpy as np
from decimal import Decimal
import re
from scipy import optimize
from astropy.modeling import models, fitting
import random as rnd
from copy import deepcopy
import pandas as pd
import networkx as nx
import astropy.units as u
import astropy.constants.si as _si
from astropy.modeling.blackbody import blackbody_nu

import matplotlib as mpl
#mpl.use('Agg') # Disables displaying plots
import matplotlib.pyplot as plt
#plt.ioff() # Turn off the interactive mode

# To print numbers in scientific form
def format_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

def dif_btw_listelemnt(lista):
    return [j-i for i, j in zip(lista[:-1], lista[1:])]

def ratio_btw_listelemnt(lista):
    return [j/i for i, j in zip(lista[:-1], lista[1:])]

def mean_btw_listelemnt(lista):
    return [(j+i)/2. for i, j in zip(lista[:-1], lista[1:])]

def checkEqual(iterator):
    '''
    Checks if all elements in a list are the same
    '''
    return len(set(iterator)) <= 1

# Find string between characters
def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

def rounding_exp(a, decimals):
    decimals = 1
    b = np.log10(a)
    c = 10.**(b-int(b))
    d = np.round(a, abs(int(b))+1+decimals)     
    return d
        
def binnings(lista_x, lista_y, binsize):
    # el rango lo hago -zbinsize y +2*zbinsize para que al pitnarlo con step no se me corte la anhura del "bin"
    bins = np.arange(min(lista_x)-binsize, max(lista_x)+binsize+binsize, binsize)
    binned = [] # binned list
    for i, val in enumerate(bins):
        binned.append([])
        for j, lob in enumerate(lista_x):
            if val-(binsize/2.0) <= lob <= val+(binsize/2.0):
                binned[i].append(binned[j])
    return (bins, binned)

def alma_cell_size(longest_base_in_wavelengths):
    # From https://casaguides.nrao.edu/index.php/Image_Continuum
    longw = longest_base_in_wavelengths
    resolution_arcsec = 206265.0/longw
    cellsize_arcsec = resolution_arcsec/7
    return cellsize_arcsec

def val_is_outrange(lista_x, lista_y, range_min, range_max):
    inrange_x = []
    inrange_y = []
    for i, val in enumerate(lista_x):
        if val <= range_min or val >= range_max:
            inrange_x.append(lista_x[i])
            inrange_y.append(lista_y[i])
    return (inrange_x, inrange_y)
    
def val_is_inrange(lista_x, lista_y, range_min, range_max):
    inrange_x = []
    inrange_y = []
    for i, val in enumerate(lista_x):
        if val >= range_min and val <= range_max:
            inrange_x.append(lista_x[i])
            inrange_y.append(lista_y[i])
    return (inrange_x, inrange_y)

def HMS2deg(ra='', dec=''):
  RA, DEC, rs, ds = '', '', 1, 1
  if dec:
    D, M, S = [float(i) for i in dec.split()]
    if str(D)[0] == '-':
      ds, D = -1, abs(D)
    deg = D + (M/60) + (S/3600)
    DEC = '{0}'.format(deg*ds)
  
  if ra:
    H, M, S = [float(i) for i in ra.split()]
    if str(H)[0] == '-':
      rs, H = -1, abs(H)
    deg = (H*15) + (M/4) + (S/240)
    RA = '{0}'.format(deg*rs)
  
  if ra and dec:
    return (RA, DEC)
  else:
    return RA or DEC

def get_numbers_from_filename(filename):
    '''
    Return numbers in filename
    '''
    return re.search(r'\d+', filename).group(0)

def deg2HMS(ra='', dec='', round=False):
  RA, DEC, rs, ds = '', '', '', ''
  if dec:
    if str(dec)[0] == '-':
      ds, dec = '-', abs(dec)
    deg = int(dec)
    decM = abs(int((dec-deg)*60))
    if round:
      decS = int((abs((dec-deg)*60)-decM)*60)
    else:
      decS = (abs((dec-deg)*60)-decM)*60
    DEC = '{0}{1} {2} {3}'.format(ds, deg, decM, decS)
  
  if ra:
    if str(ra)[0] == '-':
      rs, ra = '-', abs(ra)
    raH = int(ra/15)
    raM = int(((ra/15)-raH)*60)
    if round:
      raS = int(((((ra/15)-raH)*60)-raM)*60)
    else:
      raS = ((((ra/15)-raH)*60)-raM)*60
    RA = '{0}{1} {2} {3}'.format(rs, raH, raM, raS)
  
  if ra and dec:
    return (RA, DEC)
  else:
    return RA or DEC

def ang_distance_btw_2points(ra1, dec1, ra2, dec2):
    '''
    Angular distance between two points
    '''
    dec1 = np.deg2rad(dec1)
    dec2 = np.deg2rad(dec2)
    ra1 = np.deg2rad(ra1)
    ra2 = np.deg2rad(ra2)
    
    num = np.sqrt(((np.cos(dec2))**2)*((np.sin(ra2-ra1))**2) + (np.cos(dec1)*np.sin(dec2) - np.sin(dec1)*np.cos(dec2)*np.cos(ra2-ra1))**2)
    den = np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra2-ra1)
    
    
    dist = (180./np.pi)*np.arctan(num/den)
    #cosdist = np.cos(90-dec1)*np.cos(90-dec2)+np.sin(90-dec1)*np.sin(90-dec2)*np.cos(ra1-ra2)
    #dist = np.arccos(cosdist)
    return dist

def ang_distance_btw_2points_v2(ra1, dec1, ra2, dec2):
    '''
    Angular distance between two points
    '''
    dec1 = np.deg2rad(dec1)
    dec2 = np.deg2rad(dec2)
    ra1 = np.deg2rad(ra1)
    ra2 = np.deg2rad(ra2)
    cosdist = np.cos(90-dec1)*np.cos(90-dec2)+np.sin(90-dec1)*np.sin(90-dec2)*np.cos(ra1-ra2)
    dist = np.arccos(cosdist)
    return dist

def fit_bootstrap(p0, datax, datay, function, yerr_systematic=0.0):

    errfunc = lambda p, x, y: function(x,p) - y

    # Fit first time
    pfit, perr = optimize.leastsq(errfunc, p0, args=(datax, datay), full_output=0)


    # Get the stdev of the residuals
    residuals = errfunc(pfit, datax, datay)
    sigma_res = np.std(residuals)

    sigma_err_total = np.sqrt(sigma_res**2 + yerr_systematic**2)

    # 100 random data sets are generated and fitted
    ps = []
    for i in range(100):

        randomDelta = np.random.normal(0., sigma_err_total, len(datay))
        randomdataY = datay + randomDelta

        randomfit, randomcov = \
            optimize.leastsq(errfunc, p0, args=(datax, randomdataY),\
                             full_output=0)

        ps.append(randomfit) 

    ps = np.array(ps)
    mean_pfit = np.mean(ps,0)

    # You can choose the confidence interval that you want for your
    # parameter estimates: 
    Nsigma = 1. # 1sigma gets approximately the same as methods above
                # 1sigma corresponds to 68.3% confidence interval
                # 2sigma corresponds to 95.44% confidence interval
    err_pfit = Nsigma * np.std(ps,0) 

    pfit_bootstrap = mean_pfit
    perr_bootstrap = err_pfit
    return pfit_bootstrap, perr_bootstrap 

def gaussian_fit(datax, datay, mean_0, stddev_0):
    """
    Fitting to a Gaussian
    datax and datay must be np.arrays
    mean_0      -> mean initial guess
    stddev_0    -> stddev initial guess
    Returns [amplitude, mean, sigma, fwhm], [amplitude_err, mean_err, sigma_err, fwhm_err]
    """
    # Checking if data type is np.array or pd.Series
    if type(datax) not in (np.ndarray, np.array, pd.Series) or type(datay) not in (np.ndarray, np.array, pd.Series):
        raise ValueError('Input is not np.ndarray or pd.Series')
        
    # Defining initial gaussian and fitting
    g_init = models.Gaussian1D(amplitude=datay.max(), mean = mean_0, stddev=stddev_0)
    #g_init.stddev.bounds = 1.e-100, None
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, datax, datay)
    # Parameters and errors
    g_params = g.parameters
    cov_matrix = fit_g.fit_info['param_cov']
    if cov_matrix is None:
        g_errors = [np.nan, np.nan, np.nan]
    else:
        g_errors = np.sqrt(np.diag(cov_matrix))
    # Full width at half maximun
    g_fwhm = 2. * np.sqrt(2* np.log(2)) * g_params[2]
    g_fwhm_err = 2. * np.sqrt(2* np.log(2)) * g_errors[2]
    
    g_params = np.append(g_params, g_fwhm)
    g_errors = np.append(g_errors, g_fwhm_err)
    return g_params, g_errors, cov_matrix

def stdev_to_fwhm(std, std_err):
    g_fwhm = 2. * np.sqrt(2* np.log(2)) * std
    g_fwhm_err = 2. * np.sqrt(2* np.log(2)) * std_err
    return g_fwhm, g_fwhm_err

def fwhm_to_stdev(fwhm, fwhm_err):
    g_std = fwhm/ (2. * np.sqrt(2* np.log(2)))
    g_std_err = fwhm_err / (2. * np.sqrt(2* np.log(2)))
    return g_std, g_std_err

def gaussian_area(amplitude, stddev, amplitude_err, stddev_err):
    """
    Calculates de area unde a Gaussian and its error using:
        Area = Ampl * stddev * sqrt(2 * pi)
    """
    g_area = amplitude * stddev * np.sqrt(2*np.pi)
    g_area_err = np.sqrt(((stddev * np.sqrt(2*np.pi) * amplitude_err)**2.) +
                ((amplitude * np.sqrt(2*np.pi) * stddev_err)**2.))
    return g_area, g_area_err
    
def fit_g_bootstrap(datax, datay, data_rms, g_params, g_errors, nboost, seed):
    """
    For bootstrapping with fitted gaussian parameters
    datax       -> x original data
    datay       -> y original data
    g_params    -> amplitude, mean, stddev, fwhm (gaussian parameters)
    g_err       -> amplitude_err, mean_err, stddev_err, fwhm_err (gaussian parameters errors)
    Mirar como definir atributos
    """
    # Setting seed
    np.random.seed(seed)
    boot_param = []
    boot_error = []
    boot_amp = []
    boot_mean = []
    boot_std = []
    boot_w = []
    print '\t\tStarting Bootsrap'
    for i in range(nboost):
        resample = np.random.choice(range(len(datay)), len(datay), replace=True)
        x_boots = [datax[i] for i in resample]
        y_boots = [datay[i] for i in resample]
        
        # Montecarlo
        #randomDelta = np.random.normal(0., data_rms, len(y_boots))
        #y_resampled = y_boots + randomDelta
        y_resampled = np.random.normal(y_boots, data_rms)
        #print '\tSimul: ' + str(i)
        gg_params, gg_errors, g_cov_mat = gaussian_fit(np.array(x_boots), np.array(y_resampled), np.mean(x_boots), np.std(x_boots))
        boot_param.append(gg_params)
        boot_error.append(gg_errors)
        
        
        if not np.any(pd.isnull(gg_params)):
            boot_amp.append(gg_params[0])
            boot_mean.append(gg_params[1])
            boot_std.append(gg_params[2])
            boot_w.append(gg_params[3])
        else:
            continue
    print '\t\tBootsrap finished'
    
    #amp_min = g_params[0] - 3.*g_errors[0]
    bins_num = int(np.round(np.sqrt(nboost/2.), decimals=0))

    # Width of the Amplitud distribution
    (n_amp, bins_amp) = np.histogram(boot_amp, bins=bins_num)
    amp_params, amp_errors, amp_cov_mat = gaussian_fit(np.array(mean_btw_listelemnt(bins_amp)), n_amp, g_params[0], g_errors[0])

    # Width of the Mean distribution
    (n_mean, bins_mean) = np.histogram(boot_mean, bins=bins_num) 
    mean_params, mean_errors, mean_cov_mat = gaussian_fit(np.array(mean_btw_listelemnt(bins_mean)), n_mean, g_params[1], g_errors[1])
    
    #Width of the Stddev distribution
    (n_std, bins_std) = np.histogram(boot_std, bins=bins_num)
    std_params, std_errors, std_cov_mat = gaussian_fit(np.array(mean_btw_listelemnt(bins_std)), n_std, g_params[2], g_errors[2])
    
    #Width of the FWHM distribution
    (n_w, bins_w) = np.histogram(boot_w, bins=bins_num)
    w_params, w_errors, w_cov_mat = gaussian_fit(np.array(mean_btw_listelemnt(bins_w)), n_w, g_params[3], g_errors[3])

    # Parameter Results
    param_boots = [amp_params, mean_params, std_params, w_params]
    param_boots_err = [amp_errors, mean_errors, std_errors, w_errors]
    
    # Parameter Distributions
    amp_dist = [mean_btw_listelemnt(bins_amp), n_amp]
    mean_dist = [mean_btw_listelemnt(bins_mean), n_mean]
    std_dist = [mean_btw_listelemnt(bins_std), n_std]
    w_dist = [mean_btw_listelemnt(bins_w), n_w]
    param_dist = [amp_dist, mean_dist, std_dist, w_dist]
    param_boot_vals = [boot_amp, boot_mean, boot_std, boot_w]
    
    return param_boots, param_boots_err, param_dist, param_boot_vals
    
def fit_g_bootstrap_plotter(g_params, g_errors, param_boots, param_boots_err, param_dist, param_boot_vals, out_fig_dir, galaxy, line):
    # Disables displaying plots
    param_names = ['Amp', 'Mean', 'Std', 'FWHM']
    for i, par in enumerate(param_boots):
        fig = plt.figure()
        ax = fig.add_subplot((111))
        ax.plot(param_dist[i][0], param_dist[i][1], linewidth=0.5, drawstyle='steps-mid', color='k', label='') 
        b_gauss_amp = models.Gaussian1D(amplitude=param_boots[i][0], mean=param_boots[i][1], stddev=param_boots[i][2])
        ax.plot(param_dist[i][0], b_gauss_amp(param_dist[i][0]), linewidth=0.8, color='r', label='', linestyle='--' )
        ax.errorbar(g_params[i], 0.3*param_boots[i][0], xerr=g_errors[i], marker='|', markersize=3, color='b', elinewidth=0.5, capsize=1., capthick=0.6, label='No Boots')
        ax.errorbar(param_boots[i][1], 0.25*param_boots[i][0], xerr=param_boots[i][2], marker='|', markersize=3, color='r', elinewidth=0.5, capsize=1., capthick=0.6,label='Boots')
        ax.set_xlabel(param_names[i])
        plt.legend(loc='best', fontsize='xx-small', facecolor=None, frameon=False)
        plt.savefig(out_fig_dir+'/'+galaxy+'_'+line+'_boots_'+param_names[i]+'.png', bbox_inches='tight', transparent=True, dpi=600)
        plt.close()
        fig = plt.figure()
        ax = fig.add_subplot((111))
        ax.scatter(np.arange(1, len(param_boot_vals[i])+1.), param_boot_vals[i], label='', facecolors='none', edgecolors='k', s=3, linewidth=0.5)
        ax.errorbar(0.4*len(param_boot_vals[i]), g_params[i], yerr=g_errors[i], marker='.', markersize=3, color='b', elinewidth=0.5, capsize=1., capthick=0.6, label='No Boots')
        ax.errorbar(0.6*len(param_boot_vals[i]), param_boots[i][1], yerr=param_boots[i][2], marker='.', markersize=3, color='r', elinewidth=0.5, capsize=1., capthick=0.6,label='Boots')
        ax.set_ylabel(param_names[i])
        plt.legend(loc='best', fontsize='xx-small', facecolor=None, frameon=False)
        plt.savefig(out_fig_dir+'/'+galaxy+'_'+line+'_boots_simval_'+param_names[i]+'.png', bbox_inches='tight', transparent=True, dpi=600)
        plt.close()
        plt.close('all')
        
def repel_labels(ax, x, y, labels, k=0.01, size_font=5):
    G = nx.DiGraph()
    data_nodes = []
    init_pos = {}
    for xi, yi, label in zip(x, y, labels):
        data_str = 'data_{0}'.format(label)
        G.add_node(data_str)
        G.add_node(label)
        G.add_edge(label, data_str)
        data_nodes.append(data_str)
        init_pos[data_str] = (xi, yi)
        init_pos[label] = (xi, yi)

    pos = nx.spring_layout(G, pos=init_pos, fixed=data_nodes, k=k)

    # undo spring_layout's rescaling
    pos_after = np.vstack([pos[d] for d in data_nodes])
    pos_before = np.vstack([init_pos[d] for d in data_nodes])
    scale, shift_x = np.polyfit(pos_after[:,0], pos_before[:,0], 1)
    scale, shift_y = np.polyfit(pos_after[:,1], pos_before[:,1], 1)
    shift = np.array([shift_x, shift_y])
    for key, val in pos.items():
        pos[key] = (val*scale) + shift

    for label, data_str in G.edges():
        #ax.annotate(label,
        #            xy=pos[data_str], xycoords='data',
        #            xytext=pos[label], textcoords='data',
        #            arrowprops=dict(arrowstyle="-",
        #                            shrinkA=0, shrinkB=0,
        #                            connectionstyle="arc3", 
        #                            color='k',  lw=0.4), size=size_font)
        ax.annotate(label,
                    xy=pos[data_str], xycoords='data',
                    xytext=pos[label], textcoords='offset pixels',
                    size=size_font)
    # expand limits
    all_pos = np.vstack(pos.values())
    x_span, y_span = np.ptp(all_pos, axis=0)
    mins = np.min(all_pos-x_span*0.15, 0)
    maxs = np.max(all_pos+y_span*0.15, 0)
    ax.set_xlim([mins[0], maxs[0]])
    ax.set_ylim([mins[1], maxs[1]])
    
def transition_temperature(wavelength):
    """
    To get temperature of the transition in K
    Wavelength in micros
    T = h*f / kB
    """
    w = u.Quantity(wavelength, u.um)
    l = w.to(u.m)
    c = _si.c.to(u.m / u.s)
    h = _si.h.to(u.eV * u.s)
    kb = _si.k_B.to(u.eV / u.K)
    f = c/l
    t = h*f/kb
    return t

def numformatter(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
    
def SFRfromLFIR(LFIR):
    """
    Kennicut 1998
    To get Star formation rate from LFIR (8-1000um)
    LFIR in erg s-1
    SFR in Msun /year
    """
    SFR = 4.5E-44 * LFIR
    return SFR

def tau_from_T(Tobs, Tkin):
    """
    Line optical depth from observed temperature and excitation temperature in Kelvin
    """
    tau = -np.log(1.-(Tobs/Tkin))
    return tau


def sourcesize_from_fit2(Tmain_beam, Tex, Size_arcsec):
    # mal
    """
    To retrieve source size (diameter) from observed temperature, extication temperature and beam size
    TMB * Size = Tex * Source_size
    """
    #Source_size_arcsec = np.sqrt((Size_arcsec)*Tmain_beam/Tex)
    Source_size_arcsec = (Tmain_beam/Tex)*Size_arcsec/(1.-(Tmain_beam/Tex))
    return Source_size_arcsec

def sourcesize_from_fit(Tmain_beam, Tex, Size_arcsec):
    """
    To retrieve source size (diameter) from observed temperature, extication temperature and beam size
    TMB * Size = Tex * Source_size
    """
    Source_size_arcsec = np.sqrt((Size_arcsec)*(Tmain_beam/Tex))
    return Source_size_arcsec


def beam_size(bmin, bmaj):
    """
    Returns beam size
    """
    beam_size = np.pi * bmin * bmaj / (4. *np.log(2))
    return beam_size
   
class data:
    """
    To create a dictionary
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        
        
def virial_mass(FWHM, vel_disp):
    """
    Calculates the Vyrial Theorem based Dynamical mass 
    FWHM is the deconvolved size of the source in pc
    vel_disp is the velocity dispersion in km/s
    http://adsabs.harvard.edu/abs/2018arXiv180402083L
    Leroy et al 2018
    """
    M = 892. * FWHM * (vel_disp**2.)
    return M


####
# ALMA
####
def alma_resol(longest_baseline_in_wavelengths):
    """
    Alma resolution in arcsec 
    https://casaguides.nrao.edu/index.php/Image_Continuum
    """
    res = 206265.0/(longest_baseline_in_wavelengths)
    return res

def cellsize(resolution):
    """
    Alma cell size in arcsec given a certain resolution in arcsec
    """
    cellsize = resolution/7.
    return cellsize
        

#####
def tau_transition(N_up, B_ul, freq, FWHM, T):
    """
    Goldsmith & Langer 1999
    """
    h_si = _si.h
    k_si = _si.k_B
    tau = (h_si/FWHM)* N_up * B_ul * (np.exp(h_si*freq/(k_si*T))-1.)
    return tau

def antenna_temperature(wave, source_size, antenna_beam, brightness):
    """
    Goldsmith & Langer 1999
    """
    k_si = _si.k_B
    Ta = ((wave**2)/(2.*k_si))*(source_size/antenna_beam)*brightness
    return Ta
        
# SED fit from Perez-Beaupuits 2018
def SED_model(nu, Tdust, Mdust, phi, D):
    """
    Perez-Beaupuits 2018
    phi = beam area filling factor
    Bnu = planck function
    angular_size = source solid angle
    T = dust temperature
    Mdust = dust mass
    D = distance to the source
    phi_cold = filling factor of the coldest component
    kd = absorption coefficient
    nu = freq in GHz
    beta = 2
    """
    D = D * u.Mpc
    nu = u.Quantity(nu, u.GHz)
    Tdust = Tdust * u.K
    Mdust = Mdust * u.Msun
    angular_size_arcsec = 17.3*9.2*(u.arcsec**2) #arcsec**2
    angular_size = angular_size_arcsec.to(u.sr)
    d_m = D.to(u.m)
    Mdust_kg = Mdust.to(u.kg)
    Tcmb = 2.73 # K
    phi_cold = 5.0e-1
    beta = 2
    kd =(u.m**2 / u.kg)* 0.04*(nu/(250.*u.GHz))**beta # m^2 / kg
    tau = kd*Mdust_kg/(angular_size.value*phi_cold*d_m**2)
    Bnu_T = blackbody_nu(nu, Tdust)
    Bnu_Tcmb = blackbody_nu(nu, Tcmb)
    Snu = (1. -np.exp(-tau))*(Bnu_T-Bnu_Tcmb)*angular_size*phi
    Snu_jy = Snu.to(u.Jy)
    return Snu_jy

# Energy of a transition
def trans_energy(freq):
    """
    return Energy in Kelvins of a transition given in frecuency
    E=h*nu
    freq in GHz
    """
    freq = 354.1
    nu = u.Quantity(freq, u.GHz)
    h_si = _si.h # Planck constant
    k_si = _si.k_B # boltzmann constant
    E = h_si*(nu.to(u.Hz))/u.Hz/u.s # Joules
    E_K = E/k_si
    return E_K.value

# Excitation Temperature (radiation excitation and collisonal excitation)
# Two level approach (Goldsmith1982)
def Tex_goldsmith(Tkin, Tcore, trans_freq, Crate, Arate, f):
    """
    Goldsmith1982
    Crate   -> Downward collison rate
    Arate   -> Downward stimulated transitions rate
    Tkin    -> Gas kinetic temperature
    Tcore   -> Source temperature
    f       -> Filling factor
    """
    # Equivalent temperature of the transition
    Tstar = trans_energy(trans_freq)
    coefs = Crate/Arate
    g = f*((np.exp(Tstar/Tcore)-1.)**-1)
    Tex = Tstar*(((Tstar/Tkin)+np.log((1.+g+coefs)/(g*np.exp(Tstar/Tkin)+coefs)))**-1)
    return Tex
    
# Excitation Temperature (purely radiation excitation C->0)
# Two level approach (Goldsmith1982)
def Tex_goldsmith_rad(Tkin, Tcore, trans_freq, f):
    """
    Goldsmith1982
    Crate   -> Downward collison rate
    Arate   -> Downward stimulated transitions rate
    Tkin    -> Gas kinetic temperature
    Tcore   -> Source temperature
    f       -> Filling factor
    """
    # Equivalent temperature of the transition
    Tstar = trans_energy(trans_freq)
    
    Tex = Tstar * (1./(np.log(np.exp(Tstar/Tcore)-1.+f)-np.log(f)))
    return Tex

# Excitation Temperature (purely collisional excitation g->0)
# Two level approach (Goldsmith1982)
def Tex_goldsmith_col(Tkin, trans_freq, Crate, Arate):
    """
    Goldsmith1982
    Crate   -> Downward collison rate
    Arate   -> Downward stimulated transitions rate
    Tkin    -> Gas kinetic temperature
    Tcore   -> Source temperature
    f       -> Filling factor
    """
    # Equivalent temperature of the transition
    Tstar = trans_energy(trans_freq)
    coefs = Arate/Crate
    print 'ncrit=%1.2E' % coefs
    Tex = Tstar /((Tstar/Tkin)+np.log(1.+coefs))
    # print Tex_goldsmith_col(300., 219.17, 6E-13, 6E-4)
    return Tex

def dens_from_Tex_col(Tkin, Tex, trans_freq, Crate, Arate):
    # Equivalent temperature of the transition
    Tstar = trans_energy(trans_freq)
    coefs = Arate/Crate
    dens = coefs/(np.exp((Tstar/Tex) - (Tstar/Tkin))-1)
    #print '%1.2E' % dens_from_Tex_col(300., 139., 219.17, 6E-10, 4E-4)
    #print Tex_goldsmith(300., 300., (45.E-6 *u.m).to(u.GHz, equivalencies=u.spectral()).value, 1E-13, 6E-4, 0.5)

    return dens

# Deesxcitation collisional cross section between vib. states (Goldsmith1982)
def vib_col(tkin):
    csec = (3E-12)*np.exp(-4.8*(tkin**-1./3))
    #n = (6E-4)/vib_col(300)
    #print 'ncrit=%1.2E' % n
    return csec

