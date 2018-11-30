#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 09:54:02 2018

@author: frico
"""

#==============================================================================
#  Plotting - HC3N*
#==============================================================================

version = '8'

# Things to plot
plot_panel_fig = False
plot_cont = True
plot_v7_over_cont = True
plot_CS_over_cont = False
plot_v0_over_cont = False
plot_v0_contours = False
plot_v7_contours = True
plot_rot_diag = False
plot_panel_fig2 = False
plot_HC_spec = False

# Calculate luminosities
luminosidades = False

import statcont as scont
import argparse

from astropy.io import fits
from astropy.wcs import WCS
from astropy import wcs
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.analytic_functions import blackbody_lambda, blackbody_nu
import astropy.constants.si as _si
 
import scipy
import numpy as np
import os
import pandas as pd


import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

from copy import deepcopy

from ferpy import clear_variables
from ferpy import redshift_scale
from ferpy import u_conversion
from ferpy import utiles
#from ferpy.molecules import HC3N_analysis

from astropy.wcs.utils import wcs_to_celestial_frame, custom_frame_mappings

from radio_beam import Beam

mpl.rc('xtick', color='k', direction='in', labelsize=6)
mpl.rc('ytick', color='k', direction='in', labelsize=6)
mpl.rc('xtick.major', size=6)
mpl.rc('ytick.major', size=6)


# Starting workdir
workdir =  os.getcwd()

# Data workdir
dworkdir_cubes = '/Users/frico/Documents/data/NGC253_H3O+'
dworkdir_spec = dworkdir_cubes+'/Hotcores_v4_all'


# Out dir
out_dir = dworkdir_spec+'/Results_v'+version+'/'
if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
# Out figura 1
out_fig1_dir = out_dir+'/figura1/'
if not os.path.exists(out_fig1_dir):
        os.makedirs(out_fig1_dir)

# Out figura 2
out_fig2_dir = out_dir+'/figura2/'
if not os.path.exists(out_fig2_dir):
        os.makedirs(out_fig2_dir)
        
# Out figura 3
out_fig3_dir = out_dir+'/figura3/'
if not os.path.exists(out_fig3_dir):
        os.makedirs(out_fig3_dir)
# Out figura 4
out_fig4_dir = out_dir+'/figura4/'
if not os.path.exists(out_fig4_dir):
        os.makedirs(out_fig4_dir)
        
        
# Out figura 5
out_fig5_dir = out_dir+'/figura5/'
if not os.path.exists(out_fig5_dir):
        os.makedirs(out_fig5_dir)
        
        
class data:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        
# HC3N properties
hc3nv7 = data(lmin='hc3nv7',
              freq=[218.8608*u.GHz, 219.173757*u.GHz, 364.676275*u.GHz, 365.195232*u.GHz],
              lgint=[-1.76, -1.759, -1.419, -1.418],
              elo=[306.961, 307.081, 460.252, 460.590]
              )
              

#==============================================================================
# Tau_dust and Lum
#==============================================================================

D = 3.5 # Mpc to NGC253
vsis = 258 #km/s

nucleus_RA = '00:47:32.94' #https://arxiv.org/pdf/1509.00330.pdf
nucleus_Dec = '-25:17:19.70'

# kinematic center
nucleus_RA2=  '00:47:33.14' #http://iopscience.iop.org/article/10.1088/0004-637X/716/2/1166/pdf
nucleus_Dec2= '-25:17:17.52'



#https://watermark.silverchair.com/392-1-L16.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAekwggHlBgkqhkiG9w0BBwagggHWMIIB0gIBADCCAcsGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMp1g3wVjMYtKjcZHLAgEQgIIBnLqf7nou770j3y8xyjshRH-NKyynu0IfmHFVorbk1WIDvlL1wf3MdHWhYseESqzTz7RPpVLrzDQxWCwQQL0z-oddXqYqazbrThu69DPIZCq_BHFcNlKWT3N6FTSqNBBfpp6nwFl7bSNST5uMFarsflzDUEWwaG6y1lLH5zxcHVNteMFn5JLvGgqi4hshT_4ANtRh8o1SAjsHxHGYxl3LUJE8ucc0arlyeuhBQxt8gNif3NYilz-ZjcM6kuqcVQg_MfQdJnD8zBn-9E5sINvYlT9vkOA6qCWILi4_KHU-MlOSM7q1t-KQBa8OGesDd-oSNd26oTbvb_AxBH2NnJGVlqIVCCAsIGuZkwJMMo8W5SptP5P2BTl42vYHqlHoxYD_dIASU_iaNhFt0TiFb1gG6IEwc09gEqLkjl55T1aS3wT32Qrc8a4T9uxLH7s4S12Q5Bp3R9Blsvs1Ad3r31Rct8mXshsfNO_hF5OIpacCl4ZhVafBtcuGdxr7aVUwPnn2ZaV7JbJrqDwKKxCczlXIhpEhoRjv1P7LhtcluCg
# TH2 #http://iopscience.iop.org/article/10.1088/0004-637X/716/2/1166/pdf
TH2_RA = '00:47:33.179'
TH2_Dec = '-25:17:17.13'

IR32_RA = '00:47:33.4'
IR32_Dec= '-25:17:14.0'

#TH6 http://iopscience.iop.org/article/10.1086/304739/pdf
TH6_RA = '00:47:33.10653'
TH6_Dec = '-25:17:17.9405'
TH7_RA = '00:47:33.10942'
TH7_Dec = '-25:17:19.2206'

#SSC in NGC253 Watson et al 1996
# http://adsabs.harvard.edu/abs/1996AJ....112..534W
bright_blob_RA = '00:47:32.964'
bright_blob_Dec= '-25:17:19.26'
a_RA = '00:47:33.410'
a_Dec= '-25:17:13.60'
i_RA = '00:47:33.013'
i_Dec= '-25:17:18.12'
h_RA = '00:47:33.063' # es n en realidad
h_Dec= '-25:17:17.22'

# Alejandro
a10_RA = '00:47:33.155'
a10_Dec= '-25:17:17.12'
a11_RA = '00:47:33.172'
a11_Dec= '-25:17:17.42'
a13_RA = '00:47:33.193'
a13_Dec= '-25:17:16.77'

# Ontiveiros (a ojo, no ponen posiciones!!)
# Strongest IR peak: knot4 = TH7
O4_RA = '00:47:33.10942'
O4_Dec = '-25:17:19.2206'
O5_RA = '00:47:32.9333'
O5_Dec= '-25:17:19.5'

O16_RA = '00:47:33.13333'
O16_Dec= '-25:17:17.75'
O18_RA = '00:47:33.081333174'
O18_Dec= '-25:17:17.57'


# Positions (me and Leroy2018)
hc3n_positions = pd.read_csv(dworkdir_spec+'/HC_positions_REAL.txt', delim_whitespace= True, header=0, comment='#')
    
                             
# Free-Free  (Alejandro positions)
H_positions = pd.read_csv(dworkdir_spec+'/Alejandro_positions.txt', delim_whitespace= True, header=0, comment='#')
                                 

# Observation data
hc3n_prop = pd.read_csv(dworkdir_spec+'/tabla_hc3n_map_plot.txt', delim_whitespace= True, header=0, comment='#')
    
# Table v9 column names
hc3n_prop.columns = ['Ind', 'Ind_2', 'RA',	'Dec', 'VLSR', 'VLSR_err', 'v', 'v_err', 
                     'N(HC3N)', 'N(HC3N)_err', 'N(HC3N)_LR','N(HC3N)_LR_err', 'N(rot)','N(rot)_err',
                     'Tvib', 'Tvib_err', 'Tvib_LR', 'Tvib_LR_err', 'Trot', 'Trot_err',
                     'Size', 'hc3nv0_peak_mJy/beam', 'hc3nv7_peak_mJy/beam', 'hc3nv6_peak_mJy/beam',
                     'tau_v0', 'tau_v0_err', 'tau_v7', 'tau_v7_err', 'tau_v6', 'tau_v6_err']

# Table v9 Index
hc3n_prop['Ind_ok'] = ['14', '13', '10', '11', '12', '8', '9', '5', '4', '3', '2', '1']
    
hc3n_prop['Size_m'] = u_conversion.lin_size(D,hc3n_prop['Size'])
hc3n_prop['Radio_m'] = hc3n_prop['Size_m']/2.

header_cont = fits.getheader(dworkdir_cubes+'/MAD_CUB_NGC253_TE_219_220_0.19X0.29_briggs_v7.pbcor.fits')#'/MAD_CUB_NGC253_TE_cont_218_HR_briggs.pbcor.fits') 
beam_cont = Beam.from_fits_header(header_cont) 
bmin = beam_cont.minor.to(u.arcsec)
bmaj = beam_cont.major.to(u.arcsec)


# Brightness temperature
T_B = u_conversion.Jybeam_to_T(0.046, 218.373, bmin.value, bmaj.value)

Source_size = utiles.sourcesize_from_fit(T_B, 342, 0.1)


# Line Luminosities
#hc3n_prop['L_Watt'] = u_conversion.stef_boltz(hc3n_prop['Radio_m'], hc3n_prop['Tex'])
#hc3n_prop['L_Lsun'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt'])


tam_m = u_conversion.lin_size(3.5,0.03)
l_w = u_conversion.stef_boltz(tam_m, 342)
l_s = u_conversion.watt_to_lsun(l_w)

Jybeam = 0.046
freq = 218.373




# Optical depths
tau = utiles.tau_from_T(T_B, hc3n_prop['Tvib'])


        
#==============================================================================
# MADCUBA NGC253 (ALMA) Continuo a 218 GHz
#==============================================================================
        
os.chdir(dworkdir_spec)

cont_fits_name = 'MAD_CUB_CROP_NGC253_TE_cont_218_HR_briggs.pbcor.fits'
#cont_fits_name = 'MAD_CUB_NGC253_TE_cont_358_natural.pbcor.fits'
ngc253_cont_fits = fits.open(cont_fits_name)#'MAD_CUB_NGC253_TE_cont_358_natural.pbcor.fits')
#ngc253_cont_fits = fits.open('MAD_CUB_NGC253_TE_cont_358_natural.pbcor.fits')
ngc253_cont_data = ngc253_cont_fits[0]
ngc253_cont_shape = ngc253_cont_data.data.shape
ngc253_cont_header = ngc253_cont_data.header
#ngc253_cont_header['CDELT3'] = 1
ngc253_cont_header['CDELT4'] = 1
cont_stdev  = np.nanstd(ngc253_cont_data.data[0,0,:,:])
cont_max    = np.nanmax(ngc253_cont_data.data[0,0,:,:])
cont_min    = np.nanmin(ngc253_cont_data.data[0,0,:,:])

ngc253_cont_header['CTYPE3'] = 'FREQ' # En el header aparace 'FREQ-TOP'
wcs_1 = WCS(ngc253_cont_header)
wcs_1.wcs.ctype = ['RA---SIN', 'DEC--SIN', 'FREQ', 'STOKES']
wcs_1 = wcs_1.dropaxis(3)
wcs_1 = wcs_1.dropaxis(2)

#cont_fits_name = 'MAD_CUB_CROP_NGC253_TE_cont_218_HR_briggs.pbcor.fits'
cont358_fits_name = 'MAD_CUB_NGC253_TE_cont_358_natural.pbcor.fits'
ngc253_cont358_fits = fits.open(cont358_fits_name)#'MAD_CUB_NGC253_TE_cont_358_natural.pbcor.fits')
#ngc253_cont_fits = fits.open('MAD_CUB_NGC253_TE_cont_358_natural.pbcor.fits')
ngc253_cont358_data = ngc253_cont358_fits[0]
ngc253_cont358_shape = ngc253_cont358_data.data.shape
ngc253_cont358_header = ngc253_cont358_data.header
#ngc253_cont_header['CDELT3'] = 1
ngc253_cont358_header['CDELT4'] = 1
cont358_stdev  = np.nanstd(ngc253_cont358_data.data[0,0,:,:])
cont358_max    = np.nanmax(ngc253_cont358_data.data[0,0,:,:])
cont358_min    = np.nanmin(ngc253_cont358_data.data[0,0,:,:])

ngc253_cont358_header['CTYPE3'] = 'FREQ' # En el header aparace 'FREQ-TOP'
wcs_1_358 = WCS(ngc253_cont358_header)
wcs_1_358.wcs.ctype = ['RA---SIN', 'DEC--SIN', 'FREQ', 'STOKES']
wcs_1_358 = wcs_1_358.dropaxis(3)
wcs_1_358 = wcs_1_358.dropaxis(2)

# Create a new WCS object.  The number of axes must be set
# from the start
w = wcs.WCS(naxis=2)
# Some pixel coordinates of interest.
pixcrd = np.array([[0, 0], [50, 50], [10, 12], [100, 120]], np.float_)
# Convert pixel coordinates to world coordinates
world= wcs_1.wcs_world2pix(pixcrd, 1)
# Convert the same coordinates back to pixel coordinates.
pixcrd2 = wcs_1.wcs_pix2world(world, 1)

if plot_cont == True or plot_v7_over_cont==True:
    #from regions import CircleSkyRegion
    #from regions import EllipseSkyRegion
    from matplotlib.patches import Ellipse
    from matplotlib.patches import Rectangle

    
    
    #ellipse1 = EllipseSkyRegion()
    
    fig = plt.figure()
    ax1 = fig.add_subplot((111), aspect='equal', projection=wcs_1)
    ax1.tick_params(labelsize=6)
        
    cont_stdev  = np.nanstd(ngc253_cont_data.data[0,0,:,:])
    cont_max    = np.nanmax(ngc253_cont_data.data[0,0,:,:])
    cont_min    = np.nanmin(ngc253_cont_data.data[0,0,:,:])
    #levels = np.linspace(2*stdev, ngc253_cont_data.data[0,0,:,:].max(), 10)
    ax1.tick_params(direction='in', labelsize=3)
    
    
    pixsize = 0.03 #arcsec
    ell_ax1 = 0.292 #arcsec
    ell_ax2 = 0.197 # arcsec
    pa = 80.279
    r_ax1 = 0.35 
    r_ax2 = 0.35
    #r = Rectangle(xy=(140-r_ax2/pixsize/2, 175-r_ax2/pixsize/2), width=r_ax1/pixsize, height=r_ax2/pixsize, edgecolor='white', facecolor='white',
    #          transform=ax1.get_transform(wcs_1))
    #ax1.add_patch(r)
    c = Ellipse(xy=(185, 195), width=ell_ax1/pixsize, height=ell_ax2/pixsize, angle=90-pa, edgecolor='w', facecolor='w', linewidth=0.5,
              transform=ax1.get_transform(wcs_1))
    ax1.add_patch(c)
    
    plt.xlabel('RA (J2000)')
    plt.ylabel('Dec (J2000)', labelpad=-1)
    
    ax1.coords[0].set_major_formatter('hh:mm:ss.ss')
    ax1.coords[0].set_ticks(size=6, width=1, color='w', exclude_overlapping=True, number = 5)#, spacing=0.2 * u.arcsec,
    ax1.coords[1].set_ticks(size=6, width=1, color='w')
    ax1.coords[0].set_separator((r'$^{\rm{h}}$', r'$^{\rm{m}}$', r'$^{\rm{s}}$'))
    
    ax1.tick_params(labelsize=6)
    
    ax1.coords[0].display_minor_ticks(True)
    ax1.coords[1].display_minor_ticks(True)
    
    plt.imshow(ngc253_cont_data.data[0,0,:,:], origin='lower', vmax=cont_max, vmin=0.2*cont_stdev, cmap= cm.jet, interpolation="none")#, extent=(0,0,0,0))
    plt.xlim([160, 483])#[120, 523]) #170, 474
    plt.ylim([175, 462])#[155, 482]) #336, 453
    
    # Plotting the center
    pos_nucleus = utiles.HMS2deg(ra=nucleus_RA2.replace(':', ' '), dec=nucleus_Dec2.replace(':', ' '))
    px_n, py_n = wcs_1.wcs_world2pix(float(pos_nucleus[0]), float(pos_nucleus[1]), 1)
    ax1.plot(px_n, py_n, 'x', color='white')
    
    
    for hotcore in hc3n_prop.itertuples():
        pos = utiles.HMS2deg(ra=hotcore.RA.replace('_', ' '), dec=hotcore.Dec.replace('_', ' '))
        px, py = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        #ax1.plot(px, py, '.r')
        #ax1.text(px, py+py/40., str(hotcore.Ind), fontsize=10, ha='center', va='center', color='white')
        px_m = 0
        if hotcore.Ind in [4, 7, 9, 11]:
            py_m = -20
        elif hotcore.Ind==3:
            py_m = 35
            px_m = 20
        elif hotcore.Ind_2 ==11:
            px_m = -10
            py_m = -20
        else:
            py_m = 20
        ax1.annotate(str(hotcore.Ind_ok), xy=(px,py), xytext=(px+px_m,py+py_m),
                arrowprops={'arrowstyle': '-', 'color': 'w'}, va='center', color = 'white')
        
    pos_6and7 = [['00 47 33.01', '-25 17 19.42'], ['00 47 33.01', '-25 17 19.02']]
    for p, posdeg in enumerate(pos_6and7):
        px_m = 0
        py_m = 20
        if p == 0:
            ind = 6
            px_m = -10
            py_m = -20
        else:
            ind = 7
        pos = utiles.HMS2deg(ra=posdeg[0], dec=posdeg[1])
        px, py = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        
        ax1.annotate(str(ind), xy=(px,py), xytext=(px+px_m,py+py_m),
                arrowprops={'arrowstyle': '-', 'color': 'w'}, va='center', color = 'white')
    
    
        
        
    
    ax1.coords.frame.set_color('w') # Frame color
    ax1.coords.frame.set_linewidth(0)
    ax1.tick_params(labelsize=6)
    #fig.savefig(out_fig1_dir+'/ngc253_cont_plasma_v'+version+'.eps', bbox_inches='tight', transparent=True, dpi=600)
    plt.close()

#==============================================================================
# MADCUBA NGC253 (ALMA) Vel integrada HC3N v7=1
#==============================================================================

#### HC3N a 0.19" de resolucion

# HC3N v0
ngc253_hc3nv0_gr_fits = fits.open('MAD_CUB_CROP_II_NGC253_TE_HC3Nv0_20kms_218_219_0.19X0.29_briggs_v7.pbcor.fits')
ngc253_hc3nv0_gr_data = ngc253_hc3nv0_gr_fits[0]
ngc253_hc3nv0_gr_shape = ngc253_hc3nv0_gr_data.data.shape
ngc253_hc3nv0_gr_header = ngc253_hc3nv0_gr_data.header
#ngc253_hc3nv0_gr_header['CDELT3'] = 1
ngc253_hc3nv0_gr_header['CDELT4'] = 1

wcs_hc3nv0_gr_2 = WCS(ngc253_hc3nv0_gr_header)
wcs_hc3nv0_gr_2.wcs.ctype = ['RA---SIN', 'DEC--SIN', 'FREQ', 'STOKES']
wcs_hc3nv0_gr_2 = wcs_hc3nv0_gr_2.dropaxis(3)
#wcs_hc3nv0_gr_2 = wcs_hc3nv0_gr_2.dropaxis(2)

# HC3N v7=1
ints1 = 'MAD_CUB_CROP_II_NGC253_TE_HC3NV7_20KMS_219_220_0.19X0.29_briggs_v7.pbcor.fits'
ints2 = 'MAD_CUB_II_CROP_NGC253_TE_219_220_25KMS.pbcor.fits'
intsall=  'MAD_CUB_II_CROP_NGC253_TE_219_220_ALL.pbcor.fits'


ngc253_hc3nv7_gr_fits = fits.open(ints1)
ngc253_hc3nv7_gr_data = ngc253_hc3nv7_gr_fits[0]
ngc253_hc3nv7_gr_shape = ngc253_hc3nv7_gr_data.data.shape
ngc253_hc3nv7_gr_header = ngc253_hc3nv7_gr_data.header
#ngc253_hc3nv7_gr_header['CDELT4'] = 1
#ngc253_hc3nv7_gr_header['WCSAXES'] = 2
wcs_hc3nv7_gr_2 = WCS(ngc253_hc3nv7_gr_header, naxis=2)
#wcs_hc3nv7_gr_2.wcs.ctype = ['RA---SIN', 'DEC--SIN', 'FREQ', 'STOKES']
#wcs_hc3nv7_gr_2 = wcs_hc3nv7_gr_2.dropaxis(3)
#wcs_hc3nv7_gr_2 = wcs_hc3nv7_gr_2.dropaxis(2)

# HC3N v7=1
ints1_3938 = 'MAD_CUB_II_20kms_BL_NGC253_TE_354_356_0.30X0.25_briggs_v1.pbcor.fits'

ngc253_hc3nv7_gr_fits_3938 = fits.open(ints1_3938)
ngc253_hc3nv7_gr_data_3938 = ngc253_hc3nv7_gr_fits_3938[0]
ngc253_hc3nv7_gr_shape_3938 = ngc253_hc3nv7_gr_data_3938.data.shape
ngc253_hc3nv7_gr_header_3938 = ngc253_hc3nv7_gr_data_3938.header
#ngc253_hc3nv7_gr_header['CDELT4'] = 1
#ngc253_hc3nv7_gr_header['WCSAXES'] = 2
wcs_hc3nv7_gr_2_3938 = WCS(ngc253_hc3nv7_gr_header_3938, naxis=2)
#wcs_hc3nv7_gr_2.wcs.ctype = ['RA---SIN', 'DEC--SIN', 'FREQ', 'STOKES']
#wcs_hc3nv7_gr_2 = wcs_hc3nv7_gr_2.dropaxis(3)
#wcs_hc3nv7_gr_2 = wcs_hc3nv7_gr_2.dropaxis(2)


# Plotting HC3Nv7=1 over cont
if plot_v7_over_cont == True:
    from matplotlib.patches import Ellipse
    from matplotlib.patches import Rectangle
    fig = plt.figure()
    ax1 = fig.add_subplot((111), aspect='equal', projection=wcs_1)
    ax1.tick_params(direction='in')
    cm_max = np.nanmax(ngc253_cont_data.data[0,0,:,:])
    cm_min = np.abs(np.nanmin(ngc253_cont_data.data[0,0,:,:]))
    cm_std = np.nanstd(ngc253_cont_data.data[0,0,:,:])
    ax1.imshow(ngc253_cont_data.data[0,0,:,:], origin='lower', vmax=cont_max, vmin=0.2*cont_stdev, cmap=cm.gray, interpolation="none")# vmax=cm_max, vmin=cm_std*2, interpolation="none")
    plt.xlim([160, 483])#[120, 523]) #170, 474
    plt.ylim([175, 462])#[155, 482]) #336, 453
    plt.ylabel('Dec (J2000)')
    plt.xlabel('RA (J2000)')
    
    ax1.coords[0].set_major_formatter('hh:mm:ss.ss')
    ax1.coords[0].set_ticks(size=8,color='k', exclude_overlapping=True, number = 5, width=2)#, spacing=0.2 * u.arcsec,
    ax1.coords[1].set_ticks(size=8,color='k', width=2)
    ax1.coords[0].set_separator((r'$^{\rm{h}}$', r'$^{\rm{m}}$', r'$^{\rm{s}}$'))
    
    
    
    stdev_all_ngc253 = []
    max_all_ngc253 = []
    min_all_ngc253 = []
    #for i, vel in enumerate(colores_1):
    for i in range(ngc253_hc3nv7_gr_shape[1]): 
        #if i > 4:
            std = np.nanstd(ngc253_hc3nv7_gr_data.data[0,i,:,:])
            stdev_all_ngc253.append(std)
            maxx = np.nanmax(ngc253_hc3nv7_gr_data.data[0,i,:,:])
            max_all_ngc253.append(maxx)
            minn = np.nanmin(ngc253_hc3nv7_gr_data.data[0,i,:,:])
            min_all_ngc253.append(minn)
            levels_ngc253 = np.linspace(3.3*std, maxx, 3)
            
            ax1.contour(ngc253_hc3nv7_gr_data.data[0,i,:,:], colors='white', levels = levels_ngc253, linewidths=0.7, transform=ax1.get_transform(wcs_hc3nv7_gr_2))
    
    for hotcore in hc3n_prop.itertuples():
        pos = utiles.HMS2deg(ra=hotcore.RA.replace('_', ' '), dec=hotcore.Dec.replace('_', ' '))
        px, py = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        #ax1.plot(px, py, '.r')
        #ax1.text(px, py+py/40., str(hotcore.Ind), fontsize=10, ha='center', va='center', color='white')
        px_m = 0
        if hotcore.Ind in [4, 7, 9, 11]:
            py_m = -20
        elif hotcore.Ind==3:
            py_m = 35
            px_m = 20
        elif hotcore.Ind_2 ==11:
            px_m = -10
            py_m = -20
        else:
            py_m = 20
        ax1.annotate(str(hotcore.Ind_ok), xy=(px,py), xytext=(px+px_m,py+py_m),
                arrowprops={'arrowstyle': '-', 'color': 'w'}, va='center', color = 'white')
        
    pos_6and7 = [['00 47 33.01', '-25 17 19.42'], ['00 47 33.01', '-25 17 19.02']]
    for p, posdeg in enumerate(pos_6and7):
        px_m = 0
        py_m = 20
        if p == 0:
            ind = 6
            px_m = -10
            py_m = -20
        else:
            ind = 7
        pos = utiles.HMS2deg(ra=posdeg[0], dec=posdeg[1])
        px, py = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        
        ax1.annotate(str(ind), xy=(px,py), xytext=(px+px_m,py+py_m),
                arrowprops={'arrowstyle': '-', 'color': 'w'}, va='center', color = 'white')
        
        
    
    #fig.savefig(out_fig1_dir+'/ngc253_hc3nv7_219GHz_v'+version+'.eps', bbox_inches='tight', transparent=True, dpi=300)
    plt.close()
    
    
    # Other v7=1 at 355 GHz
    fig = plt.figure()
    ax1 = fig.add_subplot((111), aspect='equal', projection=wcs_1)
    ax1.tick_params(direction='in')
    cm_max = np.nanmax(ngc253_cont_data.data[0,0,:,:])
    cm_min = np.abs(np.nanmin(ngc253_cont_data.data[0,0,:,:]))
    cm_std = np.nanstd(ngc253_cont_data.data[0,0,:,:])
    ax1.imshow(ngc253_cont_data.data[0,0,:,:], origin='lower', vmax=cont_max, vmin=0.2*cont_stdev, cmap=cm.jet, interpolation="none")# vmax=cm_max, vmin=cm_std*2, interpolation="none")
    plt.xlim([160, 483])#[120, 523]) #170, 474
    plt.ylim([175, 462])#[155, 482]) #336, 453
    plt.ylabel('Dec (J2000)')
    plt.xlabel('RA (J2000)')
    
    ax1.coords[0].set_major_formatter('hh:mm:ss.ss')
    ax1.coords[0].set_ticks(size=8,color='k', exclude_overlapping=True, number = 5, width=2)#, spacing=0.2 * u.arcsec,
    ax1.coords[1].set_ticks(size=8,color='k', width=2)
    ax1.coords[0].set_separator((r'$^{\rm{h}}$', r'$^{\rm{m}}$', r'$^{\rm{s}}$'))
    
    
    
    stdev_all_ngc253_3938 = []
    max_all_ngc253_3938 = []
    min_all_ngc253_3938 = []
    #for i, vel in enumerate(colores_1):
    for i in range(ngc253_hc3nv7_gr_shape_3938[1]): 
        #if i > 4:
            std = np.nanstd(ngc253_hc3nv7_gr_data_3938.data[0,i,:,:])
            stdev_all_ngc253_3938.append(std)
            maxx = np.nanmax(ngc253_hc3nv7_gr_data_3938.data[0,i,:,:])
            max_all_ngc253_3938.append(maxx)
            minn = np.nanmin(ngc253_hc3nv7_gr_data_3938.data[0,i,:,:])
            min_all_ngc253_3938.append(minn)
            levels_ngc253_3938 = np.linspace(3.3*std, maxx, 3)
            
            ax1.contour(ngc253_hc3nv7_gr_data_3938.data[0,i,:,:], colors='white', levels = levels_ngc253_3938, linewidths=0.7, transform=ax1.get_transform(wcs_hc3nv7_gr_2_3938))
    
    for hotcore in hc3n_prop.itertuples():
        pos = utiles.HMS2deg(ra=hotcore.RA.replace('_', ' '), dec=hotcore.Dec.replace('_', ' '))
        px, py = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        #ax1.plot(px, py, '.r')
        #ax1.text(px, py+py/40., str(hotcore.Ind), fontsize=10, ha='center', va='center', color='white')
        px_m = 0
        if hotcore.Ind in [4, 7, 9, 11]:
            py_m = -20
        elif hotcore.Ind==3:
            py_m = 35
            px_m = 20
        elif hotcore.Ind_2 ==11:
            px_m = -10
            py_m = -20
        else:
            py_m = 20
        ax1.annotate(str(hotcore.Ind_ok), xy=(px,py), xytext=(px+px_m,py+py_m),
                arrowprops={'arrowstyle': '-', 'color': 'w'}, va='center', color = 'white')
        
    pos_6and7 = [['00 47 33.01', '-25 17 19.42'], ['00 47 33.01', '-25 17 19.02']]
    for p, posdeg in enumerate(pos_6and7):
        px_m = 0
        py_m = 20
        if p == 0:
            ind = 6
            px_m = -10
            py_m = -20
        else:
            ind = 7
        pos = utiles.HMS2deg(ra=posdeg[0], dec=posdeg[1])
        px, py = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        
        ax1.annotate(str(ind), xy=(px,py), xytext=(px+px_m,py+py_m),
                arrowprops={'arrowstyle': '-', 'color': 'w'}, va='center', color = 'white')
        
        
    
    #ax1.xaxis.set_tick_params(width=20)
    #ax1.yaxis.set_tick_params(width=20)
    #fig.savefig(out_fig1_dir+'/ngc253_hc3nv7_354GHz_v'+version+'.eps', bbox_inches='tight', transparent=True, dpi=300)
    plt.close()
    
    
    # Plotting HC3N J39-38 over 24-23
    ints1_2423 = 'MAD_CUB_II_ALL_NGC253_TE_219_220_0.19X0.29_briggs_v7.pbcor.fits'
    
    ngc253_hc3nv7_gr_fits_2423 = fits.open(ints1_2423)
    ngc253_hc3nv7_gr_data_2423 = ngc253_hc3nv7_gr_fits_2423[0]
    ngc253_hc3nv7_gr_shape_2423 = ngc253_hc3nv7_gr_data_2423.data.shape
    ngc253_hc3nv7_gr_header_2423 = ngc253_hc3nv7_gr_data_2423.header
    #ngc253_hc3nv7_gr_header['CDELT4'] = 1
    #ngc253_hc3nv7_gr_header['WCSAXES'] = 2
    wcs_hc3nv7_gr_2_2423 = WCS(ngc253_hc3nv7_gr_header_2423, naxis=2)
    
    fig = plt.figure()
    ax1 = fig.add_subplot((111), aspect='equal', projection=wcs_1)
    ax1.tick_params(direction='in')
    cm_max = np.nanmax(ngc253_hc3nv7_gr_data_2423.data[0,0,:,:])
    cm_min = np.abs(np.nanmin(ngc253_hc3nv7_gr_data_2423.data[0,0,:,:]))
    cm_std = np.nanstd(ngc253_hc3nv7_gr_data_2423.data[0,0,:,:])
    ax1.imshow(ngc253_cont_data.data[0,0,:,:], origin='lower', vmax=cont_max, vmin=0.2*cont_stdev, cmap=cm.jet, interpolation="none")# vmax=cm_max, vmin=cm_std*2, interpolation="none")
    if 'cont_218' in cont_fits_name:
        cont_name = 'cont_218'
        print cont_name

        plt.xlim([160, 483])#[120, 523]) #170, 474
        plt.ylim([175, 462])#[155, 482]) #336, 453
    elif 'cont_358' in cont_fits_name:
        cont_name = 'cont_358'
        print cont_name
        plt.xlim([150, 300])#[160, 483])#[120, 523]) #170, 474
        plt.ylim([160, 300])#[175, 462])#[155, 482]) #336, 453
    plt.ylabel('Dec (J2000)', labelpad=-1)
    plt.xlabel('RA (J2000)')
    
    pixsize_219 = 0.03 #arcsec
    ell_ax1_219 = 0.292 #arcsec
    ell_ax2_219 = 0.197 # arcsec
    pa_219 = 80.279
    
    pixsize_355 = 0.03 #arcsec es 0.05 pero estamos pintando sobre el continuo de 218 que tiene 0.03"
    ell_ax1_355 = 0.303 #arcsec
    ell_ax2_355 = 0.25 # arcsec
    pa_355 = -74.2

    c_219 = Ellipse(xy=(185, 195), width=ell_ax1_219/pixsize_219, height=ell_ax2_219/pixsize_219, angle=90-pa, edgecolor='w', facecolor='w', linewidth=0.5,
              transform=ax1.get_transform(wcs_1))
    c_355 = Ellipse(xy=(205, 195), width=ell_ax1_355/pixsize_355, height=ell_ax2_355/pixsize_355, angle=90-pa, edgecolor='r', facecolor='r', linewidth=0.5,
              transform=ax1.get_transform(wcs_1))
    ax1.add_patch(c_355)
    ax1.add_patch(c_219)
    
    ax1.coords[0].set_major_formatter('hh:mm:ss.ss')
    ax1.coords[0].set_ticks(size=6,color='w', exclude_overlapping=True, number = 5, width=1)#, spacing=0.2 * u.arcsec,
    ax1.coords[1].set_ticks(size=6,color='w', width=1)
    ax1.coords[0].set_separator((r'$^{\rm{h}}$', r'$^{\rm{m}}$', r'$^{\rm{s}}$'))
    
    ax1.coords[0].display_minor_ticks(True)
    ax1.coords[1].display_minor_ticks(True)
    
    stdev_all_ngc253 = []
    max_all_ngc253 = []
    min_all_ngc253 = []
    #for i, vel in enumerate(colores_1):
    for i in range(ngc253_hc3nv7_gr_shape[1]): 
        #if i > 4:
            std = np.nanstd(ngc253_hc3nv7_gr_data.data[0,i,:,:])
            stdev_all_ngc253_3938.append(std)
            maxx = np.nanmax(ngc253_hc3nv7_gr_data.data[0,i,:,:])
            max_all_ngc253_3938.append(maxx)
            minn = np.nanmin(ngc253_hc3nv7_gr_data.data[0,i,:,:])
            min_all_ngc253.append(minn)
            levels_ngc253 = np.geomspace(3.3*std, maxx, 3)
            
            ax1.contour(ngc253_hc3nv7_gr_data.data[0,i,:,:], colors='w', levels = levels_ngc253, linewidths=0.5, transform=ax1.get_transform(wcs_hc3nv7_gr_2))
    
    stdev_all_ngc253_3938 = []
    max_all_ngc253_3938 = []
    min_all_ngc253_3938 = []
    #for i, vel in enumerate(colores_1):
    for i in range(ngc253_hc3nv7_gr_shape_3938[1]): 
        #if i > 4:
            std = np.nanstd(ngc253_hc3nv7_gr_data_3938.data[0,i,:,:])
            stdev_all_ngc253_3938.append(std)
            maxx = np.nanmax(ngc253_hc3nv7_gr_data_3938.data[0,i,:,:])
            max_all_ngc253_3938.append(maxx)
            minn = np.nanmin(ngc253_hc3nv7_gr_data_3938.data[0,i,:,:])
            min_all_ngc253_3938.append(minn)
            levels_ngc253_3938 = np.geomspace(3.3*std, maxx, 3)
            ax1.contour(ngc253_hc3nv7_gr_data_3938.data[0,i,:,:], colors='r', levels = levels_ngc253_3938, linewidths=0.5, linestyles='solid', transform=ax1.get_transform(wcs_hc3nv7_gr_2_3938))
    

    for hotcore in hc3n_prop.itertuples():
        pos = utiles.HMS2deg(ra=hotcore.RA.replace('_', ' '), dec=hotcore.Dec.replace('_', ' '))
        px, py = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        #ax1.plot(px, py, '.r')
        #ax1.text(px, py+py/40., str(hotcore.Ind), fontsize=10, ha='center', va='center', color='white')
        px_m = 0
        if hotcore.Ind in [4, 7, 9, 11]:
            py_m = -20
        elif hotcore.Ind==3:
            py_m = 35
            px_m = 20
        elif hotcore.Ind_2 ==11:
            px_m = -10
            py_m = -20
        else:
            py_m = 20
        ax1.annotate(str(hotcore.Ind_ok), xy=(px,py), xytext=(px+px_m,py+py_m),
                arrowprops={'arrowstyle': '-', 'color': 'w'}, va='center', color = 'w')
        
    pos_6and7 = [['00 47 33.01', '-25 17 19.42'], ['00 47 33.01', '-25 17 19.02']]
    for p, posdeg in enumerate(pos_6and7):
        px_m = 0
        py_m = 20
        if p == 0:
            ind = 6
            px_m = -10
            py_m = -20
        else:
            ind = 7
        pos = utiles.HMS2deg(ra=posdeg[0], dec=posdeg[1])
        px, py = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        
        ax1.annotate(str(ind), xy=(px,py), xytext=(px+px_m,py+py_m),
                arrowprops={'arrowstyle': '-', 'color': 'w'}, va='center', color = 'w')
        
    # Plotting the kinematic center
    pos_nucleus = utiles.HMS2deg(ra=nucleus_RA2.replace(':', ' '), dec=nucleus_Dec2.replace(':', ' '))
    px_n, py_n = wcs_1.wcs_world2pix(float(pos_nucleus[0]), float(pos_nucleus[1]), 1)
    ax1.plot(px_n, py_n, 'x', color='#22CFAF') 
             
    # Plortting TH positions 
    # TH6 ~SSC IR-11 ~X-1
    pos_th6 = utiles.HMS2deg(ra=TH6_RA.replace(':', ' '), dec=TH6_Dec.replace(':', ' '))
    pxth6_n, pyth6_n = wcs_1.wcs_world2pix(float(pos_th6[0]), float(pos_th6[1]), 1)
    ax1.plot(pxth6_n, pyth6_n, '+', color='yellow')  
    # TH2 brightest
    pos_th2 = utiles.HMS2deg(ra=TH2_RA.replace(':', ' '), dec=TH2_Dec.replace(':', ' '))
    pxth2_n, pyth2_n = wcs_1.wcs_world2pix(float(pos_th2[0]), float(pos_th2[1]), 1)
    ax1.plot(pxth2_n, pyth2_n, '+', color='orange')  
    # TH7
    pos_th7 = utiles.HMS2deg(ra=TH7_RA.replace(':', ' '), dec=TH7_Dec.replace(':', ' '))
    pxth7_n, pyth7_n = wcs_1.wcs_world2pix(float(pos_th2[0]), float(pos_th2[1]), 1)
    ax1.plot(pxth7_n, pyth7_n, '+', color='orange') 
    # IR32
    pos_32 = utiles.HMS2deg(ra=IR32_RA.replace(':', ' '), dec=IR32_Dec.replace(':', ' '))
    px32_n, py32_n = wcs_1.wcs_world2pix(float(pos_32[0]), float(pos_32[1]), 1)
    ax1.plot(px32_n, py32_n, '+', color='lime') 
    
             
    # Bright blob and others from Watson1996
    pos_bb = utiles.HMS2deg(ra=bright_blob_RA.replace(':', ' '), dec=bright_blob_Dec.replace(':', ' '))
    pxbb_n, pybb_n = wcs_1.wcs_world2pix(float(pos_bb[0]), float(pos_bb[1]), 1)
    ax1.plot(pxbb_n, pybb_n, '+', color='green')
    
    pos_a = utiles.HMS2deg(ra=a_RA.replace(':', ' '), dec=a_Dec.replace(':', ' '))
    pxa_n, pya_n = wcs_1.wcs_world2pix(float(pos_a[0]), float(pos_a[1]), 1)
    ax1.plot(pxa_n, pya_n, '+', color='green') 
    
    pos_i = utiles.HMS2deg(ra=i_RA.replace(':', ' '), dec=i_Dec.replace(':', ' '))
    pxi_n, pyi_n = wcs_1.wcs_world2pix(float(pos_i[0]), float(pos_i[1]), 1)
    ax1.plot(pxi_n, pyi_n, '+', color='green') 
    
    pos_h = utiles.HMS2deg(ra=h_RA.replace(':', ' '), dec=h_Dec.replace(':', ' '))
    pxh_n, pyh_n = wcs_1.wcs_world2pix(float(pos_h[0]), float(pos_h[1]), 1)
    ax1.plot(pxh_n, pyh_n, '+', color='green') 
    
    #Alejandro
#    pos_a10 = utiles.HMS2deg(ra=a10_RA.replace(':', ' '), dec=a10_Dec.replace(':', ' '))
#    pxa10_n, pya10_n = wcs_1.wcs_world2pix(float(pos_a10[0]), float(pos_a10[1]), 1)
#    ax1.plot(pxa10_n, pya10_n, '+', color='cyan')
#    
#    pos_a11 = utiles.HMS2deg(ra=a11_RA.replace(':', ' '), dec=a11_Dec.replace(':', ' '))
#    pxa11_n, pya11_n = wcs_1.wcs_world2pix(float(pos_a11[0]), float(pos_a11[1]), 1)
#    ax1.plot(pxa11_n, pya11_n, '+', color='cyan')
#    
#    pos_a13 = utiles.HMS2deg(ra=a13_RA.replace(':', ' '), dec=a13_Dec.replace(':', ' '))
#    pxa13_n, pya13_n = wcs_1.wcs_world2pix(float(pos_a13[0]), float(pos_a13[1]), 1)
#    ax1.plot(pxa13_n, pya13_n, '+', color='cyan')
    
    # Fernandez Ontiveiros
    pos_o4 = utiles.HMS2deg(ra=O4_RA.replace(':', ' '), dec=O4_Dec.replace(':', ' '))
    pxo4_n, pyo4_n = wcs_1.wcs_world2pix(float(pos_o4[0]), float(pos_o4[1]), 1)
    ax1.plot(pxo4_n, pyo4_n, '+', color='pink')

    pos_o5 = utiles.HMS2deg(ra=O5_RA.replace(':', ' '), dec=O5_Dec.replace(':', ' '))
    pxo5_n, pyo5_n = wcs_1.wcs_world2pix(float(pos_o5[0]), float(pos_o5[1]), 1)
    ax1.plot(pxo5_n, pyo5_n, '+', color='pink')
    
    pos_o16 = utiles.HMS2deg(ra=O16_RA.replace(':', ' '), dec=O16_Dec.replace(':', ' '))
    pxo16_n, pyo16_n = wcs_1.wcs_world2pix(float(pos_o16[0]), float(pos_o16[1]), 1)
    ax1.plot(pxo16_n, pyo16_n, '+', color='pink')
    
    pos_o18 = utiles.HMS2deg(ra=O18_RA.replace(':', ' '), dec=O18_Dec.replace(':', ' '))
    pxo18_n, pyo18_n = wcs_1.wcs_world2pix(float(pos_o18[0]), float(pos_o18[1]), 1)
    ax1.plot(pxo18_n, pyo18_n, '+', color='pink')
    
    # My real positions
    pos_RA_yo = []
    pos_Dec_yo = []
    pos_yo = []
    dist_leroy = []
    from astropy.coordinates import SkyCoord 
    
    for i, line in hc3n_positions.iterrows():
        pos = utiles.HMS2deg(ra=line['RA_yo'].replace(':', ' '), dec=line['Dec_yo'].replace(':', ' '))
        px_n, py_n = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        ax1.plot(px_n, py_n, '.', color='k', markersize=3)      
        pos_RA_yo.append(float(pos[0]))
        pos_Dec_yo.append(float(pos[1]))
        
        pos_ra_s= line['RA_yo'].replace(':', 'h', 1).replace(':', 'm', 1) + 's'
        pos_Dec_s= line['Dec_yo'].replace(':', 'd', 1).replace(':', 'm', 1) + 's'
        pos_sky = SkyCoord(pos_ra_s, pos_Dec_s, frame='fk5', unit=(u.hourangle, u.deg), distance = 3.5*u.Mpc)
        pos_sky_leroy = SkyCoord(float(line['RA_leroy18']), float(line['Dec_leroy18']), frame='fk5', unit=(u.hourangle, u.deg), distance = 3.5*u.Mpc)
        
        
        pos_yo.append(pos)
        #Leroy positions
        #pos = utiles.HMS2deg(ra=line['RA_leroy18'].replace(':', ' '), dec=line['Dec_leroy18'].replace(':', ' '))
        px_n, py_n = wcs_1.wcs_world2pix(float(line['RA_leroy18']), float(line['Dec_leroy18']), 1)
        ax1.plot(px_n, py_n, 'x', color='lime', markersize=3)  
        
        ang_dist = utiles.ang_distance_btw_2points(float(pos[0]),float(pos[1]),float(line['RA_leroy18']),float(line['Dec_leroy18']))
        lin_dist = u_conversion.lin_size(D, ang_dist)
        dist_leroy.append(pos_sky_leroy.separation_3d(pos_sky).to(u.pc).value)
      
    # Alejandro
    pos_RA_ale = []
    pos_Dec_ale = []
    pos_ale=[]
    for i, line in H_positions.iterrows():
        if line['RA_H26'] != '-':
            pos = utiles.HMS2deg(ra=line['RA_H26'].replace(':', ' '), dec=line['Dec_H26'].replace(':', ' '))
            px_n, py_n = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
            ax1.plot(px_n, py_n, '+', color='magenta', markersize=3)   
            pos_RA_ale.append(float(pos[0]))
            pos_Dec_ale.append(float(pos[1]))
            print pos
            pos_ale.append(pos)
    
    
    # Diferencias con ALEJANDRO
    lin_dist_list = []
    for pyo,pale in zip(pos_yo,pos_ale):
        if '-' != pale[0] and '-' != pale[1]:
            dif_ra = float(pyo[0]) - float(pale[0])
            dif_dec = float(pyo[1]) - float(pale[1])
            ang_dist = utiles.ang_distance_btw_2points(float(pyo[0]),float(pyo[1]),float(pale[0]),float( pale[1]))
            lin_dist = u_conversion.lin_size(D, ang_dist)
            lin_dist_list.append(lin_dist.to(u.pc))
    print lin_dist_list
    
    
    
    ax1.coords.frame.set_color('w') # Frame color
    ax1.coords.frame.set_linewidth(0)
    
    
    
    fig.savefig(out_fig1_dir+'/ngc253_hc3nv7_'+cont_name+'_color_logspace_4_v'+version+'.eps', bbox_inches='tight', transparent=True, dpi=300)
    plt.close()
    
    
    # Plotting ALL CONTOURS
    ints1_2423 = 'MAD_CUB_II_ALL_NGC253_TE_219_220_0.19X0.29_briggs_v7.pbcor.fits'
    
    ngc253_hc3nv7_gr_fits_2423 = fits.open(ints1_2423)
    ngc253_hc3nv7_gr_data_2423 = ngc253_hc3nv7_gr_fits_2423[0]
    ngc253_hc3nv7_gr_shape_2423 = ngc253_hc3nv7_gr_data_2423.data.shape
    ngc253_hc3nv7_gr_header_2423 = ngc253_hc3nv7_gr_data_2423.header
    #ngc253_hc3nv7_gr_header['CDELT4'] = 1
    #ngc253_hc3nv7_gr_header['WCSAXES'] = 2
    wcs_hc3nv7_gr_2_2423 = WCS(ngc253_hc3nv7_gr_header_2423, naxis=2)
    
    fig = plt.figure()
    ax1 = fig.add_subplot((111), aspect='equal', projection=wcs_1)
    ax1.tick_params(direction='in')
    cm_max = np.nanmax(ngc253_hc3nv7_gr_data_2423.data[0,0,:,:])
    cm_min = np.abs(np.nanmin(ngc253_hc3nv7_gr_data_2423.data[0,0,:,:]))
    cm_std = np.nanstd(ngc253_hc3nv7_gr_data_2423.data[0,0,:,:])
    
    levels_cont218 = np.geomspace(2*cont_stdev, 1.5*cont_max, 5)
    levels_cont358 = np.geomspace(2*cont358_stdev, 1.5*cont358_max, 5)
    ax1.contour(ngc253_cont_data.data[0,0,:,:], colors='blue',  levels = levels_cont218, linewidths=0.5)
    ax1.contour(ngc253_cont358_data.data[0,0,:,:], colors='gray',  levels = levels_cont358, linewidths=0.5, transform=ax1.get_transform(wcs_1_358))

    
    if 'cont_218' in cont_fits_name:
        cont_name = 'cont_218'
        print cont_name

        plt.xlim([160, 483])#[120, 523]) #170, 474
        plt.ylim([175, 462])#[155, 482]) #336, 453
    elif 'cont_358' in cont_fits_name:
        cont_name = 'cont_358'
        print cont_name
        plt.xlim([150, 300])#[160, 483])#[120, 523]) #170, 474
        plt.ylim([160, 300])#[175, 462])#[155, 482]) #336, 453
    plt.ylabel('Dec (J2000)', labelpad=-1)
    plt.xlabel('RA (J2000)')
    
    pixsize_219 = 0.03 #arcsec
    ell_ax1_219 = 0.292 #arcsec
    ell_ax2_219 = 0.197 # arcsec
    pa_219 = 80.279
    
    pixsize_355 = 0.03 #arcsec es 0.05 pero estamos pintando sobre el continuo de 218 que tiene 0.03"
    ell_ax1_355 = 0.303 #arcsec
    ell_ax2_355 = 0.25 # arcsec
    pa_355 = -74.2

    c_219 = Ellipse(xy=(185, 195), width=ell_ax1_219/pixsize_219, height=ell_ax2_219/pixsize_219, angle=90-pa, edgecolor='orange', facecolor='orange', linewidth=0.5,
              transform=ax1.get_transform(wcs_1))
    c_355 = Ellipse(xy=(205, 195), width=ell_ax1_355/pixsize_355, height=ell_ax2_355/pixsize_355, angle=90-pa, edgecolor='r', facecolor='r', linewidth=0.5,
              transform=ax1.get_transform(wcs_1))
    ax1.add_patch(c_355)
    ax1.add_patch(c_219)
    
    ax1.coords[0].set_major_formatter('hh:mm:ss.ss')
    ax1.coords[0].set_ticks(size=6,color='k', exclude_overlapping=True, number = 5, width=1)#, spacing=0.2 * u.arcsec,
    ax1.coords[1].set_ticks(size=6,color='k', width=1)
    ax1.coords[0].set_separator((r'$^{\rm{h}}$', r'$^{\rm{m}}$', r'$^{\rm{s}}$'))
    
    ax1.coords[0].display_minor_ticks(True)
    ax1.coords[1].display_minor_ticks(True)
    
    stdev_all_ngc253 = []
    max_all_ngc253 = []
    min_all_ngc253 = []
    #for i, vel in enumerate(colores_1):
    for i in range(ngc253_hc3nv7_gr_shape[1]): 
        #if i > 4:
            std = np.nanstd(ngc253_hc3nv7_gr_data.data[0,i,:,:])
            stdev_all_ngc253_3938.append(std)
            maxx = np.nanmax(ngc253_hc3nv7_gr_data.data[0,i,:,:])
            max_all_ngc253_3938.append(maxx)
            minn = np.nanmin(ngc253_hc3nv7_gr_data.data[0,i,:,:])
            min_all_ngc253.append(minn)
            
            
            #ax1.contour(ngc253_hc3nv7_gr_data.data[0,i,:,:], colors='orange', levels = levels_ngc253, linewidths=0.5, transform=ax1.get_transform(wcs_hc3nv7_gr_2))
            # High signal-noise
            levels_ngc253_219 = np.geomspace(3.0*1.5/1000, 3*maxx, 4)
            ax1.contour(ngc253_hc3nv7_gr_data.data[0,i,:,:], colors='cyan', levels = levels_ngc253_219, linewidths=0.5, transform=ax1.get_transform(wcs_hc3nv7_gr_2))
            # Los signal-noise
            levels_ngc253_219_low = np.sort(np.concatenate((np.geomspace(3.5/1000, 45./1000, 10), np.array([4.5/1000, 4.0/1000, 5./1000, 6./1000, 10./1000]))))
            ax1.contour(ngc253_hc3nv7_gr_data.data[0,i,:,:], colors='orange', levels = levels_ngc253_219_low, linewidths=0.5, transform=ax1.get_transform(wcs_hc3nv7_gr_2))
            
    stdev_all_ngc253_3938 = []
    max_all_ngc253_3938 = []
    min_all_ngc253_3938 = []
    #for i, vel in enumerate(colores_1):
    for i in range(ngc253_hc3nv7_gr_shape_3938[1]): 
        #if i > 4:
            std = np.nanstd(ngc253_hc3nv7_gr_data_3938.data[0,i,:,:])
            stdev_all_ngc253_3938.append(std)
            maxx = np.nanmax(ngc253_hc3nv7_gr_data_3938.data[0,i,:,:])
            max_all_ngc253_3938.append(maxx)
            minn = np.nanmin(ngc253_hc3nv7_gr_data_3938.data[0,i,:,:])
            min_all_ngc253_3938.append(minn)
            levels_ngc253_3938 = np.geomspace(3.3*std, maxx, 3)
            ax1.contour(ngc253_hc3nv7_gr_data_3938.data[0,i,:,:], colors='red', levels = levels_ngc253_3938, linewidths=0.5, linestyles='solid', transform=ax1.get_transform(wcs_hc3nv7_gr_2_3938))
    

    for hotcore in hc3n_prop.itertuples():
        pos = utiles.HMS2deg(ra=hotcore.RA.replace('_', ' '), dec=hotcore.Dec.replace('_', ' '))
        px, py = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        #ax1.plot(px, py, '.r')
        #ax1.text(px, py+py/40., str(hotcore.Ind), fontsize=10, ha='center', va='center', color='white')
        px_m = 0
        if hotcore.Ind in [4, 7, 9, 11]:
            py_m = -20
        elif hotcore.Ind==3:
            py_m = 35
            px_m = 20
        elif hotcore.Ind_2 ==11:
            px_m = -10
            py_m = -20
        else:
            py_m = 20
        ax1.annotate(str(hotcore.Ind_ok), xy=(px,py), xytext=(px+px_m,py+py_m),
                arrowprops={'arrowstyle': '-', 'color': 'k'}, va='center', color = 'k')
        
    pos_6and7 = [['00 47 33.01', '-25 17 19.42'], ['00 47 33.01', '-25 17 19.02']]
    for p, posdeg in enumerate(pos_6and7):
        px_m = 0
        py_m = 20
        if p == 0:
            ind = 6
            px_m = -10
            py_m = -20
        else:
            ind = 7
        pos = utiles.HMS2deg(ra=posdeg[0], dec=posdeg[1])
        px, py = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        
        ax1.annotate(str(ind), xy=(px,py), xytext=(px+px_m,py+py_m),
                arrowprops={'arrowstyle': '-', 'color': 'k'}, va='center', color = 'k')
        
    # Plotting the kinematic center
    pos_nucleus = utiles.HMS2deg(ra=nucleus_RA2.replace(':', ' '), dec=nucleus_Dec2.replace(':', ' '))
    px_n, py_n = wcs_1.wcs_world2pix(float(pos_nucleus[0]), float(pos_nucleus[1]), 1)
    ax1.plot(px_n, py_n, 'x', color='#22CFAF') 
             
    # Plortting TH positions 
    # TH6 ~SSC IR-11 ~X-1
    pos_th6 = utiles.HMS2deg(ra=TH6_RA.replace(':', ' '), dec=TH6_Dec.replace(':', ' '))
    pxth6_n, pyth6_n = wcs_1.wcs_world2pix(float(pos_th6[0]), float(pos_th6[1]), 1)
    ax1.plot(pxth6_n, pyth6_n, '+', color='yellow')  
    # TH2 brightest
    pos_th2 = utiles.HMS2deg(ra=TH2_RA.replace(':', ' '), dec=TH2_Dec.replace(':', ' '))
    pxth2_n, pyth2_n = wcs_1.wcs_world2pix(float(pos_th2[0]), float(pos_th2[1]), 1)
    ax1.plot(pxth2_n, pyth2_n, '+', color='orange')  
    # TH7
    pos_th7 = utiles.HMS2deg(ra=TH7_RA.replace(':', ' '), dec=TH7_Dec.replace(':', ' '))
    pxth7_n, pyth7_n = wcs_1.wcs_world2pix(float(pos_th2[0]), float(pos_th2[1]), 1)
    ax1.plot(pxth7_n, pyth7_n, '+', color='orange') 
    # IR32
    pos_32 = utiles.HMS2deg(ra=IR32_RA.replace(':', ' '), dec=IR32_Dec.replace(':', ' '))
    px32_n, py32_n = wcs_1.wcs_world2pix(float(pos_32[0]), float(pos_32[1]), 1)
    ax1.plot(px32_n, py32_n, '+', color='lime') 
    
             
    # Bright blob and others from Watson1996
    pos_bb = utiles.HMS2deg(ra=bright_blob_RA.replace(':', ' '), dec=bright_blob_Dec.replace(':', ' '))
    pxbb_n, pybb_n = wcs_1.wcs_world2pix(float(pos_bb[0]), float(pos_bb[1]), 1)
    ax1.plot(pxbb_n, pybb_n, '+', color='green')
    
    pos_a = utiles.HMS2deg(ra=a_RA.replace(':', ' '), dec=a_Dec.replace(':', ' '))
    pxa_n, pya_n = wcs_1.wcs_world2pix(float(pos_a[0]), float(pos_a[1]), 1)
    ax1.plot(pxa_n, pya_n, '+', color='green') 
    
    pos_i = utiles.HMS2deg(ra=i_RA.replace(':', ' '), dec=i_Dec.replace(':', ' '))
    pxi_n, pyi_n = wcs_1.wcs_world2pix(float(pos_i[0]), float(pos_i[1]), 1)
    ax1.plot(pxi_n, pyi_n, '+', color='green') 
    
    pos_h = utiles.HMS2deg(ra=h_RA.replace(':', ' '), dec=h_Dec.replace(':', ' '))
    pxh_n, pyh_n = wcs_1.wcs_world2pix(float(pos_h[0]), float(pos_h[1]), 1)
    ax1.plot(pxh_n, pyh_n, '+', color='green') 

    
    # Fernandez Ontiveiros
    pos_o4 = utiles.HMS2deg(ra=O4_RA.replace(':', ' '), dec=O4_Dec.replace(':', ' '))
    pxo4_n, pyo4_n = wcs_1.wcs_world2pix(float(pos_o4[0]), float(pos_o4[1]), 1)
    ax1.plot(pxo4_n, pyo4_n, '+', color='pink')

    pos_o5 = utiles.HMS2deg(ra=O5_RA.replace(':', ' '), dec=O5_Dec.replace(':', ' '))
    pxo5_n, pyo5_n = wcs_1.wcs_world2pix(float(pos_o5[0]), float(pos_o5[1]), 1)
    ax1.plot(pxo5_n, pyo5_n, '+', color='pink')
    
    pos_o16 = utiles.HMS2deg(ra=O16_RA.replace(':', ' '), dec=O16_Dec.replace(':', ' '))
    pxo16_n, pyo16_n = wcs_1.wcs_world2pix(float(pos_o16[0]), float(pos_o16[1]), 1)
    ax1.plot(pxo16_n, pyo16_n, '+', color='pink')
    
    pos_o18 = utiles.HMS2deg(ra=O18_RA.replace(':', ' '), dec=O18_Dec.replace(':', ' '))
    pxo18_n, pyo18_n = wcs_1.wcs_world2pix(float(pos_o18[0]), float(pos_o18[1]), 1)
    ax1.plot(pxo18_n, pyo18_n, '+', color='pink')
    
    # My real positions
    pos_RA_yo = []
    pos_Dec_yo = []
    pos_yo = []
    dist_leroy = []
    for i, line in hc3n_positions.iterrows():
        pos = utiles.HMS2deg(ra=line['RA_yo'].replace(':', ' '), dec=line['Dec_yo'].replace(':', ' '))
        px_n, py_n = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        ax1.plot(px_n, py_n, '.', color='k', markersize=3)      
        pos_RA_yo.append(float(pos[0]))
        pos_Dec_yo.append(float(pos[1]))
        pos_yo.append(pos)
        #Leroy positions
        #pos = utiles.HMS2deg(ra=line['RA_leroy18'].replace(':', ' '), dec=line['Dec_leroy18'].replace(':', ' '))
        px_n, py_n = wcs_1.wcs_world2pix(float(line['RA_leroy18']), float(line['Dec_leroy18']), 1)
        ax1.plot(px_n, py_n, 'x', color='lime', markersize=3)  
        
        ang_dist = utiles.ang_distance_btw_2points(float(pos[0]),float(pos[1]),float(line['RA_leroy18']),float(line['Dec_leroy18']))
        lin_dist = u_conversion.lin_size(D, ang_dist)
        dist_leroy.append(lin_dist.to(u.pc))
      
    # Alejandro
    pos_RA_ale = []
    pos_Dec_ale = []
    pos_ale=[]
    for i, line in H_positions.iterrows():
        if line['RA_H26'] != '-':
            pos = utiles.HMS2deg(ra=line['RA_H26'].replace(':', ' '), dec=line['Dec_H26'].replace(':', ' '))
            px_n, py_n = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
            ax1.plot(px_n, py_n, '+', color='magenta', markersize=3)   
            pos_RA_ale.append(float(pos[0]))
            pos_Dec_ale.append(float(pos[1]))
            print pos
            pos_ale.append(pos)
    
    
    
    hc3n_path ='Hotcores_v4_all'+'/hc3n_obs_results_v0.txt'   
    hc3n_prop = pd.read_csv('/Users/frico/Documents/data/NGC253_H3O+/'+hc3n_path, delim_whitespace= True, header=0, comment='#')
                            
    L_protoSSC_LTE = [1.1E8, 5.2E8, 7.9E8, 5.3E8, 3.1E8, 0.0001E8,  0.0001E8, 1.8E8, 0.02E8, 0.07E8, 2.6E8,0.1E8, 22.7E8, 21.7E8]
    Tvib_LTE = [216., 304., 337., 326.,269.,92.,95.,217.,90.,132.,165.,140.,393.,312.]
    sizes_pc_yo = hc3n_prop['Source_size_pc_v7_219_g'].tolist()[:-1]
    sizes_pc_yo[5] = 1.7
    sizes_pc_yo[6] = 1.7
    s_si = _si.sigma_sb
    
    l9 = u_conversion.stef_boltz(sizes_pc_yo[8]*(1 * u.pc).to(u.m).value/2., 90).to(u.Lsun)/(10**8)
    L9 = (4*np.pi*s_si*(sizes_pc_yo[8]/2*(1 * u.pc).to(u.m))**2*(u.K*90)**4).to(u.Lsun)/(10**8)
    l10 = u_conversion.stef_boltz(sizes_pc_yo[9]*(1 * u.pc).to(u.m).value/2., 132).to(u.Lsun)/(10**8)
    L10 = (4*np.pi*s_si*(sizes_pc_yo[9]/2*(1 * u.pc).to(u.m))**2*(u.K*132)**4).to(u.Lsun)/10**8
    l11 = u_conversion.stef_boltz(sizes_pc_yo[10]*(1 * u.pc).to(u.m).value/2., 265).to(u.Lsun)/(10**8)
    L11 = (4*np.pi*s_si*(sizes_pc_yo[10]/2*(1 * u.pc).to(u.m))**2*(u.K*165)**4).to(u.Lsun)/10**8
    l12 = u_conversion.stef_boltz(sizes_pc_yo[11]*(1 * u.pc).to(u.m).value/2., 140).to(u.Lsun)/(10**8)
    L12 = (4*np.pi*s_si*(sizes_pc_yo[11]/2*(1 * u.pc).to(u.m))**2*(u.K*140)**4).to(u.Lsun)/10**8
    
    L7 = 4*np.pi*s_si*(sizes_pc_yo[6]*(1 * u.pc).to(u.m)**2)*(u.K*Tvib_LTE[6])**4
    
    NHC3N_yo = [10**15.3, 10**15.5, 10**15.1, 10**15.2, 10**15.6, 10**14.6, 10**14.7, 10**15.6, 10**14.7, 10**15.4, 10**15.4, 10**15.8, 10**15.5, 10**16.2]
    leroy_sizes_pc = [2.7, 1.2, 2.6, 2.5, 2.1, 2.1, 2.9, 1.9, 2.6, 3.5, 2.9, 4.3, 1.6, 1.6]
    leroy_gas_mass = [10**4.9, 10**4.7, 10**5.1, 10**5.1, 10**5.3, 10**3.6, 10**4.5, 10**5.2, 10**4.7, 10**5.2, 10**4.5, 10**4.1, 10**5.2, 10**5.7]
    leroy_star_mass = [10**4.3, 10**4.3, 10**4.1, 10**5.0, 10**5.4, 10**5.3, 10**4.5, 10**4.8, 10**5.5, 10**5.3, 10**5.6, 10**6.0, 10**4.8, 10**5.5]
    l_ir_nuc = 1.8E10
    lum_to_mass = 1000.
    hydrogen_mass_msun = 2.*(_si.m_p + _si.m_e).to(u.Msun).value
    
    nh2_leroy_list = []
    Nh2_leroy_list = []
    Xhc3n_list = []
    L_ir_leroy = 0
    L_ir_leroy_list = []
    T_cond = []
    T_cond3 = []
    Tcond_2 = []
    L_msSSC_L_protoSSC_ratio = []
    for i, mass in enumerate(leroy_star_mass):
        L_ir_leroy = L_ir_leroy + (lum_to_mass*mass)
        L_ir_leroy_list.append(lum_to_mass*mass)
        
        vol = (4./3)*np.pi*((leroy_sizes_pc[i]/2.)*(1 * u.pc).to(u.cm).value)**3
        
        nh2_leroy_list.append(leroy_gas_mass[i]/(vol*hydrogen_mass_msun))
        Nh2_leroy_list.append((leroy_sizes_pc[i]/2.)*(1 * u.pc).to(u.cm).value*nh2_leroy_list[i])
        Xhc3n_list.append(NHC3N_yo[i]/Nh2_leroy_list[i])
        
        L_msSSC_L_protoSSC_ratio.append(L_ir_leroy_list[i]/L_protoSSC_LTE[i])
        T_cond.append(Tvib_LTE[i]*np.sqrt(sizes_pc_yo[i]/leroy_sizes_pc[i]))
        T_cond3.append(Tvib_LTE[i]*np.sqrt(sizes_pc_yo[i]/leroy_sizes_pc[i]))

        lsum = ((L_ir_leroy_list[i]+L_protoSSC_LTE[i])*u.Lsun).to(u.W)
        rad_m = (leroy_sizes_pc[i] * u.pc).to(u.m)/2.
        Tcond_2.append((lsum/(4.*np.pi*s_si*(rad_m**2)))**(1./4))
                        
                
        
        
    #T = u.Quantity(T, u.K)
    #r = u.Quantity(r, u.m)
    #L = 4*pi*(r**2)*s_si*T**4
        
    L_ir_percen_leroy = 100.*L_ir_leroy/l_ir_nuc
    
    
    
    #ax1.coords.frame.set_color('w') # Frame color
    #ax1.coords.frame.set_linewidth(0)
    
    
    
    fig.savefig(out_fig1_dir+'/ngc253_ALLcontonly_4_v'+version+'_b.eps', bbox_inches='tight', transparent=True, dpi=300)
    plt.close()
    

# CS J=6-5 293.912 GHz
ints_cs_a = 'MAD_CUB_II_NGC253_TE_294_296_natural_CS_20KMS.pbcor_BL.fits'
intsall_cs_a=  'MAD_CUB_II_NGC253_TE_294_296_natural_CS_1plane.pbcor_BL.fits'
# CS J=7-6 342.883 GHz
ints_cs_b = 'MAD_CUB_II_NGC253_TE_342_344_0.31X0.25_briggs_v1_CS_20KMS.pbcor.fits'
intsall_cs_b=  'MAD_CUB_II_NGC253_TE_342_344_0.31X0.25_briggs_v1_CS_1plane.pbcor.fits'

# Sacar espectros

ngc253_cs_a_fits = fits.open(ints_cs_b)
ngc253_cs_a_data = ngc253_cs_a_fits[0]
ngc253_cs_a_shape = ngc253_cs_a_data.data.shape
ngc253_cs_a_header = ngc253_cs_a_data.header
#ngc253_hc3nv7_gr_header['CDELT3'] = 1
ngc253_cs_a_header['CDELT4'] = 1
wcs_cs_a_2 = WCS(ngc253_cs_a_header)
wcs_cs_a_2.wcs.ctype = ['RA---SIN', 'DEC--SIN', 'FREQ', 'STOKES']
wcs_cs_a_2 = wcs_cs_a_2.dropaxis(3)
wcs_cs_a_2 = wcs_cs_a_2.dropaxis(2)

## (flux, rms_noise, freq_axis, sigma_clip_threshold=1.8)
#            a = scont.cont_finding.c_sigmaclip(ngc253_cs_a_data[0,i,:,:], 0.002, 0)

#import sys
#sys.argv = ['--iname MAD_CUB_II_NGC253_TE_342_344_0.31X0.25_briggs_v1_CS_20KMS', '--noise 0.002','--continuum']
#scont.main()

# Plotting CS over cont
if plot_CS_over_cont == True:
    fig = plt.figure()
    ax1 = fig.add_subplot((111), aspect='equal', projection=wcs_1)
    ax1.tick_params(direction='in')
    cm_max = np.nanmax(ngc253_cont_data.data[0,0,:,:])
    cm_min = np.abs(np.nanmin(ngc253_cont_data.data[0,0,:,:]))
    cm_std = np.nanstd(ngc253_cont_data.data[0,0,:,:])
    ax1.imshow(ngc253_cont_data.data[0,0,:,:], origin='lower', vmax=cont_max, vmin=0.2*cont_stdev, cmap=cm.jet, interpolation="none")# vmax=cm_max, vmin=cm_std*2, interpolation="none")
    plt.xlim([120, 523]) #170, 474
    plt.ylim([155, 482]) #336, 453
    plt.ylabel('Dec (J2000)')
    plt.xlabel('RA (J2000)')
    
    ax1.coords[0].set_major_formatter('hh:mm:ss.ss')
    ax1.coords[0].set_ticks(size=8,color='k', exclude_overlapping=True, number = 5, width=2)#, spacing=0.2 * u.arcsec,
    ax1.coords[1].set_ticks(size=8,color='k', width=2)
    ax1.coords[0].set_separator((r'$^{\rm{h}}$', r'$^{\rm{m}}$', r'$^{\rm{s}}$'))
    
    
    
    stdev_all_ngc253 = []
    max_all_ngc253 = []
    min_all_ngc253 = []
    #for i, vel in enumerate(colores_1):
    for i in range(ngc253_cs_a_shape[1]): 
        #if i > 4:
            std = np.nanstd(ngc253_cs_a_data.data[0,i,:,:])
            stdev_all_ngc253.append(std)
            maxx = np.nanmax(ngc253_cs_a_data.data[0,i,:,:])
            max_all_ngc253.append(maxx)
            minn = np.nanmin(ngc253_cs_a_data.data[0,i,:,:])
            min_all_ngc253.append(minn)
            levels_ngc253 = np.linspace(4.*std, maxx, 3)
            
            
            ax1.contour(ngc253_cs_a_data.data[0,i,:,:], colors='white', levels = levels_ngc253, linewidths=0.7, transform=ax1.get_transform(wcs_cs_a_2))
    
    for hotcore in hc3n_prop.itertuples():
        pos = utiles.HMS2deg(ra=hotcore.RA.replace('_', ' '), dec=hotcore.Dec.replace('_', ' '))
        px, py = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        #ax1.plot(px, py, '.r')
        #ax1.text(px, py+py/40., str(hotcore.Ind), fontsize=10, ha='center', va='center', color='white')
        px_m = 0
        if hotcore.Ind in [4, 7, 9, 11]:
            py_m = -20
        elif hotcore.Ind==3:
            py_m = 35
            px_m = 20
        elif hotcore.Ind_2 ==11:
            px_m = -10
            py_m = -20
        else:
            py_m = 20
        ax1.annotate(str(hotcore.Ind_ok), xy=(px,py), xytext=(px+px_m,py+py_m),
                arrowprops={'arrowstyle': '-', 'color': 'w'}, va='center', color = 'white')
        
    pos_6and7 = [['00 47 33.01', '-25 17 19.42'], ['00 47 33.01', '-25 17 19.02']]
    for p, posdeg in enumerate(pos_6and7):
        px_m = 0
        py_m = 20
        if p == 0:
            ind = 6
            px_m = -10
            py_m = -20
        else:
            ind = 7
        pos = utiles.HMS2deg(ra=posdeg[0], dec=posdeg[1])
        px, py = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        
        ax1.annotate(str(ind), xy=(px,py), xytext=(px+px_m,py+py_m),
                arrowprops={'arrowstyle': '-', 'color': 'w'}, va='center', color = 'white')
        
        
    
    #ax1.xaxis.set_tick_params(width=20)
    #ax1.yaxis.set_tick_params(width=20)
    fig.savefig(out_fig1_dir+'/ngc253_CS_b_v'+version+'.eps', bbox_inches='tight', transparent=True, dpi=300)
    plt.close()



# Plotting HC3Nv0 over cont
if plot_v0_over_cont == True:
    fig = plt.figure()
    ax1 = fig.add_subplot((111), aspect='equal', projection=wcs_1)
    ax1.tick_params(direction='in')
    cm_max = np.nanmax(ngc253_cont_data.data[0,0,:,:])
    cm_min = np.abs(np.nanmin(ngc253_cont_data.data[0,0,:,:]))
    cm_std = np.nanstd(ngc253_cont_data.data[0,0,:,:])
    ax1.imshow(ngc253_cont_data.data[0,0,:,:], origin='lower', vmax=cont_max, vmin=0.2*cont_stdev, cmap=cm.jet, interpolation="none")# vmax=cm_max, vmin=cm_std*2, interpolation="none")
    plt.xlim([120, 523]) #170, 474
    plt.ylim([155, 482]) #336, 453
    
    plt.ylabel('Dec (J2000)')
    plt.xlabel('RA (J2000)')
    
    ax1.coords[0].set_major_formatter('hh:mm:ss.ss')
    ax1.coords[0].set_ticks(size=4,color='k', exclude_overlapping=True, number = 5)#, spacing=0.2 * u.arcsec,
    ax1.coords[0].set_separator((r'$^{\rm{h}}$', r'$^{\rm{m}}$', r'$^{\rm{s}}$'))
    
    
    
    stdev_all_ngc253 = []
    max_all_ngc253 = []
    min_all_ngc253 = []
    for i in range(ngc253_hc3nv0_gr_shape[1]): 
        std = np.nanstd(ngc253_hc3nv0_gr_data.data[0,i,:,:])
        stdev_all_ngc253.append(std)
        maxx = np.nanmax(ngc253_hc3nv0_gr_data.data[0,i,:,:])
        max_all_ngc253.append(maxx)
        minn = np.nanmin(ngc253_hc3nv0_gr_data.data[0,i,:,:])
        min_all_ngc253.append(minn)
        levels_ngc253 = np.linspace(3.2*std, maxx, 3)
        ax1.contour(ngc253_hc3nv0_gr_data.data[0,i,:,:], colors='white', levels = levels_ngc253, linewidths=0.7, transform=ax1.get_transform(wcs_hc3nv0_gr_2))
    
    for hotcore in hc3n_prop.itertuples():
        pos = utiles.HMS2deg(ra=hotcore.RA.replace('_', ' '), dec=hotcore.Dec.replace('_', ' '))
        px, py = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        #ax1.plot(px, py, '.r')
        #ax1.text(px, py+py/40., str(hotcore.Ind), fontsize=10, ha='center', va='center', color='white')
        if hotcore.Ind in [4, 7, 9, 11]:
            py_m = -20
        else:
            py_m = 20
        ax1.annotate(str(hotcore.Ind_ok), xy=(px,py), xytext=(px,py+py_m),
                arrowprops={'arrowstyle': '-', 'color': 'k'}, va='center', color = 'white')
    
    fig.savefig(out_fig1_dir+'/ngc253_hc3nv0_v'+version+'.eps', bbox_inches='tight', transparent=True, dpi=300)
    plt.close()
    
    
    
    

# Plotting HC3Nv7=1 contours only
if plot_v7_contours == True:
    fig = plt.figure()
    ax1 = fig.add_subplot((111), aspect='equal', projection=wcs_hc3nv7_gr_2)
    ax1.tick_params(direction='in')
    plt.xlim([120, 523]) #170, 474
    plt.ylim([155, 482]) #336, 453
    plt.ylabel('Dec (J2000)')
    plt.xlabel('RA (J2000)')
    
    ax1.coords[0].set_major_formatter('hh:mm:ss.ss')
    ax1.coords[0].set_ticks(size=4,color='k', exclude_overlapping=True, number = 5)#, spacing=0.2 * u.arcsec,
    ax1.coords[0].set_separator((r'$^{\rm{h}}$', r'$^{\rm{m}}$', r'$^{\rm{s}}$'))
    
    
    
    stdev_all_ngc253 = []
    max_all_ngc253 = []
    min_all_ngc253 = []
    #for i, vel in enumerate(colores_1):
    for i in range(ngc253_hc3nv7_gr_shape[1]): 
        #if i > 4:
            std = np.nanstd(ngc253_hc3nv7_gr_data.data[0,i,:,:])
            stdev_all_ngc253.append(std)
            maxx = np.nanmax(ngc253_hc3nv7_gr_data.data[0,i,:,:])
            max_all_ngc253.append(maxx)
            minn = np.nanmin(ngc253_hc3nv7_gr_data.data[0,i,:,:])
            min_all_ngc253.append(minn)
            levels_ngc253 = np.linspace(3.2*std, maxx, 3)
            ax1.contour(ngc253_hc3nv7_gr_data.data[0,i,:,:], colors='white', levels = levels_ngc253, linewidths=0.7)
    
    #for hotcore in hc3n_prop.itertuples():
    #    pos = utiles.HMS2deg(ra=hotcore.RA.replace('_', ' '), dec=hotcore.Dec.replace('_', ' '))
    #    px, py = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
    #    #ax1.plot(px, py, '.r')
    #    #ax1.text(px, py+py/40., str(hotcore.Ind), fontsize=10, ha='center', va='center', color='white')
    #    if hotcore.Ind in [4, 7, 11]:
    #        py_m = -20
    #    else:
    #        py_m = 20
    #    ax1.annotate(str(hotcore.Ind_ok), xy=(px,py), xytext=(px,py+py_m),
    #            arrowprops={'arrowstyle': '-', 'color': 'k'}, va='center', color = 'white')
    
    fig.savefig(out_fig2_dir+'/ngc253_hc3nv7_only_v'+version+'.png', bbox_inches='tight', transparent=True, dpi=300)
    plt.close()


# Plotting HC3Nv0 over cont
if plot_v0_contours == True:
    fig = plt.figure()
    ax1 = fig.add_subplot((111), aspect='equal', projection=wcs_hc3nv0_gr_2)
    ax1.tick_params(direction='in')
    plt.xlim([120, 523]) #170, 474
    plt.ylim([155, 482]) #336, 453
    
    plt.ylabel('Dec (J2000)')
    plt.xlabel('RA (J2000)')
    
    ax1.coords[0].set_major_formatter('hh:mm:ss.ss')
    ax1.coords[0].set_ticks(size=4,color='k', exclude_overlapping=True, number = 5)#, spacing=0.2 * u.arcsec,
    ax1.coords[0].set_separator((r'$^{\rm{h}}$', r'$^{\rm{m}}$', r'$^{\rm{s}}$'))

    stdev_all_ngc253 = []
    max_all_ngc253 = []
    min_all_ngc253 = []
    for i in range(ngc253_hc3nv0_gr_shape[1]): 
        std = np.nanstd(ngc253_hc3nv0_gr_data.data[0,i,:,:])
        stdev_all_ngc253.append(std)
        maxx = np.nanmax(ngc253_hc3nv0_gr_data.data[0,i,:,:])
        max_all_ngc253.append(maxx)
        minn = np.nanmin(ngc253_hc3nv0_gr_data.data[0,i,:,:])
        min_all_ngc253.append(minn)
        levels_ngc253 = np.linspace(3.2*std, maxx, 3)
        ax1.contour(ngc253_hc3nv0_gr_data.data[0,i,:,:], colors='white', levels = levels_ngc253, linewidths=0.7)
    
    for hotcore in hc3n_prop.itertuples():
        pos = utiles.HMS2deg(ra=hotcore.RA.replace('_', ' '), dec=hotcore.Dec.replace('_', ' '))
        px, py = wcs_1.wcs_world2pix(float(pos[0]), float(pos[1]), 1)
        #ax1.plot(px, py, '.r')
        #ax1.text(px, py+py/40., str(hotcore.Ind), fontsize=10, ha='center', va='center', color='white')
        if hotcore.Ind in [4, 7, 11]:
            py_m = -20
        else:
            py_m = 20
        ax1.annotate(str(hotcore.Ind), xy=(px,py), xytext=(px,py+py_m),
                arrowprops={'arrowstyle': '-', 'color': 'k'}, va='center', color = 'white')
    
    fig.savefig(out_fig2_dir+'/ngc253_hc3nv0_only_v'+version+'.eps', bbox_inches='tight', transparent=True, dpi=300)
    plt.close()

#==============================================================================
# MADCUBA NGC253 (ALMA) RotDiagram HC3N v7=1
#==============================================================================

slopes = []
intercepts = []
if plot_rot_diag == True:
    os.chdir(dworkdir_spec)  
    hot_cores = range(len(hc3n_prop))
    ngc253_rotdiag = []
    for j, hc in enumerate(hc3n_prop['Ind_ok']):
        path = 'HC_'+str(hc)+'_HR_v2'
        rotdiag = pd.read_csv(path+'/HC_'+str(j+1)+'_rotdiag_all', delim_whitespace= True, header=None)
        rotdiag.columns = ['E', 'N']
        ngc253_rotdiag.append(rotdiag)
        fig = plt.figure()
        ax1 = fig.add_subplot((111))
        #plt.scatter(rotdiag['E'], rotdiag['N'], marker='.', color='k')
        m, b = np.polyfit(rotdiag['E'], rotdiag['N'], 1)
        slopes.append(m)
        intercepts.append(b)
        fit = np.polyfit(rotdiag['E'], rotdiag['N'], 1)
        fit_fn = np.poly1d(fit) 
        xx = np.arange(0, np.max(rotdiag['E'])+100)
        ax1.plot(rotdiag['E'], rotdiag['N'], 'k.', xx, fit_fn(xx), '-r', linewidth=0.75)
        ax1.xaxis.set_tick_params(top ='on', labeltop='off')
        ax1.yaxis.set_tick_params(right='on', labelright='off')
        plt.xlim(0, np.max(xx))
        plt.ylim(np.min(fit_fn(xx)), np.max(fit_fn(xx)))
        if hc == '14':
            plt.ylim(11.0, np.max(fit_fn(xx)))
            #ax5.text(local_max_x[h], hc_max+hc_max*0.1, max_order[h], ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=5)
            ax1.text(0.6, 0.9,r'$\rm{T}_{ex}=320.16 \rm{K} \quad \log{N}=16.52$', ha='left', va='center', transform=ax1.transAxes, color='k', fontsize=8)

            
        plt.ylabel(r'log (N$_{\rm u}$/g$_{\rm u}$)')
        plt.xlabel(r'E$_{\rm u}$/k (K)')
        fig.savefig(out_fig3_dir+'/ngc253_rotgiag_HC3Nv7_HC'+str(hc)+'_v'+version+'.png', bbox_inches='tight', transparent=True, dpi=600)
        plt.close()
    

#line_hc3nv7 = np.arange(0, 8+1)
#for i, hc in enumerate(hot_cores):
os.chdir(dworkdir_spec)   
#redshift_scale.velocities_fromfreq(freq_obs, rest_freq, z_obs)

if plot_HC_spec == True:  
    hot_cores = range(len(hc3n_prop))
    local_max_y_list = []
    local_max_x_list = []
    for ii,i in enumerate(hc3n_prop['Ind_ok']):
        path = 'HC_'+str(i)+'_HR_v2'
        local_max_y_list.append([])
        local_max_x_list.append([])
        for j in range(4):
            spec = pd.read_csv(path+'/'+str(j+1)+'_MADCUBA_1-1_0-0_0-0_v4', delim_whitespace= True, header=None)
            spec.columns = ['vel', 'int', 'fit']
            
            
            a = np.array(spec['fit'].tolist())
            local_max_bool = np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True]
            local_max_y = spec['fit'][local_max_bool].tolist()
            local_max_x = spec['vel'][local_max_bool].tolist()
            local_max_y_list[ii].append(local_max_y)
            local_max_x_list[ii].append(local_max_x)
            plt.figure()
            plt.plot(spec['vel'], spec['int'], linewidth=0.8, drawstyle='steps-mid', color='k')
            plt.plot(spec['vel'], spec['fit'], color='red', linewidth=0.8)
            
            if j == 3:
                lines = ['HC3Nv0=v7_218.3']
                max_order = ['HC3N,v=0']
            elif j==2:
                lines = ['HC3Nv7=2_219.7_&_HC3Nv7=1_219.5_&_HC3Nv6=v7_219.4']
                max_order = ['HC3N,v7=2','HC3N,v7=2', 'HC3N,v7=2', 'HC3N,v6=v7=1','HC3N,v6=v7=1','HC3N,v6=v7=1','HC3N,v6=v7=1']
            elif j==1:
                lines = ['HC3Nv7=1_218.8_&_HC3Nv6=1_218.6']
                max_order = ['HC3N,v7=1','HC3N,v=6=1','HC3N,v6=1']
            elif j==0:
                lines = ['HC3Nv6=v7_219.4_&_HC3Nv7=1_219.1']
                max_order = ['HC3N,v6=v7=1','HC3N,v6=v7=1','HC3N,v7=1']
                
            plt.xlabel(r'V$_{\rm rad}$ (km/s)')
            plt.ylabel(r'Jy')
            plt.savefig(out_fig4_dir+'/ngc253_espec_HC'+str(i)+'_'+str(j+1)+'_v'+version+'.png', bbox_inches='tight', transparent=True, dpi=600)
            plt.close()
            
            
            
# =============================================================================
# Spec lines
# =============================================================================

if plot_panel_fig2 == True:
    hot_cores = range(len(hc3n_prop))
    local_max_y_list = []
    local_max_x_list = []
    hc3nv7_max = []
    for ii,i in enumerate(hc3n_prop['Ind_ok']):
        local_max_y_list.append([])
        local_max_x_list.append([])
        
        path = 'HC_'+str(i)+'_HR_v2'
        fig = plt.figure(figsize=(20, 5))
        gs1 = gridspec.GridSpec(1, 4)#, width_ratios=[2, 4 , 4, 4], height_ratios=[1])
        gs1.update(wspace = 0.3, top=0.75, bottom = 0.05)
        
        ax2 = fig.add_subplot(gs1[0])
        ax2.set_ylabel(r'Jy')
        ax3 = fig.add_subplot(gs1[1], sharex=ax2)
        ax4 = fig.add_subplot(gs1[2], sharex=ax2)
        ax5 = fig.add_subplot(gs1[3], sharex=ax2)
        hc_max_list = []
        hc_min_list = []
        for j in range(4):
            spec = pd.read_csv(path+'/'+str(j+1)+'_MADCUBA_1-1_0-0_0-0_v4', delim_whitespace= True, header=None)
            spec.columns = ['vel', 'int', 'fit']
            
            hc_max_list.append(np.nanmax(spec['int'].tolist()))
            hc_min_list.append(np.nanmin(spec['int'].tolist()))
            
        hc_max = np.max(hc_max_list)
        hc_min = np.min(hc_min_list)
        for j in range(4):
            spec = pd.read_csv(path+'/'+str(j+1)+'_MADCUBA_1-1_0-0_0-0_v4', delim_whitespace= True, header=None)
            spec.columns = ['vel', 'int', 'fit']
            
            
            a = np.array(spec['fit'].tolist())
            local_max_bool = np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True]
            local_max_y = spec['fit'][local_max_bool].tolist()
            local_max_x = spec['vel'][local_max_bool].tolist()
            
            local_max_x_list[ii].append(local_max_x)
            local_max_y_list[ii].append(local_max_y)
            
            
            if j == 3:
                lines = ['HC3Nv0=1_218.3']
                max_order = ['HC3N,v=0']
                ax2.plot(spec['vel'], spec['int'], linewidth=0.8, drawstyle='steps-mid', color='k')
                ax2.plot(spec['vel'], spec['fit'], color='red', linewidth=0.8)
                
                for h, vel in enumerate(local_max_x):
                    ax2.axvline(local_max_x[h], color='k', linestyle='--', lw=0.5)
                    ax2.text(local_max_x[h], hc_max+hc_max*0.2, max_order[h], ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=5)
            
            elif j==2:
                lines = ['HC3Nv7=2_219.7_&_HC3Nv7=1_219.5_&_HC3Nv6=v7_219.4']
                max_order = ['HC3N,v7=2','HC3N,v7=2', 'HC3N,v7=2', 'HC3N,v6=v7=1','HC3N,v6=v7=1','HC3N,v6=v7=1','HC3N,v6=v7=1']
                ax5.plot(spec['vel'], spec['int'], linewidth=0.8, drawstyle='steps-mid', color='k')
                ax5.plot(spec['vel'], spec['fit'], color='red', linewidth=0.8)
                if ii==(3-1):
                    local_max_x = [181.247, 228.133, 274.612, 487.308, 520.942, 599.394, 633.900]
                elif ii==(4-1) or ii==(5-1) or ii==(6-1) or ii==(7-1):
                    local_max_x = [181.247+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   228.133+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   274.612+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   487.308+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   520.942+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   599.394+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   633.900+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2]]
                elif ii==(1-1): 
                    local_max_x = [181.247+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   228.133+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   274.612+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   487.308+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   520.942+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   599.394+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   633.900+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2]]
                
                for h, vel in enumerate(local_max_x): 
                    ax5.axvline(local_max_x[h], color='k', linestyle='--', lw=0.5)
                    ax5.text(local_max_x[h], hc_max+hc_max*0.1, max_order[h], ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=5)
            
            elif j==1:
                lines = ['HC3Nv7=1_218.8_&_HC3Nv6=1_218.6']
                max_order = ['HC3N,v7=1','HC3N,v6=1','HC3N,v6=1']
                ax3.plot(spec['vel'], spec['int'], linewidth=0.8, drawstyle='steps-mid', color='k')
                ax3.plot(spec['vel'], spec['fit'], color='red', linewidth=0.8)
                if ii==(3-1):
                    local_max_x = [275.020, 283.767, 518.957]
                elif ii==(4-1) or ii==(5-1) or ii==(6-1) or ii==(7-1):
                    local_max_x = [275.020+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   283.767+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   518.957+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2]]
                    
                else:
                    local_max_x = [275.020+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   283.767+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   518.957+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2]]
                
                for h, vel in enumerate(local_max_x): 
                    ax3.axvline(local_max_x[h], color='k', linestyle='--', lw=0.5)
                    if h < 2:
                        ax3.text(local_max_x[h]-local_max_x[0]*0.12, hc_max+hc_max*0.1, max_order[0], ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=5)
                        ax3.text(local_max_x[h]+local_max_x[0]*0.1, hc_max+hc_max*0.1, max_order[1], ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=5)
                    else:
                        ax3.text(local_max_x[h], hc_max+hc_max*0.1, max_order[h], ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=5)
                        ax3.axvline(local_max_x[h], color='k', linestyle='--', lw=0.5)
            elif j==0:
                lines = ['HC3Nv6=v7_219.4_&_HC3Nv7=1_219.1']
                max_order = ['HC3N,v6=v7=1','HC3N,v6=v7=1','HC3N,v7=1']
                ax4.plot(spec['vel'], spec['int'], linewidth=0.8, drawstyle='steps-mid', color='k')
                ax4.plot(spec['vel'], spec['fit'], color='red', linewidth=0.8)
                hc3nv7_max.append(local_max_y[-1])
                if ii==(3-1):
                    local_max_x = [-84.86, -50.437, 276.305]
                    
                elif ii==(4-1) or ii==(5-1) or ii==(6-1) or ii==(7-1):
                    local_max_x = [-84.86+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   -50.437+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   276.305+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2]]
                elif ii==(1-1): 
                    local_max_xx = deepcopy(local_max_x)
                    local_max_x = [-84.86+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   -50.437+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   local_max_xx[-1]]
                
                for h, vel in enumerate(local_max_x): 
                    ax4.axvline(local_max_x[h], color='k', linestyle='--', lw=0.5)
                    ax4.text(local_max_x[h], hc_max+hc_max*0.1, max_order[h], ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=5)
            
                    
        axis_spec = [ax2, ax3, ax4, ax5]
        for ax in axis_spec:
            ax.tick_params(direction='in')
            ax.xaxis.set_tick_params(top ='on')
            ax.yaxis.set_tick_params(right='on', labelright='off')
            #ax.set_xlim([0, 700])
            ax.set_ylim([hc_min, hc_max+hc_max*0.4])
            ax.set_xlabel(r'v$_{rad}$ (km/s)')
                
            
        
        fig.savefig(out_fig5_dir+'/ngc253_grid_HC'+str(i)+'_v'+version+'.png', bbox_inches='tight', transparent=True, dpi=400)
        plt.close()

if plot_panel_fig == True:
    hot_cores = range(len(hc3n_prop))
    local_max_y_list = []
    local_max_x_list = []
    hc3nv7_max = []
    for ii,i in enumerate(hc3n_prop['Ind_ok']):
        local_max_y_list.append([])
        local_max_x_list.append([])
        
        path = 'HC_'+str(i)+'_HR_v2'
        fig = plt.figure(figsize=(11.0, 10))
        
        gs1 = gridspec.GridSpec(2, 2)#, width_ratios=[2, 4 , 4, 4], height_ratios=[1])
        gs1.update(wspace = 0.0, hspace=0.0, top=0.95, bottom = 0.05)
        #fig.set_tight_layout({'rect': [0, 0, 1, 0.95], 'pad': 0.3, 'h_pad': 0.3})
        
        ax2 = fig.add_subplot(gs1[0])
        
        ax3 = fig.add_subplot(gs1[1])#, sharex=ax2)
        ax4 = fig.add_subplot(gs1[2])#, sharex=ax2)
        ax5 = fig.add_subplot(gs1[3])#, sharex=ax2)
        hc_max_list = []
        hc_min_list = []
        for j in range(4):
            spec = pd.read_csv(path+'/'+str(j+1)+'_MADCUBA_1-1_0-0_0-0_v4', delim_whitespace= True, header=None)
            spec.columns = ['vel', 'int', 'fit']
            
            hc_max_list.append(np.nanmax(spec['int'].tolist()))
            hc_min_list.append(np.nanmin(spec['int'].tolist()))
            
        hc_max = np.max(hc_max_list)
        hc_min = np.min(hc_min_list)
        
        for j in range(4):
            spec = pd.read_csv(path+'/'+str(j+1)+'_MADCUBA_1-1_0-0_0-0_v4', delim_whitespace= True, header=None)
            spec.columns = ['vel', 'int', 'fit']
            
            
            a = np.array(spec['fit'].tolist())
            local_max_bool = np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True]
            local_max_y = spec['fit'][local_max_bool].tolist()
            local_max_x = spec['vel'][local_max_bool].tolist()
            
            local_max_x_list[ii].append(local_max_x)
            local_max_y_list[ii].append(local_max_y)
            
            
            if j == 3:
                lines = ['HC3Nv0=1_218.3']
                max_order = ['HC3N,v=0']
                ax2.plot(spec['vel'], spec['int'], linewidth=0.8, drawstyle='steps-mid', color='k')
                ax2.plot(spec['vel'], spec['fit'], color='red', linewidth=0.8)
                
                for h, vel in enumerate(local_max_x):
                    ax2.axvline(local_max_x[h], color='k', linestyle='--', lw=0.5)
                    ax2.text(local_max_x[h], hc_max+hc_max*0.2, max_order[h], ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=8)
            
            elif j==2:
                lines = ['HC3Nv7=2_219.7_&_HC3Nv7=1_219.5_&_HC3Nv6=v7_219.4']
                max_order = ['HC3N,v7=2','HC3N,v7=2', 'HC3N,v7=2', 'HC3N,v6=v7=1','HC3N,v6=v7=1','HC3N,v6=v7=1','HC3N,v6=v7=1']
                ax5.plot(spec['vel'], spec['int'], linewidth=0.8, drawstyle='steps-mid', color='k')
                ax5.plot(spec['vel'], spec['fit'], color='red', linewidth=0.8)
                if ii==(3-1):
                    local_max_x = [181.247, 228.133, 274.612, 487.308, 520.942, 599.394, 633.900]
                elif ii==(4-1) or ii==(5-1) or ii==(6-1) or ii==(7-1):
                    local_max_x = [181.247+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   228.133+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   274.612+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   487.308+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   520.942+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   599.394+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   633.900+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2]]
                elif ii==(1-1): 
                    local_max_x = [181.247+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   228.133+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   274.612+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   487.308+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   520.942+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   599.394+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   633.900+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2]]
                
                for h, vel in enumerate(local_max_x): 
                    ax5.axvline(local_max_x[h], color='k', linestyle='--', lw=0.5)
                    ax5.text(local_max_x[h], hc_max+hc_max*0.1, max_order[h], ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=9)
            
            elif j==1:
                lines = ['HC3Nv7=1_218.8_&_HC3Nv6=1_218.6']
                max_order = ['HC3N,v7=1','HC3N,v6=1','HC3N,v6=1']
                ax3.plot(spec['vel'], spec['int'], linewidth=0.8, drawstyle='steps-mid', color='k')
                ax3.plot(spec['vel'], spec['fit'], color='red', linewidth=0.9)
                if ii==(3-1):
                    local_max_x = [275.020, 283.767, 518.957]
                elif ii==(4-1) or ii==(5-1) or ii==(6-1) or ii==(7-1):
                    local_max_x = [275.020+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   283.767+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   518.957+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2]]
                    
                else:
                    local_max_x = [275.020+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   283.767+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   518.957+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2]]
                
                for h, vel in enumerate(local_max_x): 
                    ax3.axvline(local_max_x[h], color='k', linestyle='--', lw=0.5)
                    if h < 2:
                        ax3.text(local_max_x[h]-local_max_x[0]*0.12, hc_max+hc_max*0.1, max_order[0], ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=9)
                        ax3.text(local_max_x[h]+local_max_x[0]*0.1, hc_max+hc_max*0.1, max_order[1], ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=9)
                    else:
                        ax3.text(local_max_x[h], hc_max+hc_max*0.1, max_order[h], ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=9)
                        ax3.axvline(local_max_x[h], color='k', linestyle='--', lw=0.5)
            elif j==0:
                lines = ['HC3Nv6=v7_219.4_&_HC3Nv7=1_219.1']
                max_order = ['HC3N,v6=v7=1','HC3N,v6=v7=1','HC3N,v7=1']
                ax4.plot(spec['vel'], spec['int'], linewidth=0.8, drawstyle='steps-mid', color='k')
                ax4.plot(spec['vel'], spec['fit'], color='red', linewidth=0.8)
                hc3nv7_max.append(local_max_y[-1])
                if ii==(3-1):
                    local_max_x = [-84.86, -50.437, 276.305]
                    
                elif ii==(4-1) or ii==(5-1) or ii==(6-1) or ii==(7-1):
                    local_max_x = [-84.86+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   -50.437+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   276.305+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2]]
                elif ii==(1-1): 
                    local_max_xx = deepcopy(local_max_x)
                    local_max_x = [-84.86+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   -50.437+hc3n_prop['VLSR'][ii]-hc3n_prop['VLSR'][2],
                                   local_max_xx[-1]]
                    
                for h, vel in enumerate(local_max_x): 
                    ax4.axvline(local_max_x[h], color='k', linestyle='--', lw=0.5)
                    ax4.text(local_max_x[h], hc_max+hc_max*0.1, max_order[h], ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=9)
            if (ii==1-1):
                if j==3:
                    # H2CO
                    va_h2co = 5.246
                    vb_h2co = 352.338
                    va_cc3h2 = 28.725
                    vb_cc3h2 = 425.539
                    
                    ax2.axvline(va_h2co, color='k', linestyle='--', lw=0.5)
                    ax2.axvline(vb_h2co, color='k', linestyle='--', lw=0.5)
                    ax2.text(va_h2co-10, hc_max+hc_max*0.1, 'H2CO', ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=9)
                    ax2.text(vb_h2co, hc_max+hc_max*0.1, 'H2CO', ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=9)
                                
                    ax2.axvline(va_cc3h2, color='k', linestyle='--', lw=0.5)
                    ax2.text(va_cc3h2+10, hc_max+hc_max*0.1, 'C-C3H2', ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=9)
                    ax2.axvline(vb_cc3h2, color='k', linestyle='--', lw=0.5)
                    ax2.text(vb_cc3h2, hc_max+hc_max*0.1, 'C-C3H2', ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=9)
                
                if j==1:
                    vc_h2co = 350.035
                    ax3.axvline(vc_h2co, color='k', linestyle='--', lw=0.5)
                    ax3.text(vc_h2co, hc_max+hc_max*0.1, 'H2CO', ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=9)
                    v_hnco = 24.755
                    ax3.axvline(v_hnco, color='k', linestyle='--', lw=0.5)
                    ax3.text(v_hnco, hc_max+hc_max*0.1, 'HNCO', ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=9)
                    v_ocs = 141.064
                    ax3.axvline(v_ocs, color='k', linestyle='--', lw=0.5)
                    ax3.text(v_ocs, hc_max+hc_max*0.1, 'OCS', ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=9)
                    
                    
                if j==2:
                    v_co18 = 351.289
                    ax5.axvline(v_co18, color='k', linestyle='--', lw=0.5)
                    ax5.text(v_co18, hc_max+hc_max*0.1, 'CO-18', ha='center', va='center',rotation='vertical', backgroundcolor='white', fontsize=9)
                    
                    
        axis_spec = [ax2, ax3, ax4, ax5]
        
        for ax in axis_spec:
            ax.tick_params(direction='in')
            ax.xaxis.set_tick_params(top ='on')
            ax.yaxis.set_tick_params(right='on', labelright='off')
            ax.set_xlim([-250, 650])
            
            
            ax.set_ylim([hc_min, hc_max+hc_max*0.4])
            
            if ii==(1-1):
                ax.set_ylim([hc_min, 0.08])
            
            if ax == ax2 or ax == ax4:
                ax.set_ylabel(r'Jy/Beam', fontsize=14)
            if ax == ax4 or ax == ax5:
                ax.set_xlabel(r'v$_{rad}$ (km/s)', fontsize=14)
        
            ax.yaxis.set_label_coords(-0.12,0.5)
            
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.tick_params(axis='both', which='minor', labelsize=12)
            if ax == ax3 or ax == ax5:
                ax.set_yticklabels([])
            #if ax == ax2 or ax == ax3:
            #    ax.set_xticklabels([])
        ax3.set_xticklabels([])
        ax2.set_xticklabels([])
                
        fig.savefig(out_fig5_dir+'/vertical_ngc253_grid_HC'+str(i)+'_v'+version+'.png', bbox_inches='tight', transparent=True, dpi=400)
        plt.close()
            
    
#==============================================================================
# Luminosidades
#==============================================================================
    
# Simulated peaks HC3Nv7=1 219GHz
#hc3n_peak = []
#hc3nv6_peak = []
#for i in hot_cores:
#    hc3n_peak.append(local_max_y_list[i][0][-1])
#    hc3nv6_peak.append(local_max_y_list[i][1][-1])
        
if luminosidades == True:
        
    hc3n_prop['HC3Nv0_peak_JyBeam'] = hc3n_prop['hc3nv0_peak_mJy/beam']/1000.
    hc3n_prop['HC3Nv7_peak_JyBeam'] = hc3n_prop['hc3nv7_peak_mJy/beam']/1000. 
    hc3n_prop['HC3Nv6_peak_JyBeam'] = hc3n_prop['hc3nv6_peak_mJy/beam']/1000. 
       
    freq_v0 =   218.324723 # GHz
    freq_v7 =   219.173757 # GHz
    freq_v6 =   218.854392 # GHz
    
    # Brightness temperature
    hc3n_prop['T_B_v0'] = u_conversion.Jybeam_to_T(hc3n_prop['HC3Nv0_peak_JyBeam'], freq_v0, bmin.value, bmaj.value)
    hc3n_prop['T_B_v7'] = u_conversion.Jybeam_to_T(hc3n_prop['HC3Nv7_peak_JyBeam'], freq_v7, bmin.value, bmaj.value)
    hc3n_prop['T_B_v6'] = u_conversion.Jybeam_to_T(hc3n_prop['HC3Nv6_peak_JyBeam'], freq_v6, bmin.value, bmaj.value)
    
    # Source size
    hc3n_prop['Source_size_v0'] = utiles.sourcesize_from_fit(hc3n_prop['T_B_v0'], hc3n_prop['Tvib'], bmin.value*bmaj.value)
    hc3n_prop['Source_size_m_v0'] = u_conversion.lin_size(D,hc3n_prop['Source_size_v0'])
    hc3n_prop['Source_size_v7'] = utiles.sourcesize_from_fit(hc3n_prop['T_B_v7'], hc3n_prop['Tvib'], bmin.value*bmaj.value)
    hc3n_prop['Source_size_m_v7'] = u_conversion.lin_size(D,hc3n_prop['Source_size_v7'])
    hc3n_prop['Source_size_v6'] = utiles.sourcesize_from_fit(hc3n_prop['T_B_v6'], hc3n_prop['Tvib'], bmin.value*bmaj.value)
    hc3n_prop['Source_size_m_v6'] = u_conversion.lin_size(D,hc3n_prop['Source_size_v6'])
    
    # Line Luminosities
    hc3n_prop['L_Watt_v0'] = u_conversion.stef_boltz(hc3n_prop['Source_size_m_v0']/2., hc3n_prop['Tvib'])
    hc3n_prop['L_Watt_v0_err'] = u_conversion.stef_boltz_error(hc3n_prop['Source_size_m_v0']/2., 0, hc3n_prop['Tvib'], hc3n_prop['Tvib_err'])
    hc3n_prop['L_Lsun_v0'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v0'])
    hc3n_prop['L_Lsun_v0_err'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v0_err'])
    
    
    hc3n_prop['L_Watt_v7'] = u_conversion.stef_boltz(hc3n_prop['Source_size_m_v7']/2., hc3n_prop['Tvib'])
    hc3n_prop['L_Watt_v7_err'] = u_conversion.stef_boltz_error(hc3n_prop['Source_size_m_v7']/2., 0, hc3n_prop['Tvib'], hc3n_prop['Tvib_err'])
    hc3n_prop['L_Lsun_v7'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v7'])
    hc3n_prop['L_Lsun_v7_err'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v7_err'])
    
    
    hc3n_prop['L_Watt_v6'] = u_conversion.stef_boltz(hc3n_prop['Source_size_m_v6']/2., hc3n_prop['Tvib'])
    hc3n_prop['L_Lsun_v6'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v6'])
    
    # Size in parsecs
    hc3n_prop['Source_size_pc_v7'] = hc3n_prop['Source_size_m_v7'] * 3.2407793E-17
    
    # Virial Mass
    hc3n_prop['M_vir'] = utiles.virial_mass(hc3n_prop['Source_size_pc_v7'], hc3n_prop['v'])
    # Rotational temperatures HC3N_v7=1
    
    # HC3N_v7=1 upper
    #Eu = 460.252 # Elo cm-1
    #Ju = 40.
    #Nu = 
    #freq_u = 364.676275 # GHz
    
    # HC3N_v7=1 lower
    #El = 307.081 # Elo cm-1
    #Jl = 24. 
    #Nl = 
    #freq_l = 219.173757 # GHz
    
    
    
    #Trotational(16.3, 16.6, 460.252, 307.081, 40, 24)
    
    
    
    # Optical depths
    hc3n_prop['tau_v7'] = utiles.tau_from_T(hc3n_prop['T_B_v7'], hc3n_prop['Tvib'])
    hc3n_prop['tau_v0'] = utiles.tau_from_T(hc3n_prop['T_B_v0'], hc3n_prop['Tvib'])
    
    table_out = ['Ind_2', 'RA', 'Dec', 'VLSR', 'VLSR_err', 'v', 'v_err',
                 'N(HC3N)', 'N(HC3N)_err', 'Tvib', 'Tvib_err', 'Trot', 'Trot_err',
                 'Source_size_v7', 'hc3nv7_peak_mJy/beam', 'hc3nv6_peak_mJy/beam',
                 'L_Lsun_v7', 'L_Lsun_v7_err']
    
    hc3n_prop.to_latex(buf=None, columns=table_out)

#hc3n_prop['tau_MC'] = [0.182, 0.104, 0.019, 0.019, 0.011, 0.057, 0.015, 0.061, 0.040, 0.045, 0.046, 0.026] 


## Mass from column density
#distance = 3.5 * u.Mpc
#d_cm = distance.to(u.cm)
#muh = 2.8
#mh = 1.6733E-24

#M_t = hc3n_prop['Source_size'] * (d_cm**2.) * hc3n_prop['NH2'] * muh * _si.u #muh*mh
#msun_kg = (1 * u.M_sun).to(u.kg)
#hc3n_prop['M_tot'] = M_t * u.M_sun / msun_kg

#hc3n_prop.to_csv(out_dir+'/hot_cores_v'+version+'.txt', sep='\t', index=False)

#import sys
#sys.exit("Error message")



    








