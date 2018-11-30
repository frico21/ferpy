#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  12 13:27:06 2018

@author: frico
"""
version = '2'

"""
He comentado las lineas 79-82 de colis/hc3n_h2_faure16.f
"""
#%matplotlib inline
#%matplotlib qt

from ferpy import radtrasf_hc3n_runner_v6 as radtransf_hc3n
from ferpy import radtransf_modelpars_hc3n_v6 as radtransf_modelpars
from ferpy import radtransf_HC3N_plot_and_luminosity_v6 as hc3n_plotter
from ferpy import u_conversion
from ferpy import utiles
from ferpy import radtransf_model_calculations as m_calc

from copy import deepcopy
import scipy
import numpy as np
import os
import pandas as pd
import astropy.units as u
import astropy.constants.si as _si
from astropy.io import fits
from astropy.wcs import WCS
from astropy import wcs
#from radio_beam import Beam
import resource # to control memory usage
import gc

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits import mplot3d
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.colors



from glob import glob





#import plotly
#plotly.tools.set_credentials_file(username='igreen21', api_key='cu4ClmelAvMBdD1hhAwl')
import plotly.plotly as py
import plotly.offline as offline
#import plotly.graph_objs as go
import colorlover as cl
import scipy.stats as st


#import plotly.offline as py
from plotly.offline import download_plotlyjs, init_notebook_mode, iplot

import plotly.graph_objs as go
init_notebook_mode()

#import seaborn as sns

class data:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        

plot_nh_ch = False
plot_td_ch = True   



# Observed data workdir
dworkdir_cubes = '/Users/frico/Documents/data/NGC253_H3O+'
dworkdir_spec = dworkdir_cubes+'/Hotcores_v4_all'

# Out dir
out_dir = dworkdir_spec+'/Results_v'+version+'/'
if not os.path.exists(out_dir):
        os.makedirs(out_dir)

# NGC253 properties
D = 3.5 # Mpc to NGC253
vsis = 258 #km/s


def read_hc3n_obs(hc3n_path):
    # Observation data
    hc3n_prop = pd.read_csv('/Users/frico/Documents/data/NGC253_H3O+/'+hc3n_path, delim_whitespace= True, header=0, comment='#')
    
    rms= 0.001383090342853783*1000 # mJy/beam
    rms_219 = 1.567
    rms_350 = 1.1256
    rms_307 = 0.956
    hc3n_prop['ratio_24_23'] =hc3n_prop['HC3Nv0_218_mJyBeam']/hc3n_prop['HC3Nv7_219_mJyBeam']
    hc3n_prop['ratio_24_23_err'] = np.sqrt((rms_219/hc3n_prop['HC3Nv7_219_mJyBeam'])**2+ (hc3n_prop['HC3Nv0_218_mJyBeam']*rms_219/hc3n_prop['HC3Nv7_219_mJyBeam']**2)**2)
    
    hc3n_prop['ratio_39_38'] = hc3n_prop['HC3Nv0_354_mJyBeam']/hc3n_prop['HC3Nv7_355_mJyBeam']
    hc3n_prop['ratio_39_38_err'] = np.sqrt((rms_350/hc3n_prop['HC3Nv7_355_mJyBeam'])**2+ (hc3n_prop['HC3Nv0_354_mJyBeam']*rms_350/hc3n_prop['HC3Nv7_355_mJyBeam']**2)**2)
    
    hc3n_prop['ratio_v7'] = hc3n_prop['HC3Nv7_355_mJyBeam']/hc3n_prop['HC3Nv7_219_mJyBeam']
    hc3n_prop['ratio_v7_err'] = np.sqrt((rms_350/hc3n_prop['HC3Nv7_219_mJyBeam'])**2+ (hc3n_prop['HC3Nv7_355_mJyBeam']*rms_219/hc3n_prop['HC3Nv7_219_mJyBeam']**2)**2)
    
    hc3n_prop['ratio_v0'] = hc3n_prop['HC3Nv0_354_mJyBeam']/hc3n_prop['HC3Nv0_218_mJyBeam']
    hc3n_prop['ratio_v0_err'] = np.sqrt((rms_219/hc3n_prop['HC3Nv0_218_mJyBeam'])**2+ (hc3n_prop['HC3Nv0_354_mJyBeam']*rms_350/hc3n_prop['HC3Nv0_218_mJyBeam']**2)**2)
    # To have same name as modelled data 
    hc3n_prop[['peak_v0_24_23', 'peak_v7_24_23', 'peak_v0_39_38', 'peak_v7_39_38']] = hc3n_prop[['HC3Nv0_218_mJyBeam', 'HC3Nv7_219_mJyBeam', 'HC3Nv0_354_mJyBeam', 'HC3Nv7_355_mJyBeam']]
    
    
   
    hc3n_prop.to_csv('/Users/frico/Documents/data/NGC253_H3O+/'+'Hotcores_v4_all/hc3n_obs_results_v0_excel.txt', header=True, sep='\t', index=False, na_rep=-1)

    return hc3n_prop

hc3n_path ='Hotcores_v4_all'+'/hc3n_obs_results_v0.txt'
hc3n_prop = read_hc3n_obs(hc3n_path)

# Plotting General Pars
fig_format = 'png'
factor_names = ['07', '035']
linestyles_list = ['-', '--','--', '-.','-.', ':', ':']

 

def hc3n_plot_starter(m, n, hc3n_prop, y_axis_obs, hc_pos, yaxis, xaxis):
    #m = 2
    #n = 2
    # Figure Size
    size_rat = float(n)/float(m)
    size_x = 10.*size_rat
    size_y = 10.
    fig = plt.figure(figsize=(size_x, size_y))
    gs1 = gridspec.GridSpec(m, n)    
    gs1.update(wspace = 0.0, hspace=0.0, top=0.95, bottom = 0.05)
    # Generating specified number of axis
    axis = []
    axis_ind = []
    for i in range(m*n):
        axis.append(fig.add_subplot(gs1[i]))
        if xaxis == 'central_nh':
            axis[i].set_xscale('log')
            
        if yaxis != 'taus':
            # Plotting obs ratios
            for l, line in hc3n_prop.iterrows():
                if line['Ind_ok'] not in [6, 7]:#[9, 10, 11, 12]:
                    if xaxis == 'central_nh':
                        xpos = hc_pos
                    elif xaxis == 'central_T':
                        xpos = line['Tvib']
                    if not pd.isnull(line[y_axis_obs[i]]):
                        axis[i].text(xpos, line[y_axis_obs[i]], str(l+1), ha='center', va='center', color='k', fontsize=8)

    return fig, axis
            

def hc3n_plot_starter_cont(m, n, hc3n_prop, y_axis_obs, z_axis_obs, hc_pos, yaxis, xaxis):
    #m = 2
    #n = 2
    # Figure Size
    size_rat = float(n)/float(m)
    size_x = 10.*size_rat
    size_y = 10.
    fig = plt.figure(figsize=(size_x, size_y))
    gs1 = gridspec.GridSpec(m, n)    
    gs1.update(wspace = 0.0, hspace=0.0, top=0.95, bottom = 0.05)
    # Generating specified number of axis
    axis = []
    axis_ind = []
    for i in range(m*n):
        axis.append(fig.add_subplot(gs1[i]))
        #if xaxis == 'central_nh':
        #    axis[i].set_xscale('log')
            
        #if yaxis != 'taus':
        # Plotting obs ratios
    for l, line in hc3n_prop.iterrows():
        if line['Ind_ok'] not in [9, 10, 11, 12]:
            if not pd.isnull(line[y_axis_obs[3]]):
                axis[2].text(line[z_axis_obs[0]], line[y_axis_obs[3]], str(l+1), ha='center', va='center', color='k', fontsize=8)

    return fig, axis

def hc3n_model_ratio_cal(results_path, factor):
    summary_path = '/models_summary_'+factor+'.txt'
    hc3n_model = pd.read_csv(results_path+summary_path, delim_whitespace= True, header=0, comment='#')
    hc3n_model.drop_duplicates(inplace=True)

    h = 0
    # Calculating ratios from models
    for hotcore in hc3n_prop.itertuples():
        h +=1
        ratios_obs = [hotcore.ratio_24_23, hotcore.ratio_39_38, hotcore.ratio_v7, hotcore.ratio_v0]  
        ratios_obs_err = [hotcore.ratio_24_23_err, hotcore.ratio_39_38_err, hotcore.ratio_v7_err, hotcore.ratio_v0_err]
        # Obtaining chisq
        chisq_24_23 =(ratios_obs[0] - hc3n_model['ratio_24_23'])**2/(ratios_obs_err[0]**2)
        chisq_39_38 =(ratios_obs[1] - hc3n_model['ratio_39_38'])**2/(ratios_obs_err[1]**2)
        chisq_v7 =(ratios_obs[2] - hc3n_model['ratio_v0'])**2/(ratios_obs_err[2]**2)
        chisq_v0 =(ratios_obs[3] - hc3n_model['ratio_v7'])**2/(ratios_obs_err[3]**2)
        chisq_all = (chisq_24_23 + chisq_39_38 + chisq_v7 + chisq_v0)
        hc3n_model['chisq_'+str(h)] = chisq_all
    return hc3n_model
            

def hc3n_legend_ymax_ymin(axis, pane_names, set_ylim, min_list, max_list, custom_lines, custom_lines_txt, custom_points, custom_points_txt, cont_plot):
    axis[0].legend(custom_lines, custom_lines_txt)
    axis[1].legend(custom_points, custom_points_txt)
    
    for i, axx in enumerate(axis):
        axis[i].text(0.1, 0.9, pane_names[i],
                            horizontalalignment='left',
                            verticalalignment='top',
                            transform=axis[i].transAxes)
    if cont_plot==False:
        if set_ylim != 0:
            # Setting ylims according to max values
            for mmm, maxi in enumerate(max_list):
                if maxi > set_ylim:
                    max_list[mmm] = set_ylim # Max y for line ratios
            # Setting ylims 
            for i, aax in enumerate(axis):
                if min_list[i] != np.NaN and min_list[i] != np.Inf and min_list[i] != -np.Inf:
                    minn = min_list[i]
                else:
                    minn = 0
                if max_list[i] != np.NaN and max_list[i] != np.Inf and max_list[i] != -np.Inf:
                    maxx = max_list[i]
                else:
                    maxx = 60
                axis[i].set_ylim([minn-0.5, maxx+0.5])
                    #axis[2].set_ylim([min_list[2]-0.5, max_list[2].max()+0.5])
    return max_list

def hc3n_mod_plot(n, axis, xaxis, hc3n_model, dens_colors, marker_list, linestyles_list, y_axis_obs, min_list, max_list, t_grad):
    custom_lines = []       # For custom legend
    custom_lines_txt = []   # For custom legend text
    if xaxis == 'central_nh':
        list_start = T_list
        lp_lab = 'T'
        formatter = '%1.f'
        slicer = 'central_T'
    elif xaxis == 'central_T':
        list_start = start_dens_list
        lp_lab = 'nH'
        formatter = '%1.2E'
        slicer = 'central_nh'
    
    for s, tkin in enumerate(list_start):
        print '\t\t'+lp_lab+'=\t'+formatter % tkin
        nline_color = dens_colors[s]
        # Subsetting dataframe to corresponging Dens or Temp
        ncol = hc3n_model[(hc3n_model[slicer]/tkin > 0.95) & (hc3n_model[slicer]/tkin < 1.05)]
        ncol = ncol.sort_values(by=[xaxis], ascending = True)
        # plotting
        axis[0].plot(ncol[xaxis], ncol[y_axis_obs[0]],
                color=dens_colors[s], linestyle=linestyles_list[n], lw=0.5,
                marker=marker_list[0], mfc=nline_color, mec=nline_color, markersize=4)
        axis[1].plot(ncol[xaxis], ncol[y_axis_obs[1]],
                color=dens_colors[s], linestyle=linestyles_list[n], lw=0.5,
                marker=marker_list[0], mfc=nline_color, mec=nline_color, markersize=4)
        axis[2].plot(ncol[xaxis], ncol[y_axis_obs[2]],
                color=dens_colors[s], linestyle=linestyles_list[n], lw=0.5,
                marker=marker_list[0], mfc=nline_color, mec=nline_color, markersize=4)
        axis[3].plot(ncol[xaxis], ncol[y_axis_obs[3]],
                color=dens_colors[s], linestyle=linestyles_list[n], lw=0.5,
                marker=marker_list[0], mfc=nline_color, mec=nline_color, markersize=4)
        # Comparing with obs for ylims
        for r, ratio in enumerate(y_axis_obs):
            max_now = ncol[ratio].max()
            min_now = ncol[ratio].min()
            if max_list[r] < max_now:
                max_list[r] = max_now
            if min_list[r] > min_now:
                min_list[r]= min_now
        # Custom Legend for dens or T
        custom_lines.append(Line2D([0], [0], color=dens_colors[s], lw=0.5), )
        custom_lines_txt.append(lp_lab+'='+formatter % tkin)
    return custom_lines, custom_lines_txt, min_list, max_list

def hc3n_mod_plot_contour(n, axis, xaxis, hc3n_model, dens_colors, marker_list, linestyles_list, y_axis_obs, z_axis_obs, min_list, max_list, t_grad, T_list):
    custom_lines = []       # For custom legend
    custom_lines_txt = []   # For custom legend text
    if xaxis == 'central_nh':
        list_start = T_list
        lp_lab = 'T'
        formatter = '%1.f'
        slicer = 'central_T'
    elif xaxis == 'central_T':
        list_start = start_dens_list
        lp_lab = 'nH'
        formatter = '%1.2E'
        slicer = 'central_nh'
    
    #list_start = [400]
    for s, tkin in enumerate(list_start):
        print '\t\t'+lp_lab+'=\t'+formatter % tkin
        nline_color = dens_colors[s]
        # Subsetting dataframe to corresponging Dens or Temp
        ncol = hc3n_model[(hc3n_model[slicer]/tkin > 0.95) & (hc3n_model[slicer]/tkin < 1.05)]
        ncol = ncol.sort_values(by=[xaxis], ascending = True)
        
        x = ncol[z_axis_obs[0]]
        y = ncol[y_axis_obs[3]]
        z = ncol[xaxis]
        
        # Set up a regular grid of interpolation points
        xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
        xi, yi = np.meshgrid(xi, yi)
        # Interpolate
        rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
        zi = rbf(xi, yi)
        
        
        axis[2].contourf(xi,yi,zi)
        #im = axis[0].imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
        #   extent=[x.min(), x.max(), y.min(), y.max()])
        axis[2].scatter(x, y, c='k', s=4)
        #fig.colorbar(im, ax=axis[0])
        #plt.colorbar()
       
        
        
        # plotting
#        axis[0].plot(ncol[z_axis_obs[0]], ncol[y_axis_obs[0]],
#                color=dens_colors[s], linestyle=None,
#                marker=marker_list[0], mfc=nline_color, mec=nline_color, markersize=4)
#        axis[1].plot(ncol[z_axis_obs[0]], ncol[y_axis_obs[1]],
#                color=dens_colors[s], linestyle=None,
#                marker=marker_list[0], mfc=nline_color, mec=nline_color, markersize=4)
#        axis[2].plot(ncol[z_axis_obs[0]], ncol[y_axis_obs[2]],
#                color=dens_colors[s], linestyle=None,
#                marker=marker_list[0], mfc=nline_color, mec=nline_color, markersize=4)
#        axis[3].plot(ncol[z_axis_obs[0]],  ncol[y_axis_obs[3]],
#                color=dens_colors[s], linestyle=None,
#                marker=marker_list[0], mfc=nline_color, mec=nline_color, markersize=4)
        # Comparing with obs for ylims
        for r, ratio in enumerate(y_axis_obs):
            max_now = ncol[ratio].max()
            min_now = ncol[ratio].min()
            if max_list[r] < max_now:
                max_list[r] = max_now
            if min_list[r] > min_now:
                min_list[r]= min_now
        # Custom Legend for dens or T
        custom_lines.append(Line2D([0], [0], color=dens_colors[s], lw=0.5), )
        custom_lines_txt.append(lp_lab+'='+formatter % tkin)
    return custom_lines, custom_lines_txt, min_list, max_list



def luminosity_calculator(model_path, model_name, factor, plot_profiles, plot_SED, figoutdir, fig_format, obs_df, Temp):
    # fig_format = 'png'
    linewidth  = 0.7
    intensity_limits=['NA', 'NA']
    velocity_limits=[-100, 100]
    label_y =  r"Flux density (mJy)" 
    label_x = r"V (km/s)"
    labelsize = 12
    filas = 2
    columnas = 1
    first_last_xlabel = False
    v0_color = 'red'
    v7_color = 'k'
    
    
    # =============================================================================
    # # Line profiles
    # =============================================================================
    dist_pc = 3.5e6 # pc
    luminosity_obs = 2.928e9 # lsun
    model_size = 1. * factor # pc 
    
    if plot_profiles == True or plot_SED==True:
        if not os.path.exists(figoutdir):
            os.makedirs(figoutdir)
            
    line_names = [r'HC$_3$N 24(+1)-23(-1)', r'HC$_3$N 39(+1)-38(-1)']
    
    # Loading model line profile
    model_prof = pd.read_csv(model_path+'/'+model_name+'_.prof', delim_whitespace= True, skiprows=[0])
    model_prof.columns = ['lambda_um', 'I_ergscm2um', 'I_Ic_ergscm2um', 'V_kms']
    
    # Conversion erg/s/cm^2/um to mJy
    model_prof['I_Ic_mJy'] = u_conversion.ergscmum_to_mjy(model_prof['I_Ic_ergscm2um'], model_prof['lambda_um'])
    # Multiplying intensities by factor
    model_prof['I_Ic_mJy'] = model_prof['I_Ic_mJy']*(factor**2)
    
    
    #### J=24-23
    
    ### v=0
    
    alam_v0_24_23 = 1373.1495 # v0 24(+1)-23(-1)
    prof_v0_24_23_lim = []
    for i, row in model_prof.iterrows():
        if row['V_kms'] == 0 and abs(row['lambda_um']-alam_v0_24_23)<0.1:
            central_v0_24_23_index = i
            for j in range(i-1, 1, -1):
                if model_prof['V_kms'][j]>0:
                    prof_v0_24_23_lim.append(j+1)
                    break
            for k in range(i+1, model_prof.shape[0]):
                if model_prof['V_kms'][k] < 0:
                    prof_v0_24_23_lim.append(k-1)
                    break
                
    prof_v0_24_23 = model_prof.loc[prof_v0_24_23_lim[0]:prof_v0_24_23_lim[1]]
    peak_v0_24_23 = prof_v0_24_23['I_Ic_mJy'][prof_v0_24_23.V_kms == 0].max()
                
    ### v7=1
    
    alam_v7_24_23 = 1367.8303 # v7 24(+1)-23(-1)
    prof_v7_24_23_lim = []
    for i, row in model_prof.iterrows():
        if row['V_kms'] == 0 and abs(row['lambda_um']-alam_v7_24_23)<0.1:
            central_v7_24_23_index = i
            for j in range(i-0, 1, -1):
                if model_prof['V_kms'][j]>0:
                    prof_v7_24_23_lim.append(j+1)
                    break
            for k in range(i+1, model_prof.shape[0]):
                if model_prof['V_kms'][k] < 0:
                    prof_v7_24_23_lim.append(k-1)
                    break
    
    prof_v7_24_23 = model_prof.loc[prof_v7_24_23_lim[0]:prof_v7_24_23_lim[1]]
    peak_v7_24_23 = prof_v7_24_23['I_Ic_mJy'][prof_v7_24_23.V_kms == 0].max()
    
    #### J=39-38
    
    ### v0
    alam_v0_39_38 = 845.2062 # v0 39(+1)-38(-1)
    prof_v0_39_38_lim = []
    for i, row in model_prof.iterrows():
        if row['V_kms'] == 0 and abs(row['lambda_um']-alam_v0_39_38)<0.1:
            central_v0_39_38_index = i
            for j in range(i-0, 1, -1):
                if model_prof['V_kms'][j]>0:
                    prof_v0_39_38_lim.append(j+1)
                    break
            for k in range(i+1, model_prof.shape[0]):
                if model_prof['V_kms'][k] < 0:
                    prof_v0_39_38_lim.append(k-1)
                    break
    prof_v0_39_38 = model_prof.loc[prof_v0_39_38_lim[0]:prof_v0_39_38_lim[1]]
    peak_v0_39_38 = prof_v0_39_38['I_Ic_mJy'][prof_v0_39_38.V_kms == 0].max()
    
    ### v7=1
    alam_v7_39_38 = 843.1411 # v0 37(+1)-38(-1)
    prof_v7_39_38_lim = []
    for i, row in model_prof.iterrows():
        if row['V_kms'] == 0 and abs(row['lambda_um']-alam_v7_39_38)<0.1:
            central_v7_39_38_index = i
            for j in range(i-0, 1, -1):
                if model_prof['V_kms'][j]>0:
                    prof_v7_39_38_lim.append(j+1)
                    break
            for k in range(i+1, model_prof.shape[0]):
                if model_prof['V_kms'][k] < 0:
                    prof_v7_39_38_lim.append(k-1)
                    break
    prof_v7_39_38 = model_prof.loc[prof_v7_39_38_lim[0]:prof_v7_39_38_lim[1]]
    # peak_v7_39_38 = prof_v7_39_38['I_Ic_mJy'].max() old
    peak_v7_39_38 = prof_v7_39_38['I_Ic_mJy'][prof_v7_39_38.V_kms == 0].max()
    
    # Modeled values 
    hc3n_v0_ratio = peak_v0_39_38/peak_v0_24_23
    hc3n_v7_ratio = peak_v7_39_38/peak_v7_24_23
    hc3n_v0v7_3938_ratio = peak_v0_39_38/peak_v7_39_38
    hc3n_v0v7_2423_ratio = peak_v0_24_23/peak_v7_24_23
    
    out_dir_prof = figoutdir+'/'+'Profiles'
    if not os.path.exists(out_dir_prof):
        os.makedirs(out_dir_prof)
        
    out_dir_sed = figoutdir+'/'+'SED'
    if not os.path.exists(out_dir_sed):
        os.makedirs(out_dir_sed)
    
    # Plotting line profiles
    if plot_profiles == True:
        m = filas
        n = columnas
        dx, dy = 1.2, 1
        figsize = plt.figaspect(float(dy * m) / float(dx * n))
        fig_PROF = plt.figure(figsize=figsize)
        gs1 = gridspec.GridSpec(m, n)    
        gs1.update(wspace = 0.0, hspace=0.0, top=0.95, bottom = 0.05)
        # Generating specified number of axis
        axis = []
        axis_ind = []
        for i in range(m*n):
            axis.append(fig_PROF.add_subplot(gs1[i]))
        # Generating axis index
        ind = 0
        axis_ind = []
        for i  in range(m):
            axis_ind.append([])
            for j in range(n):
                axis_ind[i].append(ind)
                ind += 1
        for j in range(m*n):
            # Plotting J=24(+1)-23(-1)
            if j == 0:
                # Plotting spectrum
                axis[j].plot(prof_v0_24_23['V_kms'], prof_v0_24_23['I_Ic_mJy'], linewidth=linewidth, color=v0_color)
                axis[j].plot(prof_v7_24_23['V_kms'], prof_v7_24_23['I_Ic_mJy'], linewidth=linewidth, color=v7_color)
                # Line name
                axis[j].text(0.1, 0.9,
                            line_names[j], ha='left', va='center',
                            transform=axis[j].transAxes,
                            rotation='horizontal', fontsize=8)
                axis[j].text(0.1, 0.85,
                            r'HC$_3$N,v7=1 ratio='+'%1.2f' % hc3n_v7_ratio, ha='left', va='center',
                            transform=axis[j].transAxes,
                            rotation='horizontal', fontsize=8)
                
                axis[j].text(0.1, 0.80,
                            r'HC$_3$N,v=0 ratio='+'%1.2f' % hc3n_v0_ratio, ha='left', va='center',
                            transform=axis[j].transAxes,
                            rotation='horizontal', fontsize=8)
                axis[j].text(0.1, 0.75,
                            r'HC$_3$N_24-23 ratio='+'%1.2f' % hc3n_v0v7_2423_ratio, ha='left', va='center',
                            transform=axis[j].transAxes,
                            rotation='horizontal', fontsize=8)
                
                # Obs
                obs_vs = np.linspace(-50, 50, 14)
                for o, orow in obs_df.iterrows():
                    #if orow['peak_v0_39_38'] != np.nan:
                    if not pd.isnull(orow['peak_v0_24_23']):
                        axis[j].text(obs_vs[o], orow['peak_v0_24_23'],
                                str(orow['Ind_ok']), ha='left', va='center', color=v0_color,
                                #transform=axis[j].transAxes,
                                rotation='horizontal', fontsize=8)
                        
                    if not pd.isnull(orow['peak_v7_24_23']):
                        axis[j].text(obs_vs[o], orow['peak_v7_24_23'],
                                str(orow['Ind_ok']), ha='left', va='center', color=v7_color,
                                #transform=axis[j].transAxes,
                                rotation='horizontal', fontsize=8)
    
            
    
            elif j == 1:
                # Plotting spectrum
                axis[j].plot(prof_v0_39_38['V_kms'], prof_v0_39_38['I_Ic_mJy'], linewidth=linewidth, color=v0_color)
                axis[j].plot(prof_v7_39_38['V_kms'], prof_v7_39_38['I_Ic_mJy'], linewidth=linewidth, color=v7_color)
                # Line name
                axis[j].text(0.1, 0.9,
                            line_names[j], ha='left', va='center',
                            transform=axis[j].transAxes,
                            rotation='horizontal', fontsize=8)
                axis[j].text(0.1, 0.75,
                            r'HC$_3$N_39-38 ratio='+'%1.2f' % hc3n_v0v7_3938_ratio, ha='left', va='center',
                            transform=axis[j].transAxes,
                            rotation='horizontal', fontsize=8)
                
                # Obs
                obs_vs = np.linspace(-50, 50, 14)
                for o, orow in obs_df.iterrows():
                    #if orow['peak_v0_39_38'] != np.nan:
                    if not pd.isnull(orow['peak_v0_39_38']):
                        axis[j].text(obs_vs[o], orow['peak_v0_39_38'],
                                str(orow['Ind_ok']), ha='left', va='center', color=v0_color,
                                #transform=axis[j].transAxes,
                                rotation='horizontal', fontsize=8)
                        
                    if not pd.isnull(orow['peak_v7_39_38']):
                        axis[j].text(obs_vs[o], orow['peak_v7_39_38'],
                                str(orow['Ind_ok']), ha='left', va='center', color=v7_color,
                                #transform=axis[j].transAxes,
                                rotation='horizontal', fontsize=8)
                
        # Left Figures
        left_ind = []
        for i in range(m):
            left_ind.append(axis_ind[i][0])
        # Bottom Figures
        bottom_figures = axis_ind[-1]
        velocity_limits_total = [np.nanmin([prof_v0_24_23_lim[0], prof_v7_24_23_lim[0], prof_v0_39_38_lim[0], prof_v7_39_38_lim[0]]),
                                 np.nanmax([prof_v0_24_23_lim[1], prof_v7_24_23_lim[1], prof_v0_39_38_lim[1], prof_v7_39_38_lim[1]])]
        
        intensity_limits_total = [-0.05,
                                 np.nanmax([prof_v0_24_23['I_Ic_mJy'].max(),
                                            prof_v0_39_38['I_Ic_mJy'].max(),
                                            prof_v7_24_23['I_Ic_mJy'].max(),
                                            prof_v7_39_38['I_Ic_mJy'].max(),
                                            obs_df['peak_v0_39_38'].max(),
                                            obs_df['peak_v7_39_38'].max(),
                                            obs_df['peak_v0_24_23'].max(),
                                            obs_df['peak_v7_24_23'].max()])]
        
        for i, ax in enumerate(axis):
            #ax.tick_params(direction='in')
            
            majorLocator = MultipleLocator(20)
            ax.yaxis.set_major_locator(majorLocator)
            # Solo ponemos titulo a los ejes y si estan a la izquierda del todo
            ax.minorticks_on()
            #ax.tick_params(axis='both', which='minor', direction='in', top=True, right=True)
            ax.tick_params(axis='both', which='major', labelsize=labelsize, length=5, direction='in')
            ax.tick_params(axis='both', which='minor', labelsize=labelsize, direction='in')
            ax.xaxis.set_tick_params(top ='on', which='both')
            ax.yaxis.set_tick_params(right='on', which='both', labelright='off')
            
            # Limite en el eje x
            if velocity_limits != ['NA', 'NA']: # Si no se especifican limites python coge el minimo y maximo de cada expectro
                ax.set_xlim(velocity_limits)    # Limites introducidos
            else:
                ax.set_xlim(velocity_limits_total) # Todos con los mismos limites
            # Limite en el eje y
            if intensity_limits != ['NA', 'NA']: # Si no especificamos limites ejey coge el maximo y minimo disponible del espectro
                ax.set_ylim(intensity_limits)
            else:
                ax.set_ylim([intensity_limits_total[0], intensity_limits_total[1]+intensity_limits_total[1]*0.1]) # Para que el maximo en el eje y sea el mismo en todos los espectros
            ax.set_ylim([intensity_limits_total[0], intensity_limits_total[1]+intensity_limits_total[1]*0.4])
            if i in left_ind:
                ax.set_ylabel(label_y, fontsize=14)
            else:
                ax.set_yticklabels([])
            # Solo ponemos titulo a los ejes x si estan abajo del todo
            if i in bottom_figures:
                ax.set_xlabel(label_x, fontsize=14)
                if first_last_xlabel == False:
                    plt.setp(ax.get_xticklabels()[0], visible=False)    
                    plt.setp(ax.get_xticklabels()[-1], visible=False)
            else:
                ax.set_xticklabels([])
                
            #majorLocator = MultipleLocator(5)
            #minorLocator = MultipleLocator(5)
            # for the minor ticks, use no labels; default NullFormatter
            #ax.xaxis.set_minor_locator(minorLocator)
    
        #fig_PROF.savefig(figoutdir+'/'+model_name+'_f'+'%1.1f' % factor +'_profile.'+fig_format, bbox_inches='tight', transparent=True, dpi=400, format=fig_format)
        #plt.close()
        fig_PROF.savefig(out_dir_prof+'/'+model_name+'_f'+'%1.1f' % factor +'_T'+'%1.f' % Temp+'_profile.'+fig_format, bbox_inches='tight', transparent=True, dpi=400, format=fig_format)
        plt.close(fig_PROF)
            
    
    
    
    # =============================================================================
    # ## SED and tau_continuum
    # =============================================================================
        
    
    ## ISO
    model_iso = pd.read_csv(model_path+'/'+model_name+'_.iso', delim_whitespace= True)
    model_iso.columns = ['lambda_um', 'I_Ic_ergscm2um', 'Ic_ergscm2um', 'Inorm_ergscm2um', 'I_ergscm2um']
    # Conversion erg/s/cm^2/um to mJy
    model_iso['Ic_Jy'] = u_conversion.ergscmum_to_mjy(model_iso['Ic_ergscm2um'], model_iso['lambda_um'])/1000
    # Multiplying intensities by factor
    model_iso['Ic_Jy'] = model_iso['Ic_Jy']*(factor**2)
    model_iso['Ic_ergscm2um'] = model_iso['Ic_ergscm2um']*(factor**2)
    
    ## Spire
    model_spire = pd.read_csv(model_path+'/'+model_name+'_.spire', delim_whitespace= True)
    model_spire.columns = ['lambda_um', 'I_Ic_ergscm2um', 'Ic_ergscm2um', 'Inorm_ergscm2um', 'I_ergscm2um']
    # Conversion erg/s/cm^2/um to mJy
    model_spire['Ic_Jy'] = u_conversion.ergscmum_to_mjy(model_spire['Ic_ergscm2um'], model_spire['lambda_um'])/1000
    # Multiplying intensities by factor
    model_spire['Ic_Jy'] = model_spire['Ic_Jy']*(factor**2)
    model_spire['Ic_ergscm2um'] = model_spire['Ic_ergscm2um']*(factor**2)
    
    # Optical Depths
    model_taudust = pd.read_csv(model_path+'/'+model_name+'_.taudust', delim_whitespace= True)
    model_taudust.columns = ['lambda_um', 'taudust']
    
    # =============================================================================
    # Luminosity 10-1200um
    # =============================================================================
        
    luminosity_iso = 0
    for i, row in model_iso.iterrows():
        if i == 0:
            continue
        else:
            luminosity_iso += (model_iso['Ic_ergscm2um'][i]*
                              (model_iso['lambda_um'][i]-model_iso['lambda_um'][i-1])*
                              4.*3.14159*(3.5e6*3.086e18)**2)/(3.8e33)
    luminosity_spire = 0
    for i, row in model_spire.iterrows():
        if i == 0:
            continue
        else:
            luminosity_spire += (model_spire['Ic_ergscm2um'][i]*
                              (model_spire['lambda_um'][i]-model_spire['lambda_um'][i-1])*
                              4.*3.14159*(3.5e6*3.086e18)**2)/(3.8e33)
        
    luminosity_total = luminosity_iso+luminosity_spire # lsun Para un tamaño de 0.1" si factor = 0.7
    luminosity_factor = luminosity_total/luminosity_obs
    
    
    orig_luminosity_factor = (luminosity_total/(factor**2))/luminosity_obs
    new_factor = np.sqrt(1./ orig_luminosity_factor)
    new_luminosity_total =(luminosity_total/(factor**2))*new_factor # lsun
    new_model_size = 1. * new_factor # pc

    
    if plot_SED == True:
        
        # observational SED from perez-beaupuits 2018
        obs_SED = pd.read_csv('/Users/frico/Documents/data/NGC253_H3O+/NGC253_SED.txt', delim_whitespace = True, comment='#')


        SED_model_pars =pd.read_csv('/Users/frico/Documents/data/NGC253_H3O+/NGC253_SED_fit_parameters.txt', delim_whitespace = True, comment='#')

        
        modeled_comp_cold =[]
        modeled_comp_warm =[]
        modeled_comp_hot =[]
        modeled_comp_all = []
        wave_mod = range(10, 1200+1)
        D = 3.52
        for i, wave in enumerate(wave_mod):
            nu = (wave*u.um).to(u.GHz, equivalencies=u.spectral())
            df_cold = SED_model_pars.loc[SED_model_pars['Component'] == 'Cold']
            df_warm = SED_model_pars.loc[SED_model_pars['Component'] == 'Warm']
            df_hot = SED_model_pars.loc[SED_model_pars['Component'] == 'Hot']
            
            
            modeled_cold = utiles.SED_model(nu, df_cold['Tdust'][0], df_cold['Mdust_1e6'][0]*1e6, df_cold['filling_fact'][0], D)
            modeled_warm = utiles.SED_model(nu, df_warm['Tdust'][1], df_warm['Mdust_1e6'][1]*1e6, df_warm['filling_fact'][1], D)
            modeled_hot = utiles.SED_model(nu, df_hot['Tdust'][2], df_hot['Mdust_1e6'][2]*1e6, df_hot['filling_fact'][2], D)
            
            modeled_comp_cold.append(modeled_cold.value)
            modeled_comp_warm.append(modeled_warm.value)
            modeled_comp_hot.append(modeled_hot.value)
            modeled_comp_all.append(modeled_cold.value+modeled_warm.value+modeled_hot.value)
            
        m = 2
        n = 1
        fig_SED = plt.figure()
        gs1 = gridspec.GridSpec(m, n)    
        gs1.update(wspace = 0.0, hspace=0.0, top=0.95, bottom = 0.05)
        # Generating specified number of axis
        axis = []
        for i in range(m*n):
            axis.append(fig_SED.add_subplot(gs1[i]))
        for j in range(m*n):
            # SED
            if j == 0:
                # SED
                axis[j].plot(model_iso['lambda_um'], model_iso['Ic_Jy'], linewidth=linewidth, color='k')
                axis[j].plot(model_spire['lambda_um'], model_spire['Ic_Jy'], linewidth=linewidth, color='k')
                # Observed SED
                axis[j].errorbar(obs_SED['Wavlength'], obs_SED['Obs_flux'],
                                yerr=obs_SED['Obs_flux_err'], fmt='.',
                                elinewidth=0.7, capsize=0.5,
                                linewidth=linewidth, color='k')
                # Modeled SED
                axis[j].plot(wave_mod, modeled_comp_cold, linewidth=linewidth, color='blue')
                axis[j].plot(wave_mod, modeled_comp_warm, linewidth=linewidth, color='green')
                axis[j].plot(wave_mod, modeled_comp_hot, linewidth=linewidth, color='red')
                axis[j].plot(wave_mod, modeled_comp_all, linewidth=linewidth, color='k', linestyle='--')
                
            # Optical depth
            elif j==1:
                axis[j].plot(model_taudust['lambda_um'], model_taudust['taudust'], linewidth=linewidth, color='k')
                # Fit opt depth from SED
                axis[j].errorbar(obs_SED['Wavlength'], obs_SED['opt_depth'],
                                yerr=obs_SED['Obs_flux_err'], fmt='.',
                                elinewidth=0.7, capsize=0.5,
                                linewidth=linewidth, color='k')
        for i, ax in enumerate(axis):
            majorLocator = MultipleLocator(20)
            ax.yaxis.set_major_locator(majorLocator)
            # Solo ponemos titulo a los ejes y si estan a la izquierda del todo
            ax.minorticks_on()
            #ax.tick_params(axis='both', which='minor', direction='in', top=True, right=True)
            ax.tick_params(axis='both', which='major', labelsize=labelsize, length=5, direction='in')
            ax.tick_params(axis='both', which='minor', labelsize=labelsize, direction='in')
            ax.xaxis.set_tick_params(top ='on', which='both')
            ax.yaxis.set_tick_params(right='on', which='both', labelright='off')
            ax.set_xlim([10,1200])
            ax.set_yscale('log')
            ax.set_xscale('log')
        
        axis[0].set_ylim([1e-2,1e4])
        axis[0].set_xticklabels([])
        axis[0].set_ylabel(r'Flux density (Jy)', fontsize=14)
        axis[1].set_ylabel(r'Optical depth $\tau$', fontsize=14)
        axis[1].set_xlabel(r'Wavelength ($\mu$m)', fontsize=14)
        
        axis[0].text(0.1, 0.9,
                            'L='+ '% 1.2E' % luminosity_total + r' L$_\odot$', ha='left', va='center',
                            transform=axis[j].transAxes,
                            rotation='horizontal', fontsize=8)
        axis[0].text(0.5, 0.9,
                            'ss='+ '% 1.3f' % model_size + r' arcsec', ha='left', va='center',
                            transform=axis[j].transAxes,
                            rotation='horizontal', fontsize=8)
        
        #fig_SED.savefig(figoutdir+'/'+model_name+'_f'+'%1.1f' % factor +'_SED.'+fig_format, bbox_inches='tight', transparent=True, dpi=400, format=fig_format)
        #plt.close()
        fig_SED.savefig(out_dir_sed+'/'+model_name+'_f'+'%1.1f' % factor +'_T'+'%1.f' % Temp+'_SED.'+fig_format, bbox_inches='tight', transparent=True, dpi=400, format=fig_format)
        plt.close(fig_SED)
            
    
    
    return luminosity_total, peak_v0_39_38, peak_v0_24_23, peak_v7_39_38, peak_v7_24_23, fig_PROF, fig_SED


def plot_nh_ch2_fun(t_grad, xaxis, yaxis, 
                    profile_list, start_dens_list, T_list,
                    Xline_list_or_Nline_list, luminosity_list,
                    version, hc3n_path, cont_plot, 
                    plot_by_X=True, plot_sep=False, plot_taus=False, plot_intens = False):
    """
    xaxis
        central_nh
        central_T
    yaxis
        ratios -> line ratios
        taus -> line opacities
        intens -> line intensities
    """
    # Reading obs data
    hc3n_prop = read_hc3n_obs(hc3n_path)
    # Default Pars
    profile_names= ['const', 'r', 'r^1.5', 'r^2']
    Xdust_list = [1E-2]
    marker_list = ['.', 's', '^', 'v', '*', 'd']
    dens_colors = ['#62a1db', '#55d180', '#e7d87d', '#dd9f40', '#b01111', 'pink', '#581845', 'k'] 
    nline_colors = dens_colors
    linestyles_list = ['-', '--', '-.', ':']
    plot_by_N = plot_sep
    # Plotting by density
    if xaxis=='central_nh':
        #plot_by_nh = True
        xaxis_name = 'nH'
        xaxis_label= 'r$n_H$'
    elif xaxis== 'central_T':
        xaxis_name = 'T'
    
    
    hc3n_prop['Nline'] = 10**hc3n_prop['logN_vib']
    
    
    

    


    # Making figure per density profile type
    """
    Plotting ratios vs density
    """
    
    # Selecting observed columns to plot
    if yaxis == 'ratios':
        ratios_list = ['ratio_24_23', 'ratio_39_38', 'ratio_v0', 'ratio_v7']
        #plot_ratios = True
        y_axis_obs = ratios_list
        hc_pos = 3e6 # Generic density for the obs data
        set_ylim = 20 # ylim
        pane_names = ['v0_24-23/v7_24-23', 'v0_39-38/v7_39-38', 'v0_39-38/v0_24-23', 'v7_39-38/v7_24-23']
    elif yaxis == 'taus':
        #plot_taus = True
        y_axis_obs = ['tau_v0_24_23', 'tau_v7_24_23', 'tau_v0_39_38', 'tau_v7_39_38']
        pane_names = ['v0_24_23', 'v7_24_23', 'v0_39_38', 'v7_39_38']
        set_ylim = 15 # ylim
        hc_pos = 0
        
    elif yaxis == 'intens':
        #plot_intens = True
        y_axis_obs = ['peak_v0_24_23', 'peak_v7_24_23', 'peak_v0_39_38', 'peak_v7_39_38']
        #y_axis_obs_mod = ['peak_v0_24_23', 'peak_v7_24_23', 'peak_v0_39_38', 'peak_v7_39_38']
        pane_names = ['v0_24_23', 'v7_24_23', 'v0_39_38', 'v7_39_38']
        hc_pos = 0
        set_ylim = 0 # No ylim set
        
    if xaxis == 'central_nh':
        hc_pos = 3e6 
        
    z_axis_obs = ['peak_v0_24_23', 'peak_v7_24_23', 'peak_v0_39_38', 'peak_v7_39_38']
    
    
                       
                                
    constant_T2 = True                  
    if constant_T2 == True:
        # Results path
        ttt= 'tcte_T'
        res_path = '/Users/frico/Documents/Ed/modelos_hc3n/v'+version+'_'+ttt+'/'
        for p, profile in enumerate(profile_list):
            print '\n\nDensprofile:\t'+profile_names[p]
            # Directories
            directory_list = glob(res_path+'r_'+str(profile)+'/*/')
            print res_path+'r_'+str(profile)+'/*/'
            temp_list = [100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600]
            if res_path+'r_'+str(profile)+'/Figures/' in directory_list: directory_list.remove(res_path+'r_'+str(profile)+'/Figures/')
            if res_path+'r_'+str(profile)+'/fig_ok/' in directory_list: directory_list.remove(res_path+'r_'+str(profile)+'/fig_ok/')
            
            Ttogether = True
            if Ttogether == True:
                model_list = []
                for a, direct in enumerate(directory_list):
                    #if 1 == 1:
                    #direct = directory_list[0]
                    # Reading models summary
                    print '\t'+direct
                    hc3n_model = hc3n_model_ratio_cal(direct, '035')
                    hc3n_model.drop_duplicates(inplace=True)
                    
                    model_list.append(hc3n_model)
                
                    
        
                hc3n_model = pd.concat(model_list, ignore_index=False)
                #print hc3n_model['central_T'].tolist()
                # Appending all models   
                # Plotting models
                #dust_abu = [5.4E19, 5.4E20, 5.4E21, 5.4E22, 5.4E23,
                #            1.0E21] 
                #if profile == 0.0:
                #dust_abu = [2.7E18, 2.7E19, 2.7E20, 3.57E20, 5.6E20, 8.1E20, 2.7E21, 5.56E21, 9.0E21]
                 
                
                
                dust_abu = [9.2E21]#[4.5E20, 9.2E20, 4.5E21, 9.2E21, 4.5E22, 9E22] # 9.2E18, 4.5E19, 9.2E19, 
                slicer = 'Ndust'
                
                hc3n_model['cond_test'] = (hc3n_model['central_T'] == 250 ) & (hc3n_model[slicer]/dust_abu[0] > 0.70) & (hc3n_model[slicer]/dust_abu[0] < 1.30)
                hc3n_model_test = hc3n_model[hc3n_model['cond_test'] == True]
                hc3n_model_test.drop_duplicates(inplace=True)
                hc3n_model_test[['ratio_v7', 'central_nh', 'central_T']] 
                
                hc3n_model.to_csv(res_path+'r_'+str(profile)+'/models_summary_035_ALL.txt', header=True, sep='\t', index=False, na_rep=-1)
                
                xvar = 'ratio_v7'
                yvar = 'ratio_24_23'#'logN'
                zvar = 'central_T'
                
                hc3n_prop['central_T'] = hc3n_prop['Tvib']
                
                
                for da, dabu in enumerate(dust_abu):
                    print da
                    #da = 2
                    #hc3n_model = hc3n_model_ratio_cal(direct, '035')
                    
                    # Selecting a same Ndust
                    hc3n_model['cond_'+str(da)] = (hc3n_model[slicer]/dust_abu[da] > 0.70) & (hc3n_model[slicer]/dust_abu[da] < 1.30)
                    ncol_dust_cont = hc3n_model[hc3n_model['cond_'+str(da)] == True]
                    ncol_dust = hc3n_model[hc3n_model['cond_'+str(da)] == True]
                    
                    
                    ncol_dust['logN'] = np.log10(ncol_dust['Nline'])
                    Ndust = ncol_dust['Ndust'].mean()
                    
                    ncol_dust_cont['logN'] = np.log10(ncol_dust_cont['Nline'])
                    Ndust = ncol_dust_cont['Ndust'].mean()
                    
                    
                    
                    name_suf = ''
                    if zvar != 'central_T':
                        temp = 250
                        # Selecting a same nH denistiy
                        ncol_dust_cont['cond_T_'+str(da)] = (ncol_dust_cont['central_T']/temp > 0.9) & (ncol_dust_cont['central_T']/temp < 1.10)
                        ncol_dust_cont = ncol_dust_cont[ncol_dust_cont['cond_T_'+str(da)] == True]
                        name_suf = name_suf+'_T_'+'%1.f' % temp
                    else:
                        name_suf = name_suf+'_T_all'
                    
#                    if zvar != 'central_nh':
#                        nh_density = 3.0e6
#                        # Selecting a same nH denistiy
#                        ncol['cond_nh_'+str(da)] = (ncol['central_nh']/nh_density > 0.9) & (ncol['central_nh']/nh_density < 1.10)
#                        ncol = ncol[ncol['cond_nh_'+str(da)] == True]
#                        name_suf = name_suf+'_nH_'+'%1.1e' % nh_density
#                    else:
#                        name_suf = name_suf+'_nH_all'
                        
                    if yvar != 'logN':
                        hc3n_col = 2.8E16
                        # Selecting a same NHC3N column
                        ncol_dust_cont['cond_Nline_'+str(da)] = (ncol_dust_cont['Nline']/hc3n_col > 0.90) & (ncol_dust_cont['Nline']/hc3n_col < 1.10)
                        ncol_dust_cont = ncol_dust_cont[ncol_dust_cont['cond_Nline_'+str(da)] == True]
                        name_suf = name_suf+'_NHC3N_'+'%1.1e' % hc3n_col
                    else:
                        name_suf = name_suf+'_NHC3N_all'
                    
                    # Dropping dens above 3e7
                    #if profile == 0.0:
                    #    ncol = ncol[ncol['central_nh'] < 3e7]
                    #ncol = ncol[ncol['central_nh'] < 3e7]
                    #for da, dabu in enumerate(dust_abu):
                    #da = 1
                    
                    #ncol = hc3n_model[(hc3n_model[slicer]/dust_abu[da] > 0.70) & (hc3n_model[slicer]/dust_abu[da] < 1.30)]
                    
                    # Selecting less temperatures
                    ncol_dust_cont['cond_temp_'+str(da)] = (ncol_dust_cont['central_T']<= 550)
                    ncol_dust_cont = ncol_dust_cont[ncol_dust_cont['cond_temp_'+str(da)] == True]
                    
                    ncol_dust_cont = ncol_dust_cont.sort_values(by=[xaxis], ascending = True)
                    ncol_dust = ncol_dust.sort_values(by=[xaxis], ascending = True)
                    if len(ncol_dust_cont)!= 0:
                        print dabu
                        
                        #ncol['logN'] = np.log10(ncol['Nline'])
                        #hc3n_prop['central_T'] = hc3n_prop['Tvib']
                        
                        x = ncol_dust_cont[xvar]
                        y = ncol_dust_cont[yvar]
                        z = ncol_dust_cont[zvar]
                        zdens = ncol_dust_cont['central_nh']
                        
                        # Check if all elements in z are the same -> cant make contours
                        no_zrange = utiles.checkEqual(ncol_dust_cont[zvar])
                        lev = np.linspace(z.min(), z.max(), 20)
                        lev = np.arange(100, 600+50, 50)
                    
                        if no_zrange != True:
                            
                            plot_good_models = True
                            if plot_good_models == True:
                        
                                fig3 = plt.figure()
                                ax = plt.axes()
                
                                
                                
                                
                                
                                ncol_dust_cont['Ltot_14'] = ncol_dust_cont['Ltot']*0.3599/0.35
                                ncol_dust_cont['Ltot_13'] = ncol_dust_cont['Ltot']*0.2317/0.35
                                ncol_dust_cont['Ltot_2']  = ncol_dust_cont['Ltot']*0.1806/0.35
                                ncol_dust_cont['Ltot_1']  = ncol_dust_cont['Ltot']*0.1366/0.35
                                
                                
                                #'central_nh'
                                
                                
                                # Subsetting model
                                Ndust = ncol_dust_cont['Ndust'].mean()
            
                                xi, yi = np.linspace(x.min(), x.max(), 1000), np.linspace(y.min(), y.max(), 1000)
                                xi, yi = np.meshgrid(xi, yi)
                                rbf = scipy.interpolate.Rbf(x, y, z, function='linear', smooth=0.5)
                                zi = rbf(xi, yi)
                                
                                xi2, yi2 = np.linspace(x.min(), x.max(), 1000), np.linspace(y.min(), y.max(), 1000)
                                xi2, yi2 = np.meshgrid(xi2, yi2)
                                rbf2 = scipy.interpolate.Rbf(x, y, zdens, function='linear')
                                zi2 = rbf2(xi2, yi2)
                                
                                color_list1 = ["#E1EEE0","#BADEDC","#8BC7C8", "#5D8EAC",  "#3B4461", "#262637"]
                                color_list2 = ["#71c7ec" , "#1ebbd7", "#189ad3", "#107dac", "#005073"]
                                colormap = matplotlib.colors.LinearSegmentedColormap.from_list("", color_list2)
            
                                #print  ncol[zvar]
                                
                                
                                #color_list3 = ["#71c7ec" , "#1ebbd7", "#189ad3", "#107dac", "#005073"]
                                #colormap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", color_list3)
                                cs = ax.contourf(xi, yi, zi,  cmap="Reds", origin='lower', levels=lev)
                                cbar = plt.colorbar(cs)
                                
                                #contours = plt.contour(xi2, yi2, zi2, 100, colors='black')
                                #plt.clabel(contours, inline=True, fontsize=8)
                                
                                
                                
                                dens_plot = [2.5e6, 7.5e6, 2.5e7, 7.5e7]
                                unique_dens = ncol_dust_cont.central_nh.unique()
                                color_ndens_lim = 0
                                aa_p = np.linspace(1e6, 1e7, 10).tolist()
                                #bb_p = np.linspace(1e7, 1e8, 10).tolist()
                                bb_p = [1e8]
                                nppplot = aa_p+bb_p#+dens_plot
                                nnnp2 = [1e6,  5e6, 1e7, 1e8]
                                for n, ndens in enumerate(nppplot):#(unique_dens):
                                    color_ndens_lim += 1
                                    
                                #73E5CF
                                cmap_blues = matplotlib.cm.get_cmap('winter_r')#('Blues')
                                color_list_b = np.linspace(0.1, 1,color_ndens_lim )#len(unique_dens))
                                custom_lines_dens = []
                                custom_lines_txt_dens = []
                                ncol_dust_cont.drop_duplicates(inplace=True)
                                s_unique_dens = np.sort(unique_dens)
                                print s_unique_dens
                                s_nppplot = np.sort(nppplot)
                                for n, ndens in enumerate(s_nppplot):#(s_unique_dens):
                                    print 'aaaaaaaa'
                                    if ndens in s_nppplot:
                                        if ndens in nnnp2:
                                            linestyle='-'
                                            lw=0.6
                                            color='0.3'
                                        else:
                                            linestyle='--'
                                            lw=0.6
                                            color='grey'
                                        ncol_dust_cont['cond_nH_'+'%1.1E' % ndens] = (ncol_dust_cont['central_nh']/ndens > 0.98) & (ncol_dust_cont['central_nh']/ndens < 1.02) & (ncol_dust_cont['central_nh'] <= 1e8)
                                        ncol_dens = ncol_dust_cont[ncol_dust_cont['cond_nH_'+'%1.1E' % ndens] == True]
                                        line = ax.plot(ncol_dens[xvar], ncol_dens[yvar], marker='None', linestyle=linestyle, color=color, markersize = 0.8, linewidth= lw, zorder=3)#=cmap_blues(color_list_b[n])
                                        label_text = '%1.1E' % ndens
                                    if ndens in nnnp2:
                                        if ndens == 1e6:
                                            xpos=1.06
                                            ypos=1.9
                                            label_text = r'$1 \times 10^6$'
                                        elif ndens == 5e6:
                                            xpos=2.5625
                                            ypos=2.1
                                            label_text = r'$5 \times 10^6$'
                                        elif ndens == 1e7:
                                            xpos=2.5625
                                            ypos=2.7
                                            label_text = r'$1 \times 10^7$'
                                        elif ndens == 5e7:
                                            xpos=2.5625
                                            ypos=3.04
                                            label_text = r'$5 \times 10^7$'
                                        elif ndens == 1e8:
                                            xpos=2.5625
                                            ypos=3.25
                                            label_text = r'$1 \times 10^8$'
                                        labelLines(ncol_dens[xvar], ncol_dens[yvar], ax, xpos=xpos, ypos=ypos, color='0.3', label=label_text, align=True)
                                        #label_line(ncol_dens[xvar], ncol_dens[yvar], ax, label_text, near_x=ncol_dens[xvar].tolist()[-4])
                                    #if ndens in nnnp2:
                                        #labelLines(ax.get_lines(),zorder=2.5)
                                        #label_line(xdata=ncol_dens[xvar], ydata=ncol_dens[yvar], ax=ax, label_text=label_text, near_i=1, near_x=None, near_y=None, rotation_offset=0, offset=(0,0))
                                    
                                    #aba = plt.gca().get_lines()
                                    #print aba
                                    
                                    #inline_legend(plt.gca().get_lines(), n_markers=3)
                                        
                                
                                if xvar == 'central_T':
                                    xvar_obs = 'Tvib'
                                else:
                                    xvar_obs = xvar
                                
                                if yvar == 'logN':
                                    yvar_obs = 'logN_vib'
                                    y_range = [14, 17]
                                    ax.set_ylim([14,17])
                                else:
                                    yvar_obs = yvar
                                    
                                    
                                if xvar == 'ratio_v7':
                                    x_range = [0,4]
                                    ax.set_xlim([0,4])
                                elif xvar == 'ratio_24_23':
                                    x_range = [1, 6]
                                    ax.set_xlim([1,6])
                                elif xvar == 'ratio_39_38':
                                    x_range = [0, 15]
                                    ax.set_xlim([0,15])
                                elif xvar == 'ratio_v0':
                                    x_range = [0.5,2.5]
                                    ax.set_xlim([0.5,2.5])
                                    
                                if yvar == 'ratio_24_23':
                                    y_range = [1,6]
                                    ax.set_ylim([1,6])
                                elif yvar == 'ratio_39_38':
                                    y_range = [0,15]
                                    ax.set_ylim([0,15])
                                    
                                # Obs
                                hc3n_prop['Ind_ok'].astype(float)
                                hc3n_prop = hc3n_prop.sort_values(by=['Ind_ok'], ascending = True)
                                for o, orow in hc3n_prop.iterrows():
                                    #if orow['peak_v0_39_38'] != np.nan:
                                    
                                    if not pd.isnull(orow[xvar_obs]) and not pd.isnull(orow[yvar_obs]):
                                        #index = .astype(float)
                                        #if orow['Ind_ok'] != 12:
                                        if orow['Ind_ok'] not in [6,7,9,10,12]:#[9, 10, 11, 12]:
                                            ax.scatter(orow[xvar_obs], orow[yvar_obs], marker='.', c='k', s=2, zorder =10)
                                            #ax.errorbar(orow[xvar_obs], orow[yvar_obs], xerr=orow[xvar_obs+'_err'], yerr=orow[yvar_obs+'_err'], fmt='.', mec='k', mfc='k', ecolor='k', elinewidth=0.5)
                                            ax.text(orow[xvar_obs], orow[yvar_obs],
                                                    str(int(orow['Ind_ok'])), ha='left', va='bottom', color='k',
                                                    #transform=axis[j].transAxes,
                                                    rotation='horizontal', fontsize=9, zorder =11)
                                        
                                        
                                #ax.scatter(hc3n_prop[xvar_obs], hc3n_prop[yvar_obs], s=4, edgecolors='k', facecolors='k')
                                
                                out_dir = res_path+'r_'+str(profile)+'/'+'Figures/'
                                if not os.path.exists(out_dir):
                                    os.makedirs(out_dir)
                                    
                                if xvar == 'ratio_v7':
                                    xlabel = 'v7=1 (39-38) / v7=1 (24-23)'
                                else:
                                    xlabel = xvar
                                if yvar == 'ratio_24_23':
                                    ylabel = 'v=0 (24-23) / v7=1 (24-23)'
                                elif yvar == 'ratio_39_38':
                                    ylabel = 'v=0 (39-38) / v7=1 (39-38)'    
                                
                                else:
                                    ylabel = yvar
                                if zvar == 'central_T':
                                    zlabel = r'T$_{\rm{dust}}$ (K)'
                                else:
                                    zlabel = zvar
                                    
                                ax.set_xlabel(xlabel)
                                ax.set_ylabel(ylabel)
                                cbar.ax.set_ylabel(zlabel)
                                
                                ax.set_xlim([1,3])
                                ax.set_ylim([1.75,5])
                                    
                                Td = ncol_dust_cont['central_T'].mean()
                                
                                custom_lines_dens = []
                                custom_lines_txt_dens = []
                                for n, ndens in enumerate(s_nppplot):
                                    if ndens in nnnp2:#dens_plot:
                                    #if 1 == 1:
                                        linestyle='-'
                                        custom_lines_dens.append(Line2D([0], [0], color='0.4', lw=0.5, linestyle=linestyle))#, color=cmap_blues(color_list_b[n]), lw=0.5, linestyle=linestyle))
                                        custom_lines_txt_dens.append('nH='+'%1.1E' % ndens)
                                #legend_dens = plt.legend(custom_lines_dens, custom_lines_txt_dens, prop={'size': 4}, loc=1)
                                #ax.add_artist(legend_dens)
                                
                                ax.tick_params(axis='both', which='both', direction='in')
                                ax.minorticks_on()
                                ax.xaxis.set_tick_params(which='both', top ='on')
                                ax.yaxis.set_tick_params(which='both', right='on', labelright='off')
                                
                                fig3_name = 'hc3n_r'+str(profile)+'_'+xvar+'VS'+yvar+'VS'+zvar+'_Nd'+'%1.1E'  % Ndust+name_suf+'_final_newint_bold'
                                fig3.savefig(out_dir+fig3_name+'_all.'+fig_format, bbox_inches='tight', transparent=True, dpi=400)
                                fig3.savefig(out_dir+fig3_name+'_all.'+'eps', bbox_inches='tight', transparent=True, dpi=400)
                                plt.close(fig3)
                            
                            plot_all_models = False
                            if plot_all_models == True:
                                fig3 = plt.figure()
                                ax = plt.axes()
                
                                
                                
                                
                                
                                ncol_dust_cont['Ltot_14'] = ncol_dust_cont['Ltot']*0.3599/0.35
                                ncol_dust_cont['Ltot_13'] = ncol_dust_cont['Ltot']*0.2317/0.35
                                ncol_dust_cont['Ltot_2']  = ncol_dust_cont['Ltot']*0.1806/0.35
                                ncol_dust_cont['Ltot_1']  = ncol_dust_cont['Ltot']*0.1366/0.35
                                
                                
                                #'central_nh'
                                
                                
                                # Subsetting model
                                Ndust = ncol_dust_cont['Ndust'].mean()
            
                                xi, yi = np.linspace(x.min(), x.max(), 1000), np.linspace(y.min(), y.max(), 1000)
                                xi, yi = np.meshgrid(xi, yi)
                                rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
                                zi = rbf(xi, yi)
                                
                                xi2, yi2 = np.linspace(x.min(), x.max(), 1000), np.linspace(y.min(), y.max(), 1000)
                                xi2, yi2 = np.meshgrid(xi2, yi2)
                                rbf2 = scipy.interpolate.Rbf(x, y, zdens, function='linear')
                                zi2 = rbf2(xi2, yi2)
                                
                                color_list1 = ["#E1EEE0","#BADEDC","#8BC7C8", "#5D8EAC",  "#3B4461", "#262637"]
                                color_list2 = ["#71c7ec" , "#1ebbd7", "#189ad3", "#107dac", "#005073"]
                                colormap = matplotlib.colors.LinearSegmentedColormap.from_list("", color_list2)
            
                                #print  ncol[zvar]
                                
                                
                                
                                cs = ax.contourf(xi, yi, zi,  cmap="Reds", origin='lower', levels=lev)
                                cbar = plt.colorbar(cs)
                                
                                #contours = plt.contour(xi2, yi2, zi2, 100, colors='black')
                                #plt.clabel(contours, inline=True, fontsize=8)
                                
                                
                                
                                dens_plot = [1e6, 2.5e6, 5e6, 7.5e6, 1e7, 2.5e7,5e7,  7.5e7, 1e8]
                                unique_dens = ncol_dust_cont.central_nh.unique()
                                color_ndens_lim = 0
                                aa_p = np.linspace(1e6, 1e7, 10).tolist()
                                #bb_p = np.linspace(1e7, 1e8, 10).tolist()
                                bb_p = [5e7,1e8]
                                nppplot = aa_p+bb_p#+dens_plot
                                nnnp2 = [1e6,  5e6, 1e7, 5e7, 1e8]
                                s_unique_dens = np.sort(unique_dens)
                                for n, ndens in enumerate(s_unique_dens):#(unique_dens):
                                    color_ndens_lim += 1
                                cmap_blues = matplotlib.cm.get_cmap('winter_r')#('Blues')
                                color_list_b = np.linspace(0.1, 1,color_ndens_lim )#len(unique_dens))
                                custom_lines_dens = []
                                custom_lines_txt_dens = []
                                ncol_dust_cont.drop_duplicates(inplace=True)
                                
                                s_nppplot = np.sort(nppplot)
                                for n, ndens in enumerate(s_unique_dens):
                                    print 'aaaaaaaa'
                                    if ndens in dens_plot:
                                        linestyle='-'
                                        lw=0.6
                                        color='k'
                                    else:
                                        linestyle='--'
                                        lw=0.6
                                        color='gray'
                                    ncol_dust_cont['cond_nH_'+'%1.1E' % ndens] = (ncol_dust_cont['central_nh']/ndens > 0.98) & (ncol_dust_cont['central_nh']/ndens < 1.02) & (ncol_dust_cont['central_nh'] <= 1e8)
                                    ncol_dens = ncol_dust_cont[ncol_dust_cont['cond_nH_'+'%1.1E' % ndens] == True]
                                    line = ax.plot(ncol_dens[xvar], ncol_dens[yvar], marker='None', linestyle=linestyle, color=cmap_blues(color_list_b[n]), markersize = 0.8, linewidth= lw, zorder=3)
                                    label_text = '%1.1E' % ndens
    
                                    #aba = plt.gca().get_lines()
                                    #print aba
                                    
                                    #inline_legend(plt.gca().get_lines(), n_markers=3)
                                        
                                
                                if xvar == 'central_T':
                                    xvar_obs = 'Tvib'
                                else:
                                    xvar_obs = xvar
                                
                                if yvar == 'logN':
                                    yvar_obs = 'logN_vib'
                                    y_range = [14, 17]
                                    ax.set_ylim([14,17])
                                else:
                                    yvar_obs = yvar
                                    
                                    
                                if xvar == 'ratio_v7':
                                    x_range = [0,4]
                                    ax.set_xlim([0,4])
                                elif xvar == 'ratio_24_23':
                                    x_range = [1, 6]
                                    ax.set_xlim([1,6])
                                elif xvar == 'ratio_39_38':
                                    x_range = [0, 15]
                                    ax.set_xlim([0,15])
                                elif xvar == 'ratio_v0':
                                    x_range = [0.5,2.5]
                                    ax.set_xlim([0.5,2.5])
                                    
                                if yvar == 'ratio_24_23':
                                    y_range = [1,6]
                                    ax.set_ylim([1,6])
                                elif yvar == 'ratio_39_38':
                                    y_range = [0,15]
                                    ax.set_ylim([0,15])
                                    
                                # Obs
                                hc3n_prop['Ind_ok'].astype(float)
                                hc3n_prop = hc3n_prop.sort_values(by=['Ind_ok'], ascending = True)
                                for o, orow in hc3n_prop.iterrows():
                                    #if orow['peak_v0_39_38'] != np.nan:
                                    
                                    if not pd.isnull(orow[xvar_obs]) and not pd.isnull(orow[yvar_obs]):
                                        #index = .astype(float)
                                        #if orow['Ind_ok'] != 12:
                                        if orow['Ind_ok'] not in [6,7, 9, 10,12]:#[9, 10, 11, 12]:
                                            ax.scatter(orow[xvar_obs], orow[yvar_obs], marker='.', c='k', s=2, zorder =10)
                                            #ax.errorbar(orow[xvar_obs], orow[yvar_obs], xerr=orow[xvar_obs+'_err'], yerr=orow[yvar_obs+'_err'], fmt='.', mec='k', mfc='k', ecolor='k', elinewidth=0.5)
                                            ax.text(orow[xvar_obs], orow[yvar_obs],
                                                    str(int(orow['Ind_ok'])), ha='left', va='bottom', color='k',
                                                    #transform=axis[j].transAxes,
                                                    rotation='horizontal', fontsize=7, zorder =11)
                                        
                                        
                                #ax.scatter(hc3n_prop[xvar_obs], hc3n_prop[yvar_obs], s=4, edgecolors='k', facecolors='k')
                                
                                out_dir = res_path+'r_'+str(profile)+'/'+'Figures/'
                                if not os.path.exists(out_dir):
                                    os.makedirs(out_dir)
                                    
                                if xvar == 'ratio_v7':
                                    xlabel = 'v7=1 (39-38) / v7=1 (24-23)'
                                else:
                                    xlabel = xvar
                                if yvar == 'ratio_24_23':
                                    ylabel = 'v=0 (24-23) / v7=1 (24-23)'
                                elif yvar == 'ratio_39_38':
                                    ylabel = 'v=0 (39-38) / v7=1 (39-38)'    
                                
                                else:
                                    ylabel = yvar
                                if zvar == 'central_T':
                                    zlabel = 'T (K)'
                                else:
                                    zlabel = zvar
                                    
                                ax.set_xlabel(xlabel)
                                ax.set_ylabel(ylabel)
                                cbar.ax.set_ylabel(zlabel)
                                
                                ax.set_xlim([1,3])
                                ax.set_ylim([1.75,5])
                                    
                                Td = ncol_dust_cont['central_T'].mean()
                                
                                custom_lines_dens = []
                                custom_lines_txt_dens = []
                                for n, ndens in enumerate(s_unique_dens):
                                    if ndens in dens_plot:
                                    #if 1 == 1:
                                        linestyle='-'
                                        custom_lines_dens.append(Line2D([0], [0], color=cmap_blues(color_list_b[n]), lw=0.5, linestyle=linestyle))
                                        custom_lines_txt_dens.append('nH='+'%1.1E' % ndens)
                                legend_dens = plt.legend(custom_lines_dens, custom_lines_txt_dens, prop={'size': 4}, loc=1)
                                ax.add_artist(legend_dens)
                                
                                ax.tick_params(axis='both', which='both', direction='in')
                                ax.minorticks_on()
                                ax.xaxis.set_tick_params(which='both', top ='on')
                                ax.yaxis.set_tick_params(which='both', right='on', labelright='off')
                                
                                fig3_name = 'hc3n_r'+str(profile)+'_'+xvar+'VS'+yvar+'VS'+zvar+'_Nd'+'%1.1E'  % Ndust+name_suf+'_TODOS'
                                fig3.savefig(out_dir+fig3_name+'_all2.'+fig_format, bbox_inches='tight', transparent=True, dpi=300)
                                fig3.savefig(out_dir+fig3_name+'_all2.'+'eps', bbox_inches='tight', transparent=True, dpi=300)
                                plt.close(fig3)
                            
                            # Selecting models
                            # HC1, 2, 3, 4, 5, 8, 13, 14
                            # models have 1pc radies
                            D = 3.5 # Mpc
                            
                            HC1_nh_cond = [3.0e6]
                            HC1_T_cond = [250, 300]
                            
                            HC1_size_fact = u_conversion.lin_size(3.5, 0.0202).to(u.pc)/2. # pc
                            
                            ncol_dust_cont['cond_HC1_T250'] = (ncol_dust_cont['central_nh'] == HC1_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC1_T_cond[0])
                            ncol_dust_cont['cond_HC1_T300'] = (ncol_dust_cont['central_nh'] == HC1_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC1_T_cond[1])
                            HC1_pos_mod1 = ncol_dust_cont[ncol_dust_cont['cond_HC1_T250'] == True] 
                            HC1_pos_mod2 = ncol_dust_cont[ncol_dust_cont['cond_HC1_T300'] == True] 
                            HC1_pos_mods = pd.concat([HC1_pos_mod1, HC1_pos_mod2], ignore_index=True)
                            HC1_pos_mods['HC_1_Ltot']  = HC1_pos_mods['Ltot']*(HC1_size_fact**2)/(0.35**2)
                            HC1_pos_mods['HC_factor'] = HC1_size_fact
                            
                            HC1_pos_mods['HC_1_peak_v0_39_38_fok']  = HC1_pos_mods['peak_v0_39_38']*(HC1_size_fact**2)/(0.35**2)
                            HC1_pos_mods['HC_1_peak_v0_24_23_fok']  = HC1_pos_mods['peak_v0_24_23']*(HC1_size_fact**2)/(0.35**2)
                            HC1_pos_mods['HC_1_peak_v7_39_38_fok']  = HC1_pos_mods['peak_v7_39_38']*(HC1_size_fact**2)/(0.35**2)
                            HC1_pos_mods['HC_1_peak_v7_24_23_fok']  = HC1_pos_mods['peak_v7_24_23']*(HC1_size_fact**2)/(0.35**2)
                            HC1_pos_mods['HC_1_size_pc'] = 2. * HC1_size_fact
                            HC1_pos_mods['HC_1_size_m'] = HC1_pos_mods['HC_1_size_pc']*(1 * u .pc).to(u.m).value
                            HC1_pos_mods['HC_1_size_asec'] =  u_conversion.ang_size(D,HC1_pos_mods['HC_1_size_m'])
                            
                            HC1_pos_mods['HC_peak_v0_39_38_times']  = hc3n_prop['peak_v0_39_38'].loc[(hc3n_prop['Ind_ok'] == 1)].tolist()[0]/HC1_pos_mods['HC_1_peak_v0_39_38_fok']
                            HC1_pos_mods['HC_peak_v0_24_23_times']  = hc3n_prop['peak_v0_24_23'].loc[(hc3n_prop['Ind_ok'] == 1)].tolist()[0]/HC1_pos_mods['HC_1_peak_v0_24_23_fok']
                            HC1_pos_mods['HC_peak_v7_39_38_times']  = hc3n_prop['peak_v7_39_38'].loc[(hc3n_prop['Ind_ok'] == 1)].tolist()[0]/HC1_pos_mods['HC_1_peak_v7_39_38_fok']
                            HC1_pos_mods['HC_peak_v7_24_23_times']  = hc3n_prop['peak_v7_24_23'].loc[(hc3n_prop['Ind_ok'] == 1)].tolist()[0]/HC1_pos_mods['HC_1_peak_v7_24_23_fok']
                            
                            print HC1_pos_mods[['HC_1_size_asec', 'HC_1_size_pc', 'HC_1_Ltot']]
                            
                            HC2_nh_cond = [2.5e6]
                            HC2_T_cond = [300]
                            
                            HC2_size_fact = u_conversion.lin_size(3.5, 0.0219).to(u.pc)/2. # pc
                            
                            ncol_dust_cont['cond_HC2_T300'] = (ncol_dust_cont['central_nh'] == HC2_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC2_T_cond[0])
                            #ncol_dust_cont['cond_HC2_T300'] = (ncol_dust_cont['central_nh'] == HC2_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC2_T_cond[1])
                            #HC2_pos_mod1 = ncol_dust_cont[ncol_dust_cont['cond_HC2_T250'] == True] 
                            HC2_pos_mods= ncol_dust_cont[ncol_dust_cont['cond_HC2_T300'] == True] 
                            HC2_pos_mods['HC_2_Ltot']  = HC2_pos_mods['Ltot']*(HC2_size_fact**2)/(0.35**2)
                            HC2_pos_mods['HC_factor'] = HC2_size_fact
                            
                            HC2_pos_mods['HC_2_peak_v0_39_38_fok']  = HC2_pos_mods['peak_v0_39_38']*(HC2_size_fact**2)/(0.35**2)
                            HC2_pos_mods['HC_2_peak_v0_24_23_fok']  = HC2_pos_mods['peak_v0_24_23']*(HC2_size_fact**2)/(0.35**2)
                            HC2_pos_mods['HC_2_peak_v7_39_38_fok']  = HC2_pos_mods['peak_v7_39_38']*(HC2_size_fact**2)/(0.35**2)
                            HC2_pos_mods['HC_2_peak_v7_24_23_fok']  = HC2_pos_mods['peak_v7_24_23']*(HC2_size_fact**2)/(0.35**2)
                            HC2_pos_mods['HC_2_size_pc'] = 2. * HC2_size_fact
                            HC2_pos_mods['HC_2_size_m'] = HC2_pos_mods['HC_2_size_pc']*(1 * u .pc).to(u.m).value
                            HC2_pos_mods['HC_2_size_asec'] =  u_conversion.ang_size(D,HC2_pos_mods['HC_2_size_m'])
                            
                            HC2_pos_mods['HC_peak_v0_39_38_times']  = hc3n_prop['peak_v0_39_38'].loc[(hc3n_prop['Ind_ok'] == 2)].tolist()[0]/HC2_pos_mods['HC_2_peak_v0_39_38_fok']
                            HC2_pos_mods['HC_peak_v0_24_23_times']  = hc3n_prop['peak_v0_24_23'].loc[(hc3n_prop['Ind_ok'] == 2)].tolist()[0]/HC2_pos_mods['HC_2_peak_v0_24_23_fok']
                            HC2_pos_mods['HC_peak_v7_39_38_times']  = hc3n_prop['peak_v7_39_38'].loc[(hc3n_prop['Ind_ok'] == 2)].tolist()[0]/HC2_pos_mods['HC_2_peak_v7_39_38_fok']
                            HC2_pos_mods['HC_peak_v7_24_23_times']  = hc3n_prop['peak_v7_24_23'].loc[(hc3n_prop['Ind_ok'] == 2)].tolist()[0]/HC2_pos_mods['HC_2_peak_v7_24_23_fok']
                            
                            print HC2_pos_mods[['HC_2_size_asec', 'HC_2_size_pc', 'HC_2_Ltot']]
                            
                            HC3_nh_cond = [2.0e6, 2.5e6]
                            HC3_T_cond = [300, 350]
                            
                            HC3_size_fact = u_conversion.lin_size(3.5, 0.0219).to(u.pc)/2. # pc
                            
                            ncol_dust_cont['cond_HC3_T250'] = (ncol_dust_cont['central_nh'] == HC3_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC3_T_cond[0])
                            ncol_dust_cont['cond_HC3_T300'] = (ncol_dust_cont['central_nh'] == HC3_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC3_T_cond[1])
                            HC3_pos_mod1 = ncol_dust_cont[ncol_dust_cont['cond_HC3_T250'] == True] 
                            HC3_pos_mod2 = ncol_dust_cont[ncol_dust_cont['cond_HC3_T300'] == True] 
                            HC3_pos_mods = pd.concat([HC3_pos_mod1, HC3_pos_mod2], ignore_index=True)
                            HC3_pos_mods['HC_3_Ltot']  = HC3_pos_mods['Ltot']*(HC3_size_fact**2)/(0.35**2)
                            HC3_pos_mods['HC_factor'] = HC3_size_fact
                            
                            HC3_pos_mods['HC_3_peak_v0_39_38_fok']  = HC3_pos_mods['peak_v0_39_38']*(HC3_size_fact**2)/(0.35**2)
                            HC3_pos_mods['HC_3_peak_v0_24_23_fok']  = HC3_pos_mods['peak_v0_24_23']*(HC3_size_fact**2)/(0.35**2)
                            HC3_pos_mods['HC_3_peak_v7_39_38_fok']  = HC3_pos_mods['peak_v7_39_38']*(HC3_size_fact**2)/(0.35**2)
                            HC3_pos_mods['HC_3_peak_v7_24_23_fok']  = HC3_pos_mods['peak_v7_24_23']*(HC3_size_fact**2)/(0.35**2)
                            HC3_pos_mods['HC_3_size_pc'] = 2. * HC3_size_fact
                            HC3_pos_mods['HC_3_size_m'] = HC3_pos_mods['HC_3_size_pc']*(1 * u .pc).to(u.m).value
                            HC3_pos_mods['HC_3_size_asec'] =  u_conversion.ang_size(D,HC3_pos_mods['HC_3_size_m'])
                            
                            HC3_pos_mods['HC_peak_v0_39_38_times']  = hc3n_prop['peak_v0_39_38'].loc[(hc3n_prop['Ind_ok'] == 3)].tolist()[0]/HC3_pos_mods['HC_3_peak_v0_39_38_fok']
                            HC3_pos_mods['HC_peak_v0_24_23_times']  = hc3n_prop['peak_v0_24_23'].loc[(hc3n_prop['Ind_ok'] == 3)].tolist()[0]/HC3_pos_mods['HC_3_peak_v0_24_23_fok']
                            HC3_pos_mods['HC_peak_v7_39_38_times']  = hc3n_prop['peak_v7_39_38'].loc[(hc3n_prop['Ind_ok'] == 3)].tolist()[0]/HC3_pos_mods['HC_3_peak_v7_39_38_fok']
                            HC3_pos_mods['HC_peak_v7_24_23_times']  = hc3n_prop['peak_v7_24_23'].loc[(hc3n_prop['Ind_ok'] == 3)].tolist()[0]/HC3_pos_mods['HC_3_peak_v7_24_23_fok']
                            
                            print HC3_pos_mods[['HC_3_size_asec', 'HC_3_size_pc', 'HC_3_Ltot']]
                            
                            HC4_nh_cond = [5.0e6]
                            HC4_T_cond = [150, 200]
                            
                            HC4_size_fact = u_conversion.lin_size(3.5, 0.0192).to(u.pc)/2. # pc
                            
                            ncol_dust_cont['cond_HC4_T250'] = (ncol_dust_cont['central_nh'] == HC4_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC4_T_cond[0])
                            ncol_dust_cont['cond_HC4_T300'] = (ncol_dust_cont['central_nh'] == HC4_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC4_T_cond[1])
                            HC4_pos_mod1 = ncol_dust_cont[ncol_dust_cont['cond_HC4_T250'] == True] 
                            HC4_pos_mod2 = ncol_dust_cont[ncol_dust_cont['cond_HC4_T300'] == True] 
                            HC4_pos_mods = pd.concat([HC4_pos_mod1, HC4_pos_mod2], ignore_index=True)
                            HC4_pos_mods['HC_4_Ltot']  = HC4_pos_mods['Ltot']*(HC4_size_fact**2)/(0.35**2)
                            HC4_pos_mods['HC_factor'] = HC4_size_fact
                            
                            HC4_pos_mods['HC_4_peak_v0_39_38_fok']  = HC4_pos_mods['peak_v0_39_38']*(HC4_size_fact**2)/(0.35**2)
                            HC4_pos_mods['HC_4_peak_v0_24_23_fok']  = HC4_pos_mods['peak_v0_24_23']*(HC4_size_fact**2)/(0.35**2)
                            HC4_pos_mods['HC_4_peak_v7_39_38_fok']  = HC4_pos_mods['peak_v7_39_38']*(HC4_size_fact**2)/(0.35**2)
                            HC4_pos_mods['HC_4_peak_v7_24_23_fok']  = HC4_pos_mods['peak_v7_24_23']*(HC4_size_fact**2)/(0.35**2)
                            HC4_pos_mods['HC_4_size_pc'] = 2. * HC4_size_fact
                            HC4_pos_mods['HC_4_size_m'] = HC4_pos_mods['HC_4_size_pc']*(1 * u .pc).to(u.m).value
                            HC4_pos_mods['HC_4_size_asec'] =  u_conversion.ang_size(D,HC4_pos_mods['HC_4_size_m'])
                            
                            HC4_pos_mods['HC_peak_v0_39_38_times']  = hc3n_prop['peak_v0_39_38'].loc[(hc3n_prop['Ind_ok'] == 4)].tolist()[0]/HC4_pos_mods['HC_4_peak_v0_39_38_fok']
                            HC4_pos_mods['HC_peak_v0_24_23_times']  = hc3n_prop['peak_v0_24_23'].loc[(hc3n_prop['Ind_ok'] == 4)].tolist()[0]/HC4_pos_mods['HC_4_peak_v0_24_23_fok']
                            HC4_pos_mods['HC_peak_v7_39_38_times']  = hc3n_prop['peak_v7_39_38'].loc[(hc3n_prop['Ind_ok'] == 4)].tolist()[0]/HC4_pos_mods['HC_4_peak_v7_39_38_fok']
                            HC4_pos_mods['HC_peak_v7_24_23_times']  = hc3n_prop['peak_v7_24_23'].loc[(hc3n_prop['Ind_ok'] == 4)].tolist()[0]/HC4_pos_mods['HC_4_peak_v7_24_23_fok']
                            
                            print HC4_pos_mods[['HC_4_size_asec', 'HC_4_size_pc', 'HC_4_Ltot']]
                            
                            HC5_nh_cond = [1.0e8]
                            HC5_T_cond = [150, 200]
                            
                            HC5_size_fact = u_conversion.lin_size(3.5, 0.0216).to(u.pc)/2. # pc
                            
                            ncol_dust_cont['cond_HC5_T250'] = (ncol_dust_cont['central_nh'] == HC5_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC5_T_cond[0])
                            ncol_dust_cont['cond_HC5_T300'] = (ncol_dust_cont['central_nh'] == HC5_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC5_T_cond[1])
                            HC5_pos_mod1 = ncol_dust_cont[ncol_dust_cont['cond_HC5_T250'] == True] 
                            HC5_pos_mod2 = ncol_dust_cont[ncol_dust_cont['cond_HC5_T300'] == True] 
                            HC5_pos_mods = pd.concat([HC5_pos_mod1, HC5_pos_mod2], ignore_index=True)
                            HC5_pos_mods['HC_5_Ltot']  = HC5_pos_mods['Ltot']*(HC5_size_fact**2)/(0.35**2)
                            HC5_pos_mods['HC_factor'] = HC5_size_fact
                            
                            HC5_pos_mods['HC_5_peak_v0_39_38_fok']  = HC5_pos_mods['peak_v0_39_38']*(HC5_size_fact**2)/(0.35**2)
                            HC5_pos_mods['HC_5_peak_v0_24_23_fok']  = HC5_pos_mods['peak_v0_24_23']*(HC5_size_fact**2)/(0.35**2)
                            HC5_pos_mods['HC_5_peak_v7_39_38_fok']  = HC5_pos_mods['peak_v7_39_38']*(HC5_size_fact**2)/(0.35**2)
                            HC5_pos_mods['HC_5_peak_v7_24_23_fok']  = HC5_pos_mods['peak_v7_24_23']*(HC5_size_fact**2)/(0.35**2)
                            HC5_pos_mods['HC_5_size_pc'] = 2. * HC5_size_fact
                            HC5_pos_mods['HC_5_size_m'] = HC5_pos_mods['HC_5_size_pc']*(1 * u .pc).to(u.m).value
                            HC5_pos_mods['HC_5_size_asec'] =  u_conversion.ang_size(D,HC5_pos_mods['HC_5_size_m'])
                            
                            HC5_pos_mods['HC_peak_v0_39_38_times']  = hc3n_prop['peak_v0_39_38'].loc[(hc3n_prop['Ind_ok'] == 5)].tolist()[0]/HC5_pos_mods['HC_5_peak_v0_39_38_fok']
                            HC5_pos_mods['HC_peak_v0_24_23_times']  = hc3n_prop['peak_v0_24_23'].loc[(hc3n_prop['Ind_ok'] == 5)].tolist()[0]/HC5_pos_mods['HC_5_peak_v0_24_23_fok']
                            HC5_pos_mods['HC_peak_v7_39_38_times']  = hc3n_prop['peak_v7_39_38'].loc[(hc3n_prop['Ind_ok'] == 5)].tolist()[0]/HC5_pos_mods['HC_5_peak_v7_39_38_fok']
                            HC5_pos_mods['HC_peak_v7_24_23_times']  = hc3n_prop['peak_v7_24_23'].loc[(hc3n_prop['Ind_ok'] == 5)].tolist()[0]/HC5_pos_mods['HC_5_peak_v7_24_23_fok']
                            
                            print HC5_pos_mods[['HC_5_size_asec', 'HC_5_size_pc', 'HC_5_Ltot']]
                            
                            HC8_nh_cond = [2.0e6]
                            HC8_T_cond = [250]
                            
                            HC8_size_fact = u_conversion.lin_size(3.5, 0.0251).to(u.pc)/2. # pc
                            
                            ncol_dust_cont['cond_HC8_T250'] = (ncol_dust_cont['central_nh'] == HC8_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC8_T_cond[0])
                            HC8_pos_mods = ncol_dust_cont[ncol_dust_cont['cond_HC8_T250'] == True] 
                            HC8_pos_mods['HC_8_Ltot']  = HC8_pos_mods['Ltot']*(HC8_size_fact**2)/(0.35**2)
                            HC8_pos_mods['HC_factor'] = HC8_size_fact
                            
                            HC8_pos_mods['HC_8_peak_v0_39_38_fok']  = HC8_pos_mods['peak_v0_39_38']*(HC8_size_fact**2)/(0.35**2)
                            HC8_pos_mods['HC_8_peak_v0_24_23_fok']  = HC8_pos_mods['peak_v0_24_23']*(HC8_size_fact**2)/(0.35**2)
                            HC8_pos_mods['HC_8_peak_v7_39_38_fok']  = HC8_pos_mods['peak_v7_39_38']*(HC8_size_fact**2)/(0.35**2)
                            HC8_pos_mods['HC_8_peak_v7_24_23_fok']  = HC8_pos_mods['peak_v7_24_23']*(HC8_size_fact**2)/(0.35**2)
                            HC8_pos_mods['HC_8_size_pc'] = 2. * HC8_size_fact
                            HC8_pos_mods['HC_8_size_m'] = HC8_pos_mods['HC_8_size_pc']*(1 * u .pc).to(u.m).value
                            HC8_pos_mods['HC_8_size_asec'] =  u_conversion.ang_size(D,HC8_pos_mods['HC_8_size_m'])
                            
                            HC8_pos_mods['HC_peak_v0_39_38_times']  = hc3n_prop['peak_v0_39_38'].loc[(hc3n_prop['Ind_ok'] == 8)].tolist()[0]/HC8_pos_mods['HC_8_peak_v0_39_38_fok']
                            HC8_pos_mods['HC_peak_v0_24_23_times']  = hc3n_prop['peak_v0_24_23'].loc[(hc3n_prop['Ind_ok'] == 8)].tolist()[0]/HC8_pos_mods['HC_8_peak_v0_24_23_fok']
                            HC8_pos_mods['HC_peak_v7_39_38_times']  = hc3n_prop['peak_v7_39_38'].loc[(hc3n_prop['Ind_ok'] == 8)].tolist()[0]/HC8_pos_mods['HC_8_peak_v7_39_38_fok']
                            HC8_pos_mods['HC_peak_v7_24_23_times']  = hc3n_prop['peak_v7_24_23'].loc[(hc3n_prop['Ind_ok'] == 8)].tolist()[0]/HC8_pos_mods['HC_8_peak_v7_24_23_fok']
                            
                            
                            HC9_nh_cond = [2.5e6]
                            HC9_T_cond = [250]
                        
                            HC9_size_fact = u_conversion.lin_size(3.5, 0.0173).to(u.pc)/2. # pc
                            
                            ncol_dust_cont['cond_HC9_T250'] = (ncol_dust_cont['central_nh'] == HC9_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC9_T_cond[0])
                            HC9_pos_mods = ncol_dust_cont[ncol_dust_cont['cond_HC9_T250'] == True] 
                            HC9_pos_mods['HC_9_Ltot']  = HC9_pos_mods['Ltot']*(HC9_size_fact**2)/(0.35**2)
                            HC9_pos_mods['HC_factor'] = HC9_size_fact
                            
                            HC9_pos_mods['HC_9_peak_v0_39_38_fok']  = HC9_pos_mods['peak_v0_39_38']*(HC9_size_fact**2)/(0.35**2)
                            HC9_pos_mods['HC_9_peak_v0_24_23_fok']  = HC9_pos_mods['peak_v0_24_23']*(HC9_size_fact**2)/(0.35**2)
                            HC9_pos_mods['HC_9_peak_v7_39_38_fok']  = HC9_pos_mods['peak_v7_39_38']*(HC9_size_fact**2)/(0.35**2)
                            HC9_pos_mods['HC_9_peak_v7_24_23_fok']  = HC9_pos_mods['peak_v7_24_23']*(HC9_size_fact**2)/(0.35**2)
                            HC9_pos_mods['HC_9_size_pc'] = 2. * HC9_size_fact
                            HC9_pos_mods['HC_9_size_m'] = HC9_pos_mods['HC_9_size_pc']*(1 * u .pc).to(u.m).value
                            HC9_pos_mods['HC_9_size_asec'] =  u_conversion.ang_size(D,HC9_pos_mods['HC_9_size_m'])
                            
                            HC9_pos_mods['HC_peak_v0_39_38_times']  = hc3n_prop['peak_v0_39_38'].loc[(hc3n_prop['Ind_ok'] == 9)].tolist()[0]/HC9_pos_mods['HC_9_peak_v0_39_38_fok']
                            HC9_pos_mods['HC_peak_v0_24_23_times']  = hc3n_prop['peak_v0_24_23'].loc[(hc3n_prop['Ind_ok'] == 9)].tolist()[0]/HC9_pos_mods['HC_9_peak_v0_24_23_fok']
                            HC9_pos_mods['HC_peak_v7_39_38_times']  = hc3n_prop['peak_v7_39_38'].loc[(hc3n_prop['Ind_ok'] == 9)].tolist()[0]/HC9_pos_mods['HC_9_peak_v7_39_38_fok']
                            HC9_pos_mods['HC_peak_v7_24_23_times']  = hc3n_prop['peak_v7_24_23'].loc[(hc3n_prop['Ind_ok'] == 9)].tolist()[0]/HC9_pos_mods['HC_9_peak_v7_24_23_fok']
                            
                            
                            print HC9_pos_mods[['HC_9_size_asec', 'HC_9_size_pc', 'HC_9_Ltot']]
                            
                            HC11_nh_cond = [5.0e6]
                            HC11_T_cond = [200]
                        
                            HC11_size_fact = u_conversion.lin_size(3.5, 0.0206).to(u.pc)/2. # pc
                            
                            ncol_dust_cont['cond_HC11_T250'] = (ncol_dust_cont['central_nh'] == HC11_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC11_T_cond[0])
                            HC11_pos_mods = ncol_dust_cont[ncol_dust_cont['cond_HC11_T250'] == True] 
                            HC11_pos_mods['HC_11_Ltot']  = HC11_pos_mods['Ltot']*(HC11_size_fact**2)/(0.35**2)
                            HC11_pos_mods['HC_factor'] = HC11_size_fact
                            
                            HC11_pos_mods['HC_11_peak_v0_39_38_fok']  = HC11_pos_mods['peak_v0_39_38']*(HC11_size_fact**2)/(0.35**2)
                            HC11_pos_mods['HC_11_peak_v0_24_23_fok']  = HC11_pos_mods['peak_v0_24_23']*(HC11_size_fact**2)/(0.35**2)
                            HC11_pos_mods['HC_11_peak_v7_39_38_fok']  = HC11_pos_mods['peak_v7_39_38']*(HC11_size_fact**2)/(0.35**2)
                            HC11_pos_mods['HC_11_peak_v7_24_23_fok']  = HC11_pos_mods['peak_v7_24_23']*(HC11_size_fact**2)/(0.35**2)
                            HC11_pos_mods['HC_11_size_pc'] = 2. * HC11_size_fact
                            HC11_pos_mods['HC_11_size_m'] = HC11_pos_mods['HC_11_size_pc']*(1 * u .pc).to(u.m).value
                            HC11_pos_mods['HC_11_size_asec'] =  u_conversion.ang_size(D,HC11_pos_mods['HC_11_size_m'])
                            
                            HC11_pos_mods['HC_peak_v0_39_38_times']  = hc3n_prop['peak_v0_39_38'].loc[(hc3n_prop['Ind_ok'] == 11)].tolist()[0]/HC11_pos_mods['HC_11_peak_v0_39_38_fok']
                            HC11_pos_mods['HC_peak_v0_24_23_times']  = hc3n_prop['peak_v0_24_23'].loc[(hc3n_prop['Ind_ok'] == 11)].tolist()[0]/HC11_pos_mods['HC_11_peak_v0_24_23_fok']
                            HC11_pos_mods['HC_peak_v7_39_38_times']  = hc3n_prop['peak_v7_39_38'].loc[(hc3n_prop['Ind_ok'] == 11)].tolist()[0]/HC11_pos_mods['HC_11_peak_v7_39_38_fok']
                            HC11_pos_mods['HC_peak_v7_24_23_times']  = hc3n_prop['peak_v7_24_23'].loc[(hc3n_prop['Ind_ok'] == 11)].tolist()[0]/HC11_pos_mods['HC_11_peak_v7_24_23_fok']
                            
                            print HC11_pos_mods[['HC_11_size_asec', 'HC_11_size_pc', 'HC_11_Ltot']]
                            
                            HC13_nh_cond = [1.5e6]
                            HC13_T_cond = [350, 400]
                            
                            HC13_size_fact = u_conversion.lin_size(3.5, 0.0273).to(u.pc)/2. # pc
                            
                            ncol_dust_cont['cond_HC13_T250'] = (ncol_dust_cont['central_nh'] == HC13_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC13_T_cond[0])
                            ncol_dust_cont['cond_HC13_T300'] = (ncol_dust_cont['central_nh'] == HC13_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC13_T_cond[1])
                            HC13_pos_mod1 = ncol_dust_cont[ncol_dust_cont['cond_HC13_T250'] == True] 
                            HC13_pos_mod2 = ncol_dust_cont[ncol_dust_cont['cond_HC13_T300'] == True] 
                            HC13_pos_mods = pd.concat([HC13_pos_mod1, HC13_pos_mod2], ignore_index=True)
                            HC13_pos_mods['HC_13_Ltot']  = HC13_pos_mods['Ltot']*(HC13_size_fact**2)/(0.35**2)
                            HC13_pos_mods['HC_factor'] = HC13_size_fact
                            
                            HC13_pos_mods['HC_13_peak_v0_39_38_fok']  = HC13_pos_mods['peak_v0_39_38']*(HC13_size_fact**2)/(0.35**2)
                            HC13_pos_mods['HC_13_peak_v0_24_23_fok']  = HC13_pos_mods['peak_v0_24_23']*(HC13_size_fact**2)/(0.35**2)
                            HC13_pos_mods['HC_13_peak_v7_39_38_fok']  = HC13_pos_mods['peak_v7_39_38']*(HC13_size_fact**2)/(0.35**2)
                            HC13_pos_mods['HC_13_peak_v7_24_23_fok']  = HC13_pos_mods['peak_v7_24_23']*(HC13_size_fact**2)/(0.35**2)
                            HC13_pos_mods['HC_13_size_pc'] = 2. * HC13_size_fact
                            HC13_pos_mods['HC_13_size_m'] = HC13_pos_mods['HC_13_size_pc']*(1 * u .pc).to(u.m).value
                            HC13_pos_mods['HC_13_size_asec'] =  u_conversion.ang_size(D,HC13_pos_mods['HC_13_size_m'])
                            
                            HC13_pos_mods['HC_peak_v0_39_38_times']  = hc3n_prop['peak_v0_39_38'].loc[(hc3n_prop['Ind_ok'] == 13)].tolist()[0]/HC13_pos_mods['HC_13_peak_v0_39_38_fok']
                            HC13_pos_mods['HC_peak_v0_24_23_times']  = hc3n_prop['peak_v0_24_23'].loc[(hc3n_prop['Ind_ok'] == 13)].tolist()[0]/HC13_pos_mods['HC_13_peak_v0_24_23_fok']
                            HC13_pos_mods['HC_peak_v7_39_38_times']  = hc3n_prop['peak_v7_39_38'].loc[(hc3n_prop['Ind_ok'] == 13)].tolist()[0]/HC13_pos_mods['HC_13_peak_v7_39_38_fok']
                            HC13_pos_mods['HC_peak_v7_24_23_times']  = hc3n_prop['peak_v7_24_23'].loc[(hc3n_prop['Ind_ok'] == 13)].tolist()[0]/HC13_pos_mods['HC_13_peak_v7_24_23_fok']
                            
                            print HC13_pos_mods[['HC_13_size_asec', 'HC_13_size_pc', 'HC_13_Ltot']]
                            
                            HC14_nh_cond = [4.00e6]
                            HC14_T_cond = [250, 300]
                            
                            HC14_size_fact = u_conversion.lin_size(3.5, 0.0424).to(u.pc)/2. # pc
                            
                            ncol_dust_cont['cond_HC14_T250'] = (ncol_dust_cont['central_nh'] == HC14_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC14_T_cond[0])
                            ncol_dust_cont['cond_HC14_T300'] = (ncol_dust_cont['central_nh'] == HC14_nh_cond[0]) & (ncol_dust_cont['central_T'] == HC14_T_cond[1])
                            HC14_pos_mod1 = ncol_dust_cont[ncol_dust_cont['cond_HC14_T250'] == True] 
                            HC14_pos_mod2 = ncol_dust_cont[ncol_dust_cont['cond_HC14_T300'] == True] 
                            HC14_pos_mods = pd.concat([HC14_pos_mod1, HC14_pos_mod2], ignore_index=True)
                            HC14_pos_mods['HC_14_Ltot']  = HC14_pos_mods['Ltot']*(HC14_size_fact**2)/(0.35**2)
                            HC14_pos_mods['HC_factor'] = HC14_size_fact
                            
                            HC14_pos_mods['HC_14_peak_v0_39_38_fok']  = HC14_pos_mods['peak_v0_39_38']*(HC14_size_fact**2)/(0.35**2)
                            HC14_pos_mods['HC_14_peak_v0_24_23_fok']  = HC14_pos_mods['peak_v0_24_23']*(HC14_size_fact**2)/(0.35**2)
                            HC14_pos_mods['HC_14_peak_v7_39_38_fok']  = HC14_pos_mods['peak_v7_39_38']*(HC14_size_fact**2)/(0.35**2)
                            HC14_pos_mods['HC_14_peak_v7_24_23_fok']  = HC14_pos_mods['peak_v7_24_23']*(HC14_size_fact**2)/(0.35**2)
                            HC14_pos_mods['HC_14_size_pc'] = 2. * HC14_size_fact
                            HC14_pos_mods['HC_14_size_m'] = HC14_pos_mods['HC_14_size_pc']*(1 * u .pc).to(u.m).value
                            HC14_pos_mods['HC_14_size_asec'] =  u_conversion.ang_size(D,HC14_pos_mods['HC_14_size_m'])
                            
                            HC14_pos_mods['HC_peak_v0_39_38_times']  = hc3n_prop['peak_v0_39_38'].loc[(hc3n_prop['Ind_ok'] == 14)].tolist()[0]/HC14_pos_mods['HC_14_peak_v0_39_38_fok']
                            HC14_pos_mods['HC_peak_v0_24_23_times']  = hc3n_prop['peak_v0_24_23'].loc[(hc3n_prop['Ind_ok'] == 14)].tolist()[0]/HC14_pos_mods['HC_14_peak_v0_24_23_fok']
                            HC14_pos_mods['HC_peak_v7_39_38_times']  = hc3n_prop['peak_v7_39_38'].loc[(hc3n_prop['Ind_ok'] == 14)].tolist()[0]/HC14_pos_mods['HC_14_peak_v7_39_38_fok']
                            HC14_pos_mods['HC_peak_v7_24_23_times']  = hc3n_prop['peak_v7_24_23'].loc[(hc3n_prop['Ind_ok'] == 14)].tolist()[0]/HC14_pos_mods['HC_14_peak_v7_24_23_fok']
                            
                            print HC14_pos_mods[['HC_14_size_asec', 'HC_14_size_pc', 'HC_14_Ltot']]
                            
                        
                            HC_mods_lists = {'HC1': HC1_pos_mods,'HC2': HC2_pos_mods, 'HC3': HC3_pos_mods, 'HC4': HC4_pos_mods,
                                             'HC5': HC5_pos_mods, 'HC8' :HC8_pos_mods, 'HC9': HC9_pos_mods, 'HC11': HC11_pos_mods, 'HC13': HC13_pos_mods, 'HC14': HC14_pos_mods}
                            
                            
                            
                            
                            # Plotting line profiles:
                            plot_profiles = False
                            if plot_profiles == True:
                                model_list = ncol_dust_cont['model'].tolist()
                                #model_list = glob(direct+'*.inp')
                                #model= model_list[-1]
                                #if 1==1:
                                #for m, model in enumerate(model_list):
                                for m, model in ncol_dust_cont.iterrows():
                                    for tt, temp in enumerate(temp_list):
                                        model_path = res_path+'r_'+str(profile)+'/'+'r_'+str(profile)+'_T'+'%1.f' % temp
                                    
                                        
                                        model_name = model['model'].split('.inp')[0]
                                        
                                        
                                            
                                        
                                        print model_name + '\tT=' + '%1.f' % temp+ '\tN=' + '%1.1E' % model['Nline']
                                        #model_name = 'mhc3n_r0.0_Nd3.09E+21_X1.10E-08_nH3.00E+06'
                                        factor = 0.35
                                        plot_profiles = True
                                        plot_SED = True
                                        out_dir2 = res_path+'r_'+str(profile)+'/'+'Figures/SED_Profiles'
                                        if not os.path.exists(out_dir2):
                                            os.makedirs(out_dir2)
                                            
                                            
                                        
                                        #fig_format = '.png'
                                        luminosity_total, peak_v0_39_38, peak_v0_24_23, peak_v7_39_38, peak_v7_24_23, fig_PROF, fig_SED = luminosity_calculator(model_path, model_name, factor, plot_profiles, plot_SED, out_dir2, fig_format, hc3n_prop, temp)
                                        
                                        #fig_PROF.savefig(out_dir_prof+'/'+model_name+'_f'+'%1.1f' % factor +'_T'+'%1.f' % model['central_T']+'_profile.'+fig_format, bbox_inches='tight', transparent=True, dpi=400, format=fig_format)
                                        #plt.close(fig_PROF)
                                        
                                        #fig_SED.savefig(out_dir_sed+'/'+model_name+'_f'+'%1.1f' % factor +'_T'+'%1.f' % model['central_T']+'_SED.'+fig_format, bbox_inches='tight', transparent=True, dpi=400, format=fig_format)
                                        #plt.close(fig_SED)
                  
    
                          
    return HC_mods_lists, s_unique_dens, ncol_dust_cont





            
from math import atan2,degrees

#Label line with line2D label data
#def labelLine(line,x,label=None,align=True,**kwargs):
def labelLine(xvals,yvals,xpos, ypos,ax, color='k', label=None,align=True,**kwargs):
    xdata = xvals.tolist()
    ydata = yvals.tolist()
    x = xpos
    

    if (x < xdata[0]) or (x > xdata[-1]):
        print('x label location is outside data range!')
        return

    #Find corresponding y co-ordinate and angle of the line
    ip = 1
    for i in range(len(xdata)):
        if x < xdata[i]:
            ip = i
            break

    y = ydata[ip-1] + (ydata[ip]-ydata[ip-1])*(x-xdata[ip-1])/(xdata[ip]-xdata[ip-1])

    if not label:
        label = ''

    if align:
        #Compute the slope
        dx = xdata[ip] - xdata[ip-1]
        dy = ydata[ip] - ydata[ip-1]
        ang = degrees(atan2(dy,dx))

        #Transform to screen co-ordinates
        pt = np.array([x,y]).reshape((1,2))
        trans_angle = ax.transData.transform_angles(np.array((ang,)),pt)[0]

    else:
        trans_angle = 0

    #Set a bunch of keyword arguments
    if 'color' not in kwargs:
        kwargs['color'] = color

    if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
        kwargs['ha'] = 'center'

    if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
        kwargs['va'] = 'center'

    if 'backgroundcolor' not in kwargs:
        kwargs['backgroundcolor'] = ax.get_facecolor()

    if 'clip_on' not in kwargs:
        kwargs['clip_on'] = True

    if 'zorder' not in kwargs:
        kwargs['zorder'] = 2.5

    t = ax.text(x,ypos,label,rotation=trans_angle, fontsize=6, verticalalignment='center', horizontalalignment='center', **kwargs)
    t.set_bbox(dict(facecolor='white', alpha=0.0))
    
    
def luminosity_calculator_allinone(version, modeldf, fig_PROF, fig_SED, axis_PROF, axis_SED,
                                   profile='0.0', linewidth  = 0.7, labelsize = 12,
                                   v0_color = 'k', v7_color = 'red'):
    #model_path, model_name, factor, plot_profiles, plot_SED, figoutdir, fig_format, obs_df, Temp, fig_PROF, fig_SED):

    
    # Model path
    ttt= 'tcte_T'
    res_path = '/Users/frico/Documents/Ed/modelos_hc3n/v'+version+'_'+ttt+'/'+'r_'+profile+'/'
    
    model_temp = int(modeldf['Tdust'])
    model_path =res_path+'r_'+profile+'_T'+str(model_temp)
    model_name = modeldf['model']
    
    # Model resizing factor
    factor = modeldf['HC_factor']
    
    
    # =============================================================================
    # # Line profiles
    # =============================================================================
    dist_pc = 3.5e6 # pc
    luminosity_obs = 2.928e9 # lsun
    model_size = 1. * factor # pc 
    
    # Loading model line profile
    model_prof = pd.read_csv(model_path+'/'+model_name+'_.prof', delim_whitespace= True, skiprows=[0])
    model_prof.columns = ['lambda_um', 'I_ergscm2um', 'I_Ic_ergscm2um', 'V_kms']
    
    # Conversion erg/s/cm^2/um to mJy
    model_prof['I_Ic_mJy'] = u_conversion.ergscmum_to_mjy(model_prof['I_Ic_ergscm2um'], model_prof['lambda_um'])
    # Multiplying intensities by factor
    model_prof['I_Ic_mJy'] = model_prof['I_Ic_mJy']*(factor**2)
    
    
    #### J=24-23
    ### v=0
    alam_v0_24_23 = 1373.1495 # v0 24(+1)-23(-1)
    prof_v0_24_23_lim = []
    for i, row in model_prof.iterrows():
        if row['V_kms'] == 0 and abs(row['lambda_um']-alam_v0_24_23)<0.1:
            central_v0_24_23_index = i
            for j in range(i-1, 1, -1):
                if model_prof['V_kms'][j]>0:
                    prof_v0_24_23_lim.append(j+1)
                    break
            for k in range(i+1, model_prof.shape[0]):
                if model_prof['V_kms'][k] < 0:
                    prof_v0_24_23_lim.append(k-1)
                    break
    prof_v0_24_23 = model_prof.loc[prof_v0_24_23_lim[0]:prof_v0_24_23_lim[1]]
    peak_v0_24_23 = prof_v0_24_23['I_Ic_mJy'][prof_v0_24_23.V_kms == 0].max()
                
    ### v7=1
    alam_v7_24_23 = 1367.8303 # v7 24(+1)-23(-1)
    prof_v7_24_23_lim = []
    for i, row in model_prof.iterrows():
        if row['V_kms'] == 0 and abs(row['lambda_um']-alam_v7_24_23)<0.1:
            central_v7_24_23_index = i
            for j in range(i-0, 1, -1):
                if model_prof['V_kms'][j]>0:
                    prof_v7_24_23_lim.append(j+1)
                    break
            for k in range(i+1, model_prof.shape[0]):
                if model_prof['V_kms'][k] < 0:
                    prof_v7_24_23_lim.append(k-1)
                    break
    prof_v7_24_23 = model_prof.loc[prof_v7_24_23_lim[0]:prof_v7_24_23_lim[1]]
    peak_v7_24_23 = prof_v7_24_23['I_Ic_mJy'][prof_v7_24_23.V_kms == 0].max()
    
    #### J=39-38
    ### v0
    alam_v0_39_38 = 845.2062 # v0 39(+1)-38(-1)
    prof_v0_39_38_lim = []
    for i, row in model_prof.iterrows():
        if row['V_kms'] == 0 and abs(row['lambda_um']-alam_v0_39_38)<0.1:
            central_v0_39_38_index = i
            for j in range(i-0, 1, -1):
                if model_prof['V_kms'][j]>0:
                    prof_v0_39_38_lim.append(j+1)
                    break
            for k in range(i+1, model_prof.shape[0]):
                if model_prof['V_kms'][k] < 0:
                    prof_v0_39_38_lim.append(k-1)
                    break
    prof_v0_39_38 = model_prof.loc[prof_v0_39_38_lim[0]:prof_v0_39_38_lim[1]]
    peak_v0_39_38 = prof_v0_39_38['I_Ic_mJy'][prof_v0_39_38.V_kms == 0].max()
    
    ### v7=1
    alam_v7_39_38 = 843.1411 # v0 37(+1)-38(-1)
    prof_v7_39_38_lim = []
    for i, row in model_prof.iterrows():
        if row['V_kms'] == 0 and abs(row['lambda_um']-alam_v7_39_38)<0.1:
            central_v7_39_38_index = i
            for j in range(i-0, 1, -1):
                if model_prof['V_kms'][j]>0:
                    prof_v7_39_38_lim.append(j+1)
                    break
            for k in range(i+1, model_prof.shape[0]):
                if model_prof['V_kms'][k] < 0:
                    prof_v7_39_38_lim.append(k-1)
                    break
    prof_v7_39_38 = model_prof.loc[prof_v7_39_38_lim[0]:prof_v7_39_38_lim[1]]
    peak_v7_39_38 = prof_v7_39_38['I_Ic_mJy'][prof_v7_39_38.V_kms == 0].max()
    
    # Modeled values 
    hc3n_v0_ratio = peak_v0_39_38/peak_v0_24_23
    hc3n_v7_ratio = peak_v7_39_38/peak_v7_24_23
    hc3n_v0v7_3938_ratio = peak_v0_39_38/peak_v7_39_38
    hc3n_v0v7_2423_ratio = peak_v0_24_23/peak_v7_24_23
    
    if plot_profiles == True:
        for j in range(len(axis_PROF)):
            if j == 0:
                # Plotting spectrum 24-23
                axis_PROF[j].plot(prof_v0_24_23['V_kms'], prof_v0_24_23['I_Ic_mJy']*modeldf['HC_peak_v0_24_23_times'], linewidth=linewidth, color=v0_color)
                axis_PROF[j].plot(prof_v7_24_23['V_kms'], prof_v7_24_23['I_Ic_mJy']*modeldf['HC_peak_v7_24_23_times'], linewidth=linewidth, color=v7_color)
            elif j == 1:
                # Plotting spectrum 39-38
                axis_PROF[j].plot(prof_v0_39_38['V_kms'], prof_v0_39_38['I_Ic_mJy']*modeldf['HC_peak_v0_39_38_times'], linewidth=linewidth, color=v0_color)
                axis_PROF[j].plot(prof_v7_39_38['V_kms'], prof_v7_39_38['I_Ic_mJy']*modeldf['HC_peak_v7_39_38_times'], linewidth=linewidth, color=v7_color)

    # =============================================================================
    # ## SED and tau_continuum
    # =============================================================================
        
    
    ## ISO
    model_iso = pd.read_csv(model_path+'/'+model_name+'_.iso', delim_whitespace= True)
    model_iso.columns = ['lambda_um', 'I_Ic_ergscm2um', 'Ic_ergscm2um', 'Inorm_ergscm2um', 'I_ergscm2um']
    # Conversion erg/s/cm^2/um to mJy
    model_iso['Ic_Jy'] = u_conversion.ergscmum_to_mjy(model_iso['Ic_ergscm2um'], model_iso['lambda_um'])/1000
    # Multiplying intensities by factor
    model_iso['Ic_Jy'] = model_iso['Ic_Jy']*(factor**2)
    model_iso['Ic_ergscm2um'] = model_iso['Ic_ergscm2um']*(factor**2)
    
    ## Spire
    model_spire = pd.read_csv(model_path+'/'+model_name+'_.spire', delim_whitespace= True)
    model_spire.columns = ['lambda_um', 'I_Ic_ergscm2um', 'Ic_ergscm2um', 'Inorm_ergscm2um', 'I_ergscm2um']
    # Conversion erg/s/cm^2/um to mJy
    model_spire['Ic_Jy'] = u_conversion.ergscmum_to_mjy(model_spire['Ic_ergscm2um'], model_spire['lambda_um'])/1000
    # Multiplying intensities by factor
    model_spire['Ic_Jy'] = model_spire['Ic_Jy']*(factor**2)
    model_spire['Ic_ergscm2um'] = model_spire['Ic_ergscm2um']*(factor**2)
    
    # Optical Depths
    model_taudust = pd.read_csv(model_path+'/'+model_name+'_.taudust', delim_whitespace= True)
    model_taudust.columns = ['lambda_um', 'taudust']
    
    # =============================================================================
    # Luminosity 10-1200um
    # =============================================================================
        
    luminosity_iso = 0
    for i, row in model_iso.iterrows():
        if i == 0:
            continue
        else:
            luminosity_iso += (model_iso['Ic_ergscm2um'][i]*
                              (model_iso['lambda_um'][i]-model_iso['lambda_um'][i-1])*
                              4.*3.14159*(3.5e6*3.086e18)**2)/(3.8e33)
    luminosity_spire = 0
    for i, row in model_spire.iterrows():
        if i == 0:
            continue
        else:
            luminosity_spire += (model_spire['Ic_ergscm2um'][i]*
                              (model_spire['lambda_um'][i]-model_spire['lambda_um'][i-1])*
                              4.*3.14159*(3.5e6*3.086e18)**2)/(3.8e33)
        
    luminosity_total = luminosity_iso+luminosity_spire # lsun Para un tamaño de 0.1" si factor = 0.7
    luminosity_factor = luminosity_total/luminosity_obs
    
    
    orig_luminosity_factor = (luminosity_total/(factor**2))/luminosity_obs
    new_factor = np.sqrt(1./ orig_luminosity_factor)
    new_luminosity_total =(luminosity_total/(factor**2))*new_factor # lsun
    new_model_size = 1. * new_factor # pc

    plot_SED = True
    if plot_SED == True:
        
        for j in range(len(axis_SED)):
            # SED
            if j == 0:
                # SED
                axis_SED[j].plot(model_iso['lambda_um'], model_iso['Ic_Jy'], linewidth=linewidth, color='0.7')
                axis_SED[j].plot(model_spire['lambda_um'], model_spire['Ic_Jy'], linewidth=linewidth, color='0.7')

                
            # Optical depth
            elif j==1:
                axis_SED[j].plot(model_taudust['lambda_um'], model_taudust['taudust'], linewidth=linewidth, color='0.7')


#        axis_SED[0].text(0.1, 0.9,
#                            'L='+ '% 1.2E' % luminosity_total + r' L$_\odot$', ha='left', va='center',
#                            transform=axis_SED[j].transAxes,
#                            rotation='horizontal', fontsize=8)
#        axis_SED[0].text(0.5, 0.9,
#                            'ss='+ '% 1.3f' % model_size + r' arcsec', ha='left', va='center',
#                            transform=axis_SED[j].transAxes,
#                            rotation='horizontal', fontsize=8)
        
    
    
    return fig_PROF, fig_SED, model_iso, model_spire


def plot_PROF_ini(obs_df, velocity_limits = ['NA', 'NA'], intensity_limits = ['NA', 'NA']):
    line_names = [r'HC$_3$N 24(+1)-23(-1)', r'HC$_3$N 39(+1)-38(-1)']
    linewidth  = 0.7
    intensity_limits=['NA', 'NA']
    velocity_limits=[-100, 100]
    label_y =  r"Flux density (mJy)" 
    label_x = r"V (km/s)"
    labelsize = 12
    filas = 2
    columnas = 1
    first_last_xlabel = False
    v0_color = 'red'
    v7_color = 'k'
    m = filas
    n = columnas
    dx, dy = 1.2, 1
    figsize = plt.figaspect(float(dy * m) / float(dx * n))
    fig_PROF = plt.figure(figsize=figsize)
    gs1 = gridspec.GridSpec(m, n)    
    gs1.update(wspace = 0.0, hspace=0.0, top=0.95, bottom = 0.05)
    # Generating specified number of axis
    axis = []
    axis_ind = []
    for i in range(m*n):
        axis.append(fig_PROF.add_subplot(gs1[i]))
    # Generating axis index
    ind = 0
    axis_ind = []
    for i  in range(m):
        axis_ind.append([])
        for j in range(n):
            axis_ind[i].append(ind)
            ind += 1
    for j in range(m*n):
        # Plotting J=24(+1)-23(-1)
        if j == 0:
            # Line name
            axis[j].text(0.1, 0.9,
                        line_names[j], ha='left', va='center',
                        transform=axis[j].transAxes,
                        rotation='horizontal', fontsize=8)
            
#            # Obs
#            obs_vs = np.linspace(-50, 50, 14)
#            for o, orow in obs_df.iterrows():
#                #if orow['peak_v0_39_38'] != np.nan:
#                if not pd.isnull(orow['peak_v0_24_23']):
#                    axis[j].text(obs_vs[o], orow['peak_v0_24_23'],
#                            str(orow['Ind_ok']), ha='left', va='center', color=v0_color,
#                            #transform=axis[j].transAxes,
#                            rotation='horizontal', fontsize=8)
#                    
#                if not pd.isnull(orow['peak_v7_24_23']):
#                    axis[j].text(obs_vs[o], orow['peak_v7_24_23'],
#                            str(orow['Ind_ok']), ha='left', va='center', color=v7_color,
#                            #transform=axis[j].transAxes,
#                            rotation='horizontal', fontsize=8)
                    
        elif j == 1:
            # Plotting spectrum
            # Line name
            axis[j].text(0.1, 0.9,
                        line_names[j], ha='left', va='center',
                        transform=axis[j].transAxes,
                        rotation='horizontal', fontsize=8)
            
            # Obs
#            obs_vs = np.linspace(-50, 50, 14)
#            for o, orow in obs_df.iterrows():
#                #if orow['peak_v0_39_38'] != np.nan:
#                if not pd.isnull(orow['peak_v0_39_38']):
#                    axis[j].text(obs_vs[o], orow['peak_v0_39_38'],
#                            str(orow['Ind_ok']), ha='left', va='center', color=v0_color,
#                            #transform=axis[j].transAxes,
#                            rotation='horizontal', fontsize=8)
#                    
#                if not pd.isnull(orow['peak_v7_39_38']):
#                    axis[j].text(obs_vs[o], orow['peak_v7_39_38'],
#                            str(orow['Ind_ok']), ha='left', va='center', color=v7_color,
#                            #transform=axis[j].transAxes,
#                            rotation='horizontal', fontsize=8)
    for i, ax in enumerate(axis):
        #ax.tick_params(direction='in')
        
        majorLocator = MultipleLocator(20)
        ax.yaxis.set_major_locator(majorLocator)
        # Solo ponemos titulo a los ejes y si estan a la izquierda del todo
        ax.minorticks_on()
        #ax.tick_params(axis='both', which='minor', direction='in', top=True, right=True)
        ax.tick_params(axis='both', which='major', labelsize=labelsize, length=5, direction='in')
        ax.tick_params(axis='both', which='minor', labelsize=labelsize, direction='in')
        ax.xaxis.set_tick_params(top ='on', which='both')
        ax.yaxis.set_tick_params(right='on', which='both', labelright='off')
        
        # Left Figures
        left_ind = []
        for i in range(m):
            left_ind.append(axis_ind[i][0])
        # Bottom Figures
        bottom_figures = axis_ind[-1]
        
        # Limite en el eje x
        if velocity_limits != ['NA', 'NA']: # Si no se especifican limites python coge el minimo y maximo de cada expectro
            ax.set_xlim(velocity_limits)    # Limites introducidos
        # Limite en el eje y
        if intensity_limits != ['NA', 'NA']: # Si no especificamos limites ejey coge el maximo y minimo disponible del espectro
            ax.set_ylim(intensity_limits)
        if i in left_ind:
            ax.set_ylabel(label_y, fontsize=14)
        else:
            ax.set_yticklabels([])
        # Solo ponemos titulo a los ejes x si estan abajo del todo
        if i in bottom_figures:
            ax.set_xlabel(label_x, fontsize=14)
            if first_last_xlabel == False:
                plt.setp(ax.get_xticklabels()[0], visible=False)    
                plt.setp(ax.get_xticklabels()[-1], visible=False)
        else:
            ax.set_xticklabels([])
    return fig_PROF, axis

def plot_SED_ini(linewidth=0.7, labelsize=12):
    m = 2
    n = 1
    # Figure Size
    size_rat = float(n)/float(m)
    size_y = 10.*size_rat
    size_x = 10.
    fig_SED = plt.figure(figsize=(size_y, size_x))
    gs1 = gridspec.GridSpec(m, n)    
    gs1.update(wspace = 0.0, hspace=0.0, top=0.95, bottom = 0.05)
    # Generating specified number of axis
    axis = []
    for i in range(m*n):
        axis.append(fig_SED.add_subplot(gs1[i]))
        
        
    # observational SED from perez-beaupuits 2018
    obs_SED = pd.read_csv('/Users/frico/Documents/data/NGC253_H3O+/NGC253_SED.txt', delim_whitespace = True, comment='#')
    SED_model_pars =pd.read_csv('/Users/frico/Documents/data/NGC253_H3O+/NGC253_SED_fit_parameters.txt', delim_whitespace = True, comment='#')
    modeled_comp_cold =[]
    modeled_comp_warm =[]
    modeled_comp_hot =[]
    modeled_comp_all = []
    wave_mod = range(10, 1200+1)
    D = 3.52
    for i, wave in enumerate(wave_mod):
        nu = (wave*u.um).to(u.GHz, equivalencies=u.spectral())
        df_cold = SED_model_pars.loc[SED_model_pars['Component'] == 'Cold']
        df_warm = SED_model_pars.loc[SED_model_pars['Component'] == 'Warm']
        df_hot = SED_model_pars.loc[SED_model_pars['Component'] == 'Hot']
        
        
        modeled_cold = utiles.SED_model(nu, df_cold['Tdust'][0], df_cold['Mdust_1e6'][0]*1e6, df_cold['filling_fact'][0], D)
        modeled_warm = utiles.SED_model(nu, df_warm['Tdust'][1], df_warm['Mdust_1e6'][1]*1e6, df_warm['filling_fact'][1], D)
        modeled_hot = utiles.SED_model(nu, df_hot['Tdust'][2], df_hot['Mdust_1e6'][2]*1e6, df_hot['filling_fact'][2], D)
        
        modeled_comp_cold.append(modeled_cold.value)
        modeled_comp_warm.append(modeled_warm.value)
        modeled_comp_hot.append(modeled_hot.value)
        modeled_comp_all.append(modeled_cold.value+modeled_warm.value+modeled_hot.value)
        
    for j in range(m*n):
        # SED
        if j == 0:
            # SED
            #axis[j].plot(model_iso['lambda_um'], model_iso['Ic_Jy'], linewidth=linewidth, color='k')
            #axis[j].plot(model_spire['lambda_um'], model_spire['Ic_Jy'], linewidth=linewidth, color='k')
            # Observed SED
            axis[j].errorbar(obs_SED['Wavlength'], obs_SED['Obs_flux'],
                            yerr=obs_SED['Obs_flux_err'], fmt='.',
                            elinewidth=0.7, capsize=0.5,
                            linewidth=linewidth, color='k')
            # Modeled SED
            axis[j].plot(wave_mod, modeled_comp_cold, linewidth=linewidth, color='blue')
            axis[j].plot(wave_mod, modeled_comp_warm, linewidth=linewidth, color='green')
            axis[j].plot(wave_mod, modeled_comp_hot, linewidth=linewidth, color='red')
            axis[j].plot(wave_mod, modeled_comp_all, linewidth=linewidth, color='k', linestyle='--')
        # Optical depth
        elif j==1:
            # Fit opt depth from SED
            axis[j].errorbar(obs_SED['Wavlength'], obs_SED['opt_depth'],
                    yerr=obs_SED['Obs_flux_err'], fmt='.',
                    elinewidth=0.7, capsize=0.5,
                    linewidth=linewidth, color='k')
            
    axis[0].set_ylim([1e-2,1e4])
    
    
    
    
    axis[0].set_ylabel(r'Flux density (Jy)', fontsize=14)
    axis[1].set_ylabel(r'Optical depth $\tau$', fontsize=14)
    axis[1].set_xlabel(r'Wavelength ($\mu$m)', fontsize=14)
    
    for i, ax in enumerate(axis):
        majorLocator = MultipleLocator(20)
        ax.yaxis.set_major_locator(majorLocator)
        # Solo ponemos titulo a los ejes y si estan a la izquierda del todo
        ax.minorticks_on()
        #ax.tick_params(axis='both', which='minor', direction='in', top=True, right=True)
        ax.tick_params(axis='both', which='major', labelsize=labelsize, length=5, direction='in')
        ax.tick_params(axis='both', which='minor', labelsize=labelsize, direction='in')
        ax.xaxis.set_tick_params(top ='on', which='both')
        ax.yaxis.set_tick_params(right='on', which='both', labelright='off')
        ax.set_xlim([1,1200])
        ax.set_yscale('log')
        ax.set_xscale('log')

    axis[0].set_xticklabels([])

    return fig_SED, axis
    
def labelLines(xvals, yvals, ax, xpos, ypos, color, label, align=True,**kwargs):
    labelLine(xvals,yvals,xpos,ypos,ax,color,label,align,**kwargs)

"""
Results version
"""    
version = '11_a'
"""
Densities and profile
"""   
profile_list = [0.0]#, 1.0, 1.5, 2.0]
start_dens_list = [5e5, 1e6, 5e7, 1e7, 5e7, 1e8, 5e8]
"""
Abundances
"""
high_xlist = [5E-6, 5E-7, 1E-7, 5E-8]
low_xlist = [1E-7,1E-8, 1E-9, 1E-10]
Xline_list = low_xlist
"""
Temperatures, for Tgrad = False
"""
T_list = [100, 150, 200, 250, 300, 350]#[100, 200, 250, 300, 350, 400, 500, 600]
"""
Luminosities, for Tgrad = True
"""
luminosity_list  = [5e7, 1e8,5e8, 1e9, 2.9e9] 
"""
Observed Data
"""
hc3n_path ='Hotcores_v4_all'+'/hc3n_obs_results_v0.txt'   


"""
    xaxis
        central_nh
        central_T
    yaxis
        ratios -> line ratios
        taus -> line opacities
        intens -> line intensities
"""
    
yaxis_list = ['ratios']#['ratios', 'taus']#, 'intens']
cont_plot = True
xaxis = 'central_T'
for y, yaxis in enumerate(yaxis_list):
    HC_mods_lists, s_unique_dens, mod_df = plot_nh_ch2_fun(t_grad= False, xaxis = xaxis, yaxis=yaxis, 
                        profile_list=profile_list, start_dens_list=start_dens_list, 
                        T_list = T_list, 
                        Xline_list_or_Nline_list=Xline_list, luminosity_list=luminosity_list,
                        version=version, hc3n_path=hc3n_path, cont_plot=cont_plot,
                        plot_by_X=True, plot_sep=False, plot_taus=False, plot_intens = False)
    
Ltotal_lowt = HC_mods_lists['HC1']['HC_1_Ltot'].tolist()[0] + HC_mods_lists['HC2']['HC_2_Ltot'].tolist()[0] + HC_mods_lists['HC3']['HC_3_Ltot'].tolist()[0] + HC_mods_lists['HC4']['HC_4_Ltot'].tolist()[0] + HC_mods_lists['HC5']['HC_5_Ltot'].tolist()[0] + HC_mods_lists['HC8']['HC_8_Ltot'].tolist()[0]+ HC_mods_lists['HC9']['HC_9_Ltot'].tolist()[0]+ HC_mods_lists['HC11']['HC_11_Ltot'].tolist()[0] + HC_mods_lists['HC13']['HC_13_Ltot'].tolist()[0] +  HC_mods_lists['HC14']['HC_14_Ltot'].tolist()[0] 
Ltotal_hight = HC_mods_lists['HC1']['HC_1_Ltot'].tolist()[1] + HC_mods_lists['HC2']['HC_2_Ltot'].tolist()[0] + HC_mods_lists['HC3']['HC_3_Ltot'].tolist()[1] + HC_mods_lists['HC4']['HC_4_Ltot'].tolist()[1] + HC_mods_lists['HC5']['HC_5_Ltot'].tolist()[1] + HC_mods_lists['HC8']['HC_8_Ltot'].tolist()[0]+ HC_mods_lists['HC9']['HC_9_Ltot'].tolist()[0]+ HC_mods_lists['HC11']['HC_11_Ltot'].tolist()[0] + + HC_mods_lists['HC13']['HC_13_Ltot'].tolist()[1] +  HC_mods_lists['HC14']['HC_14_Ltot'].tolist()[1] 



# Models ok
HC_mods_ok_highT = {}
HC_mods_ok_lowT = {}
for i, num in enumerate([1,2,3,4,5,8,9,11,13,14]):
    print num
    #print HC_mods_lists['HC'+str(num)][['HC_'+str(num)+'_peak_v7_24_23_times']]
    #print HC_mods_lists['HC'+str(num)][['HC_'+str(num)+'_size_asec', 'HC_'+str(num)+'_Ltot']]
    #print HC_mods_lists['HC'+str(num)]['HC_'+str(num)+'_size_pc']
    # If wechange the source size of the model we change everything after. (we go from 1pc to xpc size)
    # New NH column density
#    HC_mods_lists['HC'+str(num)]['HC'+str(num)+'_new_NH'] = (HC_mods_lists['HC'+str(num)]['HC_'+str(num)+'_size_pc']*3.085677E18/2.)*HC_mods_lists['HC'+str(num)]['nH']
#    # New HC3N column density
#    HC_mods_lists['HC'+str(num)]['HC'+str(num)+'_new_NHC3N'] = HC_mods_lists['HC'+str(num)]['Xline']*HC_mods_lists['HC'+str(num)]['HC'+str(num)+'_new_NH'] 
#    HC_mods_lists['HC'+str(num)]['HC'+str(num)+'_new_Ndust'] = HC_mods_lists['HC'+str(num)]['Xdust']*HC_mods_lists['HC'+str(num)]['HC'+str(num)+'_new_NH'] 
#
#    # We mantain desired columns but since we change the size, this changes the abundances
#    HC_mods_lists['HC'+str(num)]['HC'+str(num)+'_new_Xd'] = HC_mods_lists['HC'+str(num)]['Ndust']/((HC_mods_lists['HC'+str(num)]['HC_'+str(num)+'_size_pc']*3.085677E18/2.)*HC_mods_lists['HC'+str(num)]['nH'])
#    HC_mods_lists['HC'+str(num)]['HC'+str(num)+'_new_Xline'] = HC_mods_lists['HC'+str(num)]['Nline']/((HC_mods_lists['HC'+str(num)]['HC_'+str(num)+'_size_pc']*3.085677E18/2.)*HC_mods_lists['HC'+str(num)]['nH'])
#    HC_mods_lists['HC'+str(num)]['HC'+str(num)+'_Mass'] = HC_mods_lists['HC'+str(num)]['central_nh']*(4./3)*np.pi*((HC_mods_lists['HC'+str(num)]['HC_'+str(num)+'_size_pc']*3.085677E18/2.)**3)
#    HC_mods_lists['HC'+str(num)]['HC'+str(num)+'_Mass_msun'] = HC_mods_lists['HC'+str(num)]['HC'+str(num)+'_Mass'] *2.*1.6737236E-24/(1 * u.Msun).to(u.g).value
    # New NH column density
    HC_mods_lists['HC'+str(num)]['new_NH'] = (HC_mods_lists['HC'+str(num)]['HC_'+str(num)+'_size_pc']*3.085677E18/2.)*HC_mods_lists['HC'+str(num)]['nH']
    # New HC3N column density
    HC_mods_lists['HC'+str(num)]['new_NHC3N'] = HC_mods_lists['HC'+str(num)]['Xline']*HC_mods_lists['HC'+str(num)]['new_NH'] 
    HC_mods_lists['HC'+str(num)]['new_Ndust'] = HC_mods_lists['HC'+str(num)]['Xdust']*HC_mods_lists['HC'+str(num)]['new_NH'] 

    # We mantain desired columns but since we change the size, this changes the abundances
    HC_mods_lists['HC'+str(num)]['new_Xd'] = HC_mods_lists['HC'+str(num)]['Ndust']/((HC_mods_lists['HC'+str(num)]['HC_'+str(num)+'_size_pc']*3.085677E18/2.)*HC_mods_lists['HC'+str(num)]['nH'])
    HC_mods_lists['HC'+str(num)]['new_Xline'] = HC_mods_lists['HC'+str(num)]['Nline']/((HC_mods_lists['HC'+str(num)]['HC_'+str(num)+'_size_pc']*3.085677E18/2.)*HC_mods_lists['HC'+str(num)]['nH'])
    HC_mods_lists['HC'+str(num)]['Mass'] = HC_mods_lists['HC'+str(num)]['central_nh']*(4./3)*np.pi*((HC_mods_lists['HC'+str(num)]['HC_'+str(num)+'_size_pc']*3.085677E18/2.)**3)
    HC_mods_lists['HC'+str(num)]['Mass_msun'] = HC_mods_lists['HC'+str(num)]['Mass'] *2.*1.6737236E-24/(1 * u.Msun).to(u.g).value
    
    
    HC_mods_lists['HC'+str(num)]['Size_pc'] = HC_mods_lists['HC'+str(num)]['HC_'+str(num)+'_size_pc']
    HC_mods_lists['HC'+str(num)]['Size_asec'] = HC_mods_lists['HC'+str(num)]['HC_'+str(num)+'_size_asec']
    HC_mods_lists['HC'+str(num)]['SHC_Ltot'] = HC_mods_lists['HC'+str(num)]['HC_'+str(num)+'_Ltot']
    print HC_mods_lists['HC'+str(num)][['tau_v0_24_23', 'tau_v0_39_38', 'tau_v7_24_23', 'tau_v7_39_38']]
    #print HC_mods_lists['HC'+str(num)][['HC'+str(num)+'_new_NH', 'HC'+str(num)+'_new_Xline', 'HC'+str(num)+'_new_Xd']]
    
    #print HC_mods_lists['HC'+str(num)][['HC'+str(num)+'_new_NH', 'HC'+str(num)+'_new_Ndust', 'HC'+str(num)+'_new_NHC3N', 'Xline', 'Xdust']]
    #print 'Xline'
    #print HC_mods_lists['HC'+str(num)]['Xline']
    #print 'Xline_new'
    #print HC_mods_lists['HC'+str(num)]['HC'+str(num)+'_new_NHC3N']/HC_mods_lists['HC'+str(num)]['HC'+str(num)+'_new_NH']
    
    #print HC_mods_lists['HC'+str(num)]['HC'+str(num)+'_Mass_msun']/1000.
    HC_mods_ok_highT['HC'+str(num)] = HC_mods_lists['HC'+str(num)].loc[HC_mods_lists['HC'+str(num)]['Tdust'].idxmax()]
    HC_mods_ok_lowT['HC'+str(num)] = HC_mods_lists['HC'+str(num)].loc[HC_mods_lists['HC'+str(num)]['Tdust'].idxmin()]
    

df_list = []
for i, num in enumerate([1,2,3,4,5,8,9,11,13,14]):
    
    HC_mods_ok_highT['HC'+str(num)]['Ind'] = num
    df_list.append(HC_mods_ok_highT['HC'+str(num)])
    
# Merging models to a same Df
df_highT2 = pd.DataFrame(df_list)
df_highT2.set_index('Ind', inplace=True)
# Saving to latex
table_out = ['Tdust', 'nH', 'Size_pc', 'Size_asec',
             'new_NH', 'new_NHC3N','Xline', 'new_Ndust',
             'Mass_msun', 'SHC_Ltot']

df_highT2.to_csv('/Users/frico/Documents/data/NGC253_H3O+/'+'Hotcores_v4_all/model_results_summary_v1.txt', columns=table_out, header=True, sep='\t', index=True, na_rep=-1)


sizes_pc_yo = hc3n_prop['Source_size_pc_v7_219_g'].tolist()
sizes_pc_yo[5] = 1.7
sizes_pc_yo[6] = 1.7

# Leroy Stuff
leroy_df= deepcopy(df_highT2[table_out])
# Filling for 6 and 7, 10, 12 rows with empty non-LTE modeling
for n in [6,7,10,12]:
    leroy_df.loc[n] = [np.nan]*len(leroy_df.columns)
leroy_df.sort_index(inplace=True)
L_IR_TOTAL_GA15 = 3.2E10
hc3n_prop.drop([14], inplace=True)
L_protoSSC_LTE_old = hc3n_prop['L_Lsun_v7_219_g'].tolist()
#Tvib_LTE = hc3n_prop['Tvib'].tolist()
# Nuevo limites! (para los que no tengo Tvib, pongo Trot v=0)
Tvib_LTE = [216., 304., 337., 326.,269.,92.,95.,217.,88.,84.,113.,76.,393.,312.]
# L_protoSSC_LTE
L_protoSSC_LTE = []
s_si = _si.sigma_sb
for i, temp in enumerate(Tvib_LTE):
    Lum = (((4*np.pi*(((sizes_pc_yo[i]*u.pc/2.).to(u.cm))**2)*s_si*((u.K * temp)**4)).to(u.Lsun)).value)
    L_protoSSC_LTE.append(Lum)
    
    
    
### Para mantener los de la tabla:
L_protoSSC_LTE = [1.132087E8, 5.242486E8, 7.866067E8, 5.289255E8, 3.112957E8, 0.0, 0.0, 1.782151E8, 0.022827E8, 0.012235E8, 0.087843E8, 0.008873E8, 22.7424252E8, 21.68234234E8]     
Tex_nLTe = leroy_df['Tdust'].tolist()
L_nLTE = leroy_df['SHC_Ltot'].tolist()
# Give zero lum to those with no det
L_nLTE[5] = 0
L_nLTE[6] = 0
L_nLTE[9] = 0
L_nLTE[11] = 0
L_protoSSC_LTE[5] = 0
L_protoSSC_LTE[6] = 0



NHC3N_yo = (10**hc3n_prop['logN_vib']).tolist()
NHC3N_yo_nLTE = leroy_df['new_NHC3N'].tolist()
leroy_sizes_pc = [2.7, 1.2, 2.6, 2.5, 2.1, 2.1, 2.9, 1.9, 2.6, 3.5, 2.9, 4.3, 1.6, 1.6]
leroy_gas_mass = [10**4.9, 10**4.7, 10**5.1, 10**5.1, 10**5.3, 10**3.6, 10**4.5, 10**5.2, 10**4.7, 10**5.2, 10**4.5, 10**4.1, 10**5.2, 10**5.7]
leroy_star_mass = [10**4.3, 10**4.3, 10**4.1, 10**5.0, 10**5.4, 10**5.3, 10**4.5, 10**4.8, 10**5.5, 10**5.3, 10**5.6, 10**6.0, 10**4.8, 10**5.5]
l_ir_nuc = 1.8E10
lum_to_mass = 1000.
hydrogen_mass_msun = 2.*(_si.m_p + _si.m_e).to(u.Msun).value


nh2_leroy_list = []
Nh2_leroy_list = []
Xhc3n_list = []
Xhc3n_list_nLTE = []
L_ir_leroy = 0
L_ir_leroy_list = []
T_cond = []
Tcond_2 = []
Tcond_3 = []
L_msSSC_L_protoSSC_ratio = []
sizes_Ltot = []
sizes_Ltot_nLTE = []
Tcond_nLTE = []
Tcond_nLTE_2 = []
ratio_LTE= []
ratio_nLTE = []
for i, mass in enumerate(leroy_star_mass):
    L_ir_leroy = L_ir_leroy + (lum_to_mass*mass)
    L_ir_leroy_list.append(lum_to_mass*mass)
    
    vol = (4./3)*np.pi*((leroy_sizes_pc[i]/2.)*(1 * u.pc).to(u.cm).value)**3
    
    nh2_leroy_list.append(leroy_gas_mass[i]/(vol*hydrogen_mass_msun))
    Nh2_leroy_list.append((leroy_sizes_pc[i]/2.)*(1 * u.pc).to(u.cm).value*nh2_leroy_list[i])
    Xhc3n_list.append(NHC3N_yo[i]/Nh2_leroy_list[i])
    Xhc3n_list_nLTE.append(NHC3N_yo_nLTE[i]/Nh2_leroy_list[i])
    
    if L_protoSSC_LTE[i] == 0:
        L_proto_LTE = 0.0001
    else: 
        L_proto_LTE = L_protoSSC_LTE[i]
    if L_nLTE[i] == 0:
        L_n = 0.0001
    else:
        L_n = L_nLTE[i]
    L_msSSC_L_protoSSC_ratio.append(L_ir_leroy_list[i]/L_proto_LTE)
    T_cond.append(Tvib_LTE[i]*np.sqrt(sizes_pc_yo[i]/leroy_sizes_pc[i]))
    Tcond_nLTE.append(Tex_nLTe[i]*np.sqrt(sizes_pc_yo[i]/leroy_sizes_pc[i]))

    lsum = ((L_ir_leroy_list[i]+L_protoSSC_LTE[i])*u.Lsun).to(u.W)
    lsum_nLTE = ((L_ir_leroy_list[i]+L_nLTE[i])*u.Lsun).to(u.W)
    
    rad_m = (leroy_sizes_pc[i] * u.pc).to(u.m)/2. # tamaño leroy
    Tcond_2.append(((lsum/(4.*np.pi*s_si*(rad_m**2)))**(1./4)).value)
    Tcond_nLTE_2.append(((lsum_nLTE/(4.*np.pi*s_si*(rad_m**2)))**(1./4)).value)
    rad_m = (sizes_pc_yo[i] * u.pc).to(u.m)/2. # tamaño yo
    Tcond_3.append((lsum/(4.*np.pi*s_si*(rad_m**2)))**(1./4))
    
    # Tamaños a partir de la luminosidad LSHC + L_MS
    sizes_Ltot.append(np.sqrt(lsum/(4.*np.pi*s_si*((u.K*T_cond[i])**4.))).to(u.pc).value)
    sizes_Ltot_nLTE.append(np.sqrt(lsum_nLTE/(4.*np.pi*s_si*((u.K*Tcond_nLTE[i])**4.))).to(u.pc).value)
    
    ratio_LTE.append(L_proto_LTE/L_ir_leroy_list[i])
    ratio_nLTE.append(L_ir_leroy_list[i]/L_n)
#T = u.Quantity(T, u.K)
#r = u.Quantity(r, u.m)
#L = 4*pi*(r**2)*s_si*T**4
    
L_ir_percen_leroy = 100.*L_ir_leroy/l_ir_nuc



#### RMSss ####
rms_v724_23 = 1.5
rms_v739_38 = 1.1
resol_channel = 5

rms_df = pd.DataFrame({
    'v7_24_23_intens' : [7.47, 23.8999, 62.7569, 6.155],
    'v7_39_38_intens' : [20.55, 22.1407, 103.609, 23.047],
    'v7_24_23_vels' : [162-196, 256-297, 112-165, 139-165],
    'v7_39_38_vels' : [756-720, 681-601, 766-729, 770-715]
    })
rms_df['24_23_ncanales'] =    np.abs(rms_df['v7_24_23_vels']/resol_channel)
rms_df['39_38_ncanales'] =    np.abs(rms_df['v7_39_38_vels']/resol_channel)

rms_df['24_23_sigma'] = rms_v724_23/np.sqrt(rms_df['24_23_ncanales'])
rms_df['39_38_sigma'] = rms_v739_38/np.sqrt(rms_df['39_38_ncanales'])    

rms_df['24_23_3sigmaAv'] =  3*rms_df['24_23_sigma']*np.abs(rms_df['v7_24_23_vels'])
rms_df['39_38_3sigmaAv'] =  3*rms_df['39_38_sigma']*np.abs(rms_df['v7_39_38_vels'])

# LTE MADCUBA
leroy_df['LTE_NHC3N'] = (10**hc3n_prop['logN_vib']).tolist()
leroy_df['Tvib'] = Tvib_LTE
leroy_df['LTE_Ltot'] = L_protoSSC_LTE
leroy_df['LTE_Mproto'] = leroy_df['LTE_Ltot']/1000 # Lum-to-mass ratio = 1000

# Leroy
leroy_df['MS_Ltot'] = L_ir_leroy_list
leroy_df['dcond_nH'] = nh2_leroy_list
leroy_df['dcond_NH'] = Nh2_leroy_list
leroy_df['dcond_sizes_pc'] = leroy_sizes_pc
leroy_df['GasMass'] = leroy_gas_mass
leroy_df['StarMass'] = leroy_star_mass
leroy_df['dcond_sizes_pc'] =leroy_sizes_pc

# Both
leroy_df['XHC3N_LTE_leroy'] = Xhc3n_list
leroy_df['XHC3N_nLTE_leroy'] = Xhc3n_list_nLTE
leroy_df['ratio_LTE'] = ratio_LTE
leroy_df['ratio_nLTE'] = ratio_nLTE
leroy_df['Tcond_LTE'] = T_cond
leroy_df['Tcond_nLTE'] = Tcond_nLTE
leroy_df['Tcond_LTE_Lsum'] = Tcond_2
leroy_df['Tcond_nLTE_Lsum'] = Tcond_nLTE_2
leroy_df['Sizes_pc_LTE'] = sizes_Ltot
leroy_df['Sizes_pc_nLTE'] = sizes_Ltot_nLTE

leroy_df['LTE_Ltot_using'] = [1.132087, 5.242486, 7.866067, 5.289255,3.112957, 0.0, 0.0, 1.782151, 0.022827, 0.012235, 0.087843, 0.008873, 22.7424252, 21.68234234] 
leroy_df['LTE_Ltot_using'] = leroy_df['LTE_Ltot_using'] * 1e8
leroy_df['LTE_Ltot_using_1e8'] = [1.132087, 5.242486, 7.866067, 5.289255,3.112957, 0.0, 0.0, 1.782151, 0.022827, 0.012235, 0.087843, 0.008873, 22.7424252, 21.68234234] 


leroy_df['LTE_Ltot_using_surface '] = leroy_df['LTE_Ltot_using']/(np.pi * ((1.7/2)**2))

leroy_df['Size_pc'].loc[10] = 0.24
leroy_df['Size_pc'].loc[12] = 0.25
leroy_df['Size_pc'].loc[6] = 1.9
leroy_df['Size_pc'].loc[7] = 1.9

leroy_df['Mass_msun'].loc[6] = 0
leroy_df['Mass_msun'].loc[7] = 0
leroy_df['Mass_msun'].loc[10] = 0
leroy_df['Mass_msun'].loc[12] = 0
sizes_36ghz = 1.9 # pc
leroy_df['LTE_Ltotsup_using'] = leroy_df['LTE_Ltot_using']/(np.pi * ((sizes_36ghz/2)**2))
leroy_df['LTE_Ltotsup_sizeme_using'] = leroy_df['LTE_Ltot_using']/(np.pi * ((leroy_df['Size_pc']/2)**2))


# Tan et al. 2014 Massive star Formation review
leroy_df['LTE_Mgas_surface_Msunpc2'] = leroy_df['GasMass']/(np.pi * ((leroy_df['Size_pc']/2)**2)) 
leroy_df['LTE_Mgas_surface_gcm2'] = leroy_df['LTE_Mgas_surface_Msunpc2']*((1 * u.Msun / (u.pc**2)).to(u.g/(u.cm**2))).value
leroy_df['sigma'] = hc3n_prop['fwhm'].tolist()


leroy_df['alfaVIR'] = 5. * (((1 * u.km *2. *np.sqrt( 2.*np.log(2))*leroy_df['sigma']).to(u.m)/u.s)**2) * ((1* u.pc*leroy_df['Size_pc'] )/2).to(u.m)/(_si.G*(1*u.Msun *leroy_df['GasMass']).to(u.kg))
leroy_df['vesc_kms'] = 2. *np.sqrt( 2.*np.log(2))*leroy_df['sigma']*(10./leroy_df['alfaVIR'])**(1./2)
leroy_df['Pmean'] = (3.*np.pi/20)*(leroy_df['GasMass']/leroy_df['StarMass'])*leroy_df['alfaVIR'] * _si.G * ((1*u.g / u.cm**2)*leroy_df['LTE_Mgas_surface_gcm2'] ).to(u.kg/u.m**2)**2


leroy_df['dens1'] = ((1. *(u.Msun/(u.pc**3)) * leroy_df['GasMass']/((4./3)*np.pi*(leroy_df['Size_pc']/2)**3))).to(u.kg/u.m**3)
leroy_df['dens3'] = ((1. *(u.Msun/(u.pc**3)) * leroy_df['GasMass']/((4./3)*np.pi*(1.9/2)**3))).to(u.kg/u.m**3)

leroy_df['nH_kgm3'] = leroy_df['nH'] * 2.*(_si.m_p + _si.m_e) * (u.cm**3).to(u.m**3)
leroy_df['tff_mean1_yr'] =  (3.*np.pi/(32.*_si.G*leroy_df['dens1']*u.kg/u.m**3)**(1./2)).to(u.yr)
leroy_df['tff_mean3_yr'] =  (3.*np.pi/(32.*_si.G*leroy_df['dens3']*u.kg/(u.m**3))**(1./2)).to(u.yr)
leroy_df['tff_mean2_yr'] = 6.85E4  *((leroy_df['GasMass']/1000.)**(1./4)) * (leroy_df['LTE_Mgas_surface_gcm2']**(-3./4))


leroy_df['acc_rate1']= leroy_df['GasMass']/leroy_df['tff_mean1_yr']
leroy_df['acc_rate2']= leroy_df['GasMass']/leroy_df['tff_mean2_yr']
leroy_df['acc_rate3']= leroy_df['GasMass']/leroy_df['tff_mean3_yr']

# Distance to the center (i.e. TH2)
# Positions (me and Leroy2018)
hc3n_positions = pd.read_csv(dworkdir_spec+'/HC_positions_REAL.txt', delim_whitespace= True, header=0, comment='#')
         
from astropy.coordinates import SkyCoord                            
# Center 
TH2_RA = '00:47:33.179'
TH2_Dec = '-25:17:17.13'
th2_pos  = utiles.HMS2deg(ra=TH2_RA.replace(':', ' '), dec=TH2_Dec.replace(':', ' '))
TH2_RA_s= TH2_RA.replace(':', 'h', 1).replace(':', 'm', 1) + 's'
TH2_Dec_s= TH2_Dec.replace(':', 'd', 1).replace(':', 'm', 1) + 's'

th2_sky = SkyCoord(TH2_RA_s, TH2_Dec, frame='fk5', unit=(u.hourangle, u.deg), distance = D*u.Mpc)
# Kinematic center
nucleus_RA2=  '00:47:33.14' #http://iopscience.iop.org/article/10.1088/0004-637X/716/2/1166/pdf
nucleus_Dec2= '-25:17:17.52'

nucleus_RA2_s= nucleus_RA2.replace(':', 'h', 1).replace(':', 'm', 1) + 's'
nucleus_Dec2_s= nucleus_Dec2.replace(':', 'd', 1).replace(':', 'm', 1) + 's'

kin_pos  = utiles.HMS2deg(ra=nucleus_RA2.replace(':', ' '), dec=nucleus_Dec2.replace(':', ' '))
nuc_sky = SkyCoord(nucleus_RA2_s, nucleus_Dec2, frame='fk5', unit=(u.hourangle, u.deg), distance =D*u.Mpc)

dist_cent_th2 = []
dist_cent_kin = []
dist = []
dist_d = []
ang_dist_th2 = []
dist_k = []

for i, line in hc3n_positions.iterrows():
        pos = utiles.HMS2deg(ra=line['RA_yo'].replace(':', ' '), dec=line['Dec_yo'].replace(':', ' '))

        pos_ra_s= line['RA_yo'].replace(':', 'h', 1).replace(':', 'm', 1) + 's'
        pos_Dec_s= line['Dec_yo'].replace(':', 'd', 1).replace(':', 'm', 1) + 's'
        pos_sky = SkyCoord(pos_ra_s, pos_Dec_s, frame='fk5', unit=(u.hourangle, u.deg), distance = 3.5*u.Mpc)

        ang_dist = utiles.ang_distance_btw_2points_v2(float(pos[0]),float(pos[1]),float(th2_pos[0]),float(th2_pos[1]))
        lin_dist = u_conversion.lin_size(D, ang_dist*60*60.)
        dist_cent_th2.append(lin_dist.to(u.pc).value)
        ang_dist_th2.append(ang_dist)
        
        #a_d = np.sqrt((float(pos[0])- float(th2_pos[0]))**2 + (float(pos[1])- float(th2_pos[1]))**2)
        #l_d = u_conversion.lin_size(D, a_d*60.*60).to(u.pc)
        if i >=12:
            signo = -1
        else:
            signo = +1
        dist.append(signo*th2_sky.separation(pos_sky))

        dist_d.append(signo*th2_sky.separation_3d(pos_sky).to(u.pc).value)
        
        

        ang_dist = utiles.ang_distance_btw_2points_v2(float(pos[0]),float(pos[1]),float(kin_pos[0]),float(kin_pos[1]))
        lin_dist = u_conversion.lin_size(D, ang_dist*60.*60.)
        dist_cent_kin.append(lin_dist.to(u.pc).value)
        dist_k.append(signo*nuc_sky.separation_3d(pos_sky).to(u.pc).value)
        
leroy_df['dist_th2_pc'] = dist_d       
        

# L IC 860  Costagliola2011 et al
Tv7ic860 = 42
Tv0ic860 = 15
size = 29 #arcsec
dist_860 = 90.1 # Mpc
size_pc_860 = u_conversion.lin_size(dist_860, size).to(u.pc)
Lic860 = u_conversion.stef_boltz(u_conversion.lin_size(dist_860, size)/2,Tv7ic860).to(u.Lsun)
Lic860_sup = Lic860/(np.pi * ((size_pc_860/2)**2))

## Plotting 
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
from scipy.ndimage import rotate
# Cube workdir
dworkdir_cubes = '/Users/frico/Documents/data/NGC253_H3O+'
dworkdir_spec = dworkdir_cubes+'/Hotcores_v4_all'
# Out dir
out_dir = dworkdir_spec+'/Results_v8/figura3/'
if not os.path.exists(out_dir):
        os.makedirs(out_dir)
# Cont at 218GHz
cont_fits_name = dworkdir_spec+'/MAD_CUB_CROP_NGC253_TE_cont_218_HR_briggs.pbcor.fits'
ngc253_cont_fits = fits.open(cont_fits_name)#'MAD_CUB_NGC253_TE_cont_358_natural.pbcor.fits')
ngc253_cont_data = ngc253_cont_fits[0]
ngc253_cont_shape = ngc253_cont_data.data.shape
ngc253_cont_header = ngc253_cont_data.header
ngc253_cont_header['CDELT4'] = 1
cont_stdev  = np.nanstd(ngc253_cont_data.data[0,0,:,:])
cont_max    = np.nanmax(ngc253_cont_data.data[0,0,:,:])
cont_min    = np.nanmin(ngc253_cont_data.data[0,0,:,:])
ngc253_cont_header['CTYPE3'] = 'FREQ' # En el header aparace 'FREQ-TOP'
wcs_1 = WCS(ngc253_cont_header)
wcs_1.wcs.ctype = ['RA---SIN', 'DEC--SIN', 'FREQ', 'STOKES']
wcs_1 = wcs_1.dropaxis(3)
wcs_1 = wcs_1.dropaxis(2)       

pa = 47. 

#from mpdaf.obj import Image, WCS
#ima = Image(data=ngc253_cont_data.data[0,0,:,:], wcs=wcs_1)

def rotate2(degs, header):
        """Return a rotation matrix for counterclockwise rotation by ``deg`` degrees."""
        rads = np.radians(degs)
        s = np.sin(rads)
        c = np.cos(rads)
        return np.array([[c*header['CDELT1'], -s*header['CDELT2']],
                      [s*header['CDELT1'], c*header['CDELT2']]])
                           
cd_matrix = rotate2(pa-90., ngc253_cont_header)  


#w.wcs.set_pv([(2, 1, 45.0)])
#
#    
#from ferpy import fits_header_mod
#


cd11 = cd_matrix[0][0]
cd21 = cd_matrix[1][0]
cd22 = cd_matrix[1][1]
pixsize = 0.03/60./60 #deg

ra0 = 11.88835307 #261
dec0 = -25.29169382 #-69
rot = rotate(ngc253_cont_data.data[0,0,:,:], pa-90., reshape=False)
H2_RA = '00:47:33.179'
TH2_Dec = '-25:17:17.13'
# Convert pixel coordinates to world coordinates
#center = w.wcs_pix2world( len(rot[0])/2, len(rot)/2,  1)
#ang_dist_th2cent = utiles.ang_distance_btw_2points_v2(center[0],center[1],float(th2_pos[0]),float(th2_pos[1]))
#ang_dist_th2cent_pc = u_conversion.lin_size(D, ang_dist_th2cent*60.*60.).to(u.pc).value



center_px = len(rot[0])/2
center_py = len(rot)/2

#center_px = -59+len(rot[0])/2 ## Kin center
center_px = -81+len(rot[0])/2 ## Th2 center

origin_px = 0
origin_py = 0

last_px = len(rot[0])
last_py = len(rot)



new_xaxis = []
new_xaxis2 = []
new_yaxis = [] 
new_yaxis2 = []
for i, px in enumerate(rot[0]):
    dif_to_cen_px = (i)-center_px
    new_xaxis.append(dif_to_cen_px*u_conversion.lin_size(D, -cd11*60.*60.).to(u.pc).value)
    new_xaxis2.append(dif_to_cen_px*u_conversion.lin_size(D, 0.03).to(u.pc).value)
for i, pyind in enumerate(rot):
    dif_to_cen_py = (i)-center_py
    new_yaxis.append(dif_to_cen_py*u_conversion.lin_size(D, cd21*60.*60.).to(u.pc).value)
    new_yaxis2.append(dif_to_cen_py*u_conversion.lin_size(D, 0.03).to(u.pc).value)

    

origin_x = wcs_1.wcs_pix2world(0, 0,  1)
last_x = wcs_1.wcs_pix2world(0, len(rot[0]),  1)
last_y = wcs_1.wcs_pix2world(0, len(rot),  1)
ang_dist_th2orig = utiles.ang_distance_btw_2points_v2(origin_x[0],origin_x[1],float(th2_pos[0]),float(th2_pos[1]))
ang_dist_th2last = utiles.ang_distance_btw_2points_v2(last_x[0],last_x[1],float(th2_pos[0]),float(th2_pos[1]))

ang_x = utiles.ang_distance_btw_2points_v2(last_x[0],last_x[1],float(origin_x[0]),float(origin_x[1]))


dec = dec0 + cd22
ra = ra0 - cd11/np.cos(np.radians(dec))


w = WCS(naxis=2)

w.wcs.crpix = [len(rot)/2, len(rot[0])/2]
w.wcs.cdelt = [cd_matrix[0][0], -cd_matrix[0][0]]
w.wcs.crval = [11.887901459968038, -25.28834381452057]
w.wcs.ctype = ['RA---SIN', 'DEC--SIN']
TH2w_px = w.wcs_world2pix(float(th2_pos[0]), float(th2_pos[1]), 1)


header = w.to_header()
hdu = fits.PrimaryHDU(rot, header=header)
hdu_data = hdu.data
hdu_header = hdu.header
wcs_hdu = WCS(hdu_header)

TH2_RA = '00:47:33.179'
TH2_Dec = '-25:17:17.13'
# Convert pixel coordinates to world coordinates
TH2_px = wcs_1.wcs_world2pix(float(th2_pos[0]), float(th2_pos[1]), 1)
# Convert the same coordinates back to pixel coordinates.
pixcrd2 = wcs_1.wcs_pix2world(len(rot)/2, len(rot[0])/2, 1)

origin  = wcs_1.wcs_pix2world(0, 0,  1)
origin_ra = 1 * cd11 
origin_dec = 1 * cd21 

pxs = []
for i, ang in enumerate(ang_dist_th2):
    pxs.append(-ang/cd11)

plot_cont = True
if plot_cont == True:
    
    # Sorting df by distance
    leroy_df_sorted = leroy_df.sort_values(by=['dist_th2_pc'])
    #fig = plt.figure()
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, figsize=(10,20))#, sharex=True)
    #ax1 = fig.add_subplot((111), aspect='equal')#, projection=wcs)
    ax1.tick_params(labelsize=6)
        
    cont_stdev  = np.nanstd(ngc253_cont_data.data[0,0,:,:])
    cont_max    = np.nanmax(ngc253_cont_data.data[0,0,:,:])
    cont_min    = np.nanmin(ngc253_cont_data.data[0,0,:,:])
    
    axes = [ax1, ax2, ax3, ax4, ax5]
    
    pixsize = 0.03 #arcsec
    ell_ax1 = 0.292 #arcsec
    ell_ax2 = 0.197 # arcsec
    pa_e = 80.279
    r_ax1 = 0.35 
    r_ax2 = 0.35
    
    
    #rot = rotate_yo(wcs_1, pa-90)
    dx = (new_xaxis[1]-new_xaxis[2])/2.
    dy = (new_yaxis[1]-new_yaxis[2])/2.
    extent = [new_xaxis2[0], new_xaxis2[-1], new_yaxis2[0], new_yaxis2[-1]]

    ax1.imshow(rot, origin='lower', vmax=cont_max, vmin=0.2*cont_stdev, cmap= cm.jet, interpolation="none", extent=extent,  aspect="auto")#, extent=(0,0,0,0))
    
    ax1.set_ylim([-15, 15])#[155, 482]) #336, 453
    
    # Plotting the center
    pos_nucleus = utiles.HMS2deg(ra=nucleus_RA2.replace(':', ' '), dec=nucleus_Dec2.replace(':', ' '))
    px_n, py_n = wcs_hdu.wcs_world2pix(float(pos_nucleus[0]), float(pos_nucleus[1]), 1)
    ax1.plot(px_n, py_n, 'x', color='white')
    
    #ax1.plot(TH2w_px[0], TH2w_px[1], 'x', color='white')
    ax1.plot(0, 0, '+', color='yellow')
    # Plotting indices
    px_list_shc = [110, 106, 103, 77, 62.5, 54.59, 51.55, 17.65, 23.29, 4, 7, 1.267, -7.5, -37.59]
    py_list_shc = [8, 0.5, 7.5, 1, 1.75, -8, 3.8, 7, -8, 7, -8, -6, 5, 5]
    
    
    lproto_uplims = []
    for i, ind in enumerate(dist_d):
        if i in [1,5,8,10,11]:
            py_shc=-8
        else:
            py_shc=+8
        if i in [5,6,8,9,10,11]:
            lproto_uplims.append(1)
        else:
            lproto_uplims.append(0)
        ax1.annotate(str(i+1), xy=(dist_d[i],py_list_shc[i]), xytext=(dist_d[i],py_list_shc[i]),
                      va='center', ha='center', color = 'white')
                #arrowprops={'arrowstyle': '-', 'color': 'w'},)
    
    

    ax1.set_ylabel(r'$\Delta$y (pc)', fontsize=11, labelpad=1)
    
    #lproto_plot_dist =leroy_df_sorted['dist_th2_pc'].tolist()
    #del lproto_plot_dist[5:6]
    
    # plotting lines:
    # Lproto
    ax2.plot(leroy_df_sorted['dist_th2_pc'],leroy_df_sorted['LTE_Ltot']/1E8,
             c='k', marker='', linestyle='--', linewidth=0.6)
    # Lproto/L*
    ax3.plot(leroy_df_sorted['dist_th2_pc'],leroy_df_sorted['ratio_LTE'],
             c='k', marker='', linestyle='--', linewidth=0.6)
    # Mproto+M*
    ax4.plot(leroy_df_sorted['dist_th2_pc'],(leroy_df_sorted['StarMass']+leroy_df_sorted['LTE_Mproto'])/1E4,
             c='k', marker='', linestyle='--', linewidth=0.6)
    # Mproto+M*/Mgas
    ax5.plot(leroy_df_sorted['dist_th2_pc'],(leroy_df_sorted['StarMass']+leroy_df_sorted['LTE_Mproto'])/leroy_df_sorted['GasMass'],
             c='k', marker='', linestyle='--', linewidth=0.6)

    
    for l, line in leroy_df_sorted.iterrows():
        # upperlimits
        if l in [6,7,9,10,11,12]:
            if l not in [6,7]:
                # Lproto
                ax2.errorbar(line['dist_th2_pc'], line['LTE_Ltot']/1E8,
                             uplims=True,
                             yerr=2,
                             #lolims=lowerlimits,
                             capsize=2,
                             marker='', color='k',
                             markeredgecolor='k', markerfacecolor='k',
                             linewidth=0.6, linestyle='--', alpha=1)
                # Lproto/L*
                ax3.errorbar(line['dist_th2_pc'], line['ratio_LTE'],
                             uplims=True,
                             yerr=4,
                             #lolims=lowerlimits,
                             capsize=2,
                             marker='', color='k',
                             markeredgecolor='k', markerfacecolor='k',
                             linewidth=0.6, linestyle='--', alpha=1)
            # Lproto/L*
            ax4.errorbar(line['dist_th2_pc'], (line['StarMass']+line['LTE_Mproto'])/1E4,
                             uplims=True,
                             yerr=15,
                             #lolims=lowerlimits,
                             capsize=2,
                             marker='', color='k',
                             markeredgecolor='k', markerfacecolor='k',
                             linewidth=0.6, linestyle='--', alpha=1)
            # Mproto+M*/Mgas
            ax5.errorbar(line['dist_th2_pc'], (line['StarMass']+line['LTE_Mproto'])/line['GasMass'],
                             uplims=True,
                             yerr=5,
                             #lolims=lowerlimits,
                             capsize=2,
                             marker='', color='k',
                             markeredgecolor='k', markerfacecolor='k',
                             linewidth=0.6, linestyle='--', alpha=1)
        else:
            ax2.plot(line['dist_th2_pc'],line['LTE_Ltot']/1E8, c='k', marker='.',
                     markersize=6, linestyle='')
            ax3.plot(line['dist_th2_pc'],line['ratio_LTE'],
                     c='k', marker='.', markersize=6, linestyle='')
            ax4.plot(line['dist_th2_pc'],(line['StarMass']+line['LTE_Mproto'])/1E4,
             c='k', marker='.', linestyle='')
            ax5.plot(line['dist_th2_pc'],(line['StarMass']+line['LTE_Mproto'])/line['GasMass'],
             c='k', marker='.', linestyle='')
            
    
    #ax2.plot(leroy_df_sorted['dist_th2_pc'],leroy_df_sorted['LTE_Ltot']/1E8, c='k', marker='.', markersize=8, linestyle='--', linewidth=0.6)
    ax2.set_ylabel(r'L$_{\rm{p}*}$($\times 10^8$L$_\odot$)', fontsize=11, labelpad=1)
    ax2.set_ylim([-5, 30])
    # Lproto/L*
    ax3.set_ylabel(r'L$_{\rm{p}*}/$L$_*$', fontsize=11, labelpad=1)
    ax3.set_ylim([-10, 75])
    # Mproto+M*
    ax4.set_ylabel(r'M$_{\rm{p}*}+$M$_{*}$ ($\times 10^4$M$_\odot$)', fontsize=11, labelpad=3)
    ax4.set_ylim([-25, 300])
    # Mproto+M*/Mgas
    ax5.set_ylabel(r'M$_{\rm{p}*}+$M$_{*}/$M$_{\rm{gas}}$', fontsize=11, labelpad=3)
    ax5.set_xlabel(r'$\Delta$x (pc)', fontsize=12)
    ax5.set_ylim([-10, 90])
    for i, axis in enumerate(axes):
        axis.set_xlim([-50, 130])
        axis.tick_params(labelsize=12)
        #axis.tick_params(direction='in', labelsize=3)
        axis.tick_params(axis='both', which='both', direction='in', width=1.0)
        axis.minorticks_on()
        axis.xaxis.set_tick_params(which='both', top ='on',labeltop='off')
        axis.yaxis.set_tick_params(which='both', right='on', labelright='off')
        ylim = axis.get_ylim()
        axis.get_yaxis().set_label_coords(-0.05,0.5)#np.abs(ylim[0]-ylim[1])/2)
        if i != len(axes)-1 and i != 0:
            axis.set_xticklabels([])
            
        if i == 0:
            axis.tick_params(axis='both', which='both', direction='in', width=1.0, color='white')
            #axis.xaxis.set_tick_params(which='both', top ='on', labeltop='on')
            axis.xaxis.set_tick_params(labeltop='on', labelbottom='off')

        
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.savefig(out_dir+'ngc253_cont218_subplot_v'+version+'.eps', bbox_inches='tight', transparent=True, dpi=300)
    plt.show() 

    
    
table_leroy_out = ['Tdust', 'Tvib', 'Tcond_LTE', 'Tcond_nLTE', 'Tcond_LTE_Lsum', 'Tcond_nLTE_Lsum',
                   'nH', 'dcond_nH', 
                   'Size_asec', 'Size_pc', 'dcond_sizes_pc',
                   'new_NH', 'dcond_NH', 'new_NHC3N', 'LTE_NHC3N',
                   'Xline', 'XHC3N_LTE_leroy', 'XHC3N_nLTE_leroy',
                   'Mass_msun','GasMass', 'StarMass',
                   'SHC_Ltot', 'LTE_Ltot', 'MS_Ltot',
                   'ratio_LTE', 'ratio_nLTE']

leroy_df.to_csv('/Users/frico/Documents/data/NGC253_H3O+/'+'Hotcores_v4_all/model_results_and_Leroy_v2.txt', columns=table_leroy_out, header=True, sep='\t', index=True, na_rep=-1)


#v7_1_df.to_csv('/Users/frico/Documents/data/NGC253_H3O+/'+'Hotcores_v4_all/v7_1_parameters.txt', header=True, sep='\t', index=True, na_rep=-1)
#v7_1_df.loc[(v7_1_df['nup'] == 24) & (v7_1_df['nlow'] == 23)]
#v7_1_df.loc[(v7_1_df['nup'] == 39) & (v7_1_df['nlow'] == 38)]

#df_highT2.to_latex(buf=None, columns=table_out)    

# Plotting line profiles:
plot_profiles = False
if plot_profiles == True:
    figoutdir = '/Users/frico/Documents/Ed/modelos_hc3n/v'+version+'_tcte_T/r_0.0/Figures/'
    out_dir_prof = figoutdir+'/'+'Profiles'
    if not os.path.exists(out_dir_prof):
        os.makedirs(out_dir_prof)
    out_dir_sed = figoutdir+'/'+'SED'
    if not os.path.exists(out_dir_sed):
        os.makedirs(out_dir_sed)
        
    # Initializing figure
    fig_PROF, axis_PROF = plot_PROF_ini(obs_df=hc3n_prop, velocity_limits = ['NA', 'NA'], intensity_limits = ['NA', 'NA'])
    fig_SED, axis_SED = plot_SED_ini()
    
    # For each selected model:
    
    # Plotting profiles and sed
    for i, num in enumerate([1,2,3,4,5,8,9,11,13,14]):
        print 'plotting'
        modeldf = HC_mods_ok_highT['HC'+str(num)]
        fig_PROF, fig_SED, SED_ISO, SED_SPIRE = luminosity_calculator_allinone(version, modeldf, fig_PROF, fig_SED, axis_PROF, axis_SED)
        if num == 1:
            SED_ISO_tot = np.array(SED_ISO['Ic_Jy'])
            SED_SPIRE_tot = np.array(SED_SPIRE['Ic_Jy'])
        else:
            SED_ISO_tot = SED_ISO_tot + np.array(SED_ISO['Ic_Jy'])
            SED_SPIRE_tot = SED_SPIRE_tot + np.array(SED_SPIRE['Ic_Jy'])
            
    axis_SED[0].plot(SED_ISO['lambda_um'], SED_ISO_tot, linewidth=0.7, color='k')
    axis_SED[0].plot(SED_SPIRE['lambda_um'], SED_SPIRE_tot, linewidth=0.7, color='k')
    
    # Saving figures
    fig_PROF.savefig(out_dir_prof+'/r_00_Tcte_profile.'+fig_format, bbox_inches='tight', transparent=True, dpi=100, format=fig_format)
    plt.close(fig_PROF)
    
    fig_SED.savefig(out_dir_sed+'/r_00_Tcte_SED.'+fig_format, bbox_inches='tight', transparent=True, dpi=100, format=fig_format)
    plt.close(fig_SED)
    
    # Initializing figure
    fig_PROF, axis_PROF = plot_PROF_ini(obs_df=hc3n_prop, velocity_limits = ['NA', 'NA'], intensity_limits = ['NA', 'NA'])
    fig_SED, axis_SED = plot_SED_ini()
    
    # Low temperatures
    # Plotting profiles and sed
    for i, num in enumerate([1,2,3,4,5,8,9,11,13,14]):
        print 'plotting'
        modeldf = HC_mods_ok_lowT['HC'+str(num)]
        fig_PROF, fig_SED, SED_ISO, SED_SPIRE = luminosity_calculator_allinone(version, modeldf, fig_PROF, fig_SED, axis_PROF, axis_SED)
        if num == 1:
            SED_ISO_tot = np.array(SED_ISO['Ic_Jy'])
            SED_SPIRE_tot = np.array(SED_SPIRE['Ic_Jy'])
        else:
            SED_ISO_tot = SED_ISO_tot + np.array(SED_ISO['Ic_Jy'])
            SED_SPIRE_tot = SED_SPIRE_tot + np.array(SED_SPIRE['Ic_Jy'])
            
    axis_SED[0].plot(SED_ISO['lambda_um'], SED_ISO_tot, linewidth=0.7, color='k')
    axis_SED[0].plot(SED_SPIRE['lambda_um'], SED_SPIRE_tot, linewidth=0.7, color='k')
    
    # Saving figures
    fig_PROF.savefig(out_dir_prof+'/r_00_Tctelow_profile.'+fig_format, bbox_inches='tight', transparent=True, dpi=300, format=fig_format)
    plt.close(fig_PROF)
    
    fig_SED.savefig(out_dir_sed+'/r_00_Tctelow_SED.'+fig_format, bbox_inches='tight', transparent=True, dpi=300, format=fig_format)
    plt.close(fig_SED)
    

print 'Llow=\t'+'%1.3E' % Ltotal_lowt
print 'Lhig=\t'+'%1.3E' % Ltotal_hight

L_IR_TOTAL_GA15 = 3.2E10
print 'L-LTE  %=\t'+'%1.3E' % (hc3n_prop['L_Lsun_v7_219_g'].sum()*100/L_IR_TOTAL_GA15)
print 'L-NLTE %=\t'+'%1.3E' % (Ltotal_hight*100/L_IR_TOTAL_GA15)

hc3n_prop[['Ind_ok', 'Tvib', 'L_Lsun_v7_219_g', 'L_Lsun_v7_219_err_g']]
    
print s_unique_dens
print mod_df.central_T.unique()