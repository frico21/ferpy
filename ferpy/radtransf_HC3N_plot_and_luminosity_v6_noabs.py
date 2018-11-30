#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 11:47:04 2018

@author: frico
"""

# =============================================================================
# Plot Model Data
# =============================================================================


import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import pandas as pd


from ferpy import u_conversion
from ferpy import utiles
import astropy.units as u
import gc

"""
"""

# Plotting
#plot_profiles = True
#plot_SED = False

# Model paths
#model_path = '/Users/frico/Documents/Ed/modelos_hc3n/modelos2'
#model_name = 'mhc3n_2_tgas350_tdust350'

# Model Factor
#model_factor = 0.7#0.1854 #0.7

# Fig out dir
#figoutdir = model_path+'/figures'
        
# Tamano minimo 0.039" 
tam_min_pc = 0.039*17/2
factor_tam_min = tam_min_pc**2

# Tamano max 0.1" 
tam_max_pc = 0.1*17/2
factor_tam_max = tam_max_pc**2
        

# Observed values HC14
hc3n_v0_24_23_obs = 54.08 #mJy/beam
hc3n_v0_39_38_obs = 70.29 #mJy
hc3n_v0_ratio_obs = hc3n_v0_39_38_obs/hc3n_v0_24_23_obs

hc3n_v7_24_23_obs = 22.05 #mJy
hc3n_v7_39_38_obs = 44.59 #mJy
hc3n_v7_ratio_obs = hc3n_v7_39_38_obs/hc3n_v7_24_23_obs

def luminosity_calculator(model_path, model_name, factor, plot_profiles, plot_SED, figoutdir, fig_format):
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
    # Multiplying intensities by factor^2
    model_prof['I_Ic_mJy'] = model_prof['I_Ic_mJy']*(factor**2.)
    
    
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
    # To avoid ratios with absorption
    if prof_v0_24_23['I_Ic_mJy'][prof_v0_24_23.V_kms == 0].max() == prof_v0_24_23['I_Ic_mJy'].max():
        peak_v0_24_23 = prof_v0_24_23['I_Ic_mJy'][prof_v0_24_23.V_kms == 0].max()
    else: 
        peak_v0_24_23 = 0
                
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
    # To avoid ratios with absorption
    if prof_v7_24_23['I_Ic_mJy'][prof_v7_24_23.V_kms == 0].max() == prof_v7_24_23['I_Ic_mJy'].max():
        peak_v7_24_23 = prof_v7_24_23['I_Ic_mJy'][prof_v7_24_23.V_kms == 0].max()
    else: 
        peak_v7_24_23 = 0
        
        
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
    # To avoid ratios with absorption
    if prof_v0_39_38['I_Ic_mJy'][prof_v0_39_38.V_kms == 0].max() == prof_v0_39_38['I_Ic_mJy'].max():
        peak_v0_39_38 = prof_v0_39_38['I_Ic_mJy'][prof_v0_39_38.V_kms == 0].max()
    else:
        peak_v0_39_38 = 0
    
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
    # To avoid ratios with absorption
    if prof_v7_39_38['I_Ic_mJy'][prof_v7_39_38.V_kms == 0].max() == prof_v7_39_38['I_Ic_mJy'].max():
        peak_v7_39_38 = prof_v7_39_38['I_Ic_mJy'][prof_v7_39_38.V_kms == 0].max()
    else:
        peak_v7_39_38 = 0
        
    
    # Modeled values 
    hc3n_v0_ratio = peak_v0_39_38/peak_v0_24_23
    hc3n_v7_ratio = peak_v7_39_38/peak_v7_24_23
    hc3n_v0v7_3938_ratio = peak_v0_39_38/peak_v7_39_38
    hc3n_v0v7_2423_ratio = peak_v0_24_23/peak_v7_24_23
    
    # Plotting line profiles
    if plot_profiles == True:
        m = filas
        n = columnas
        dx, dy = 1.2, 1
        figsize = plt.figaspect(float(dy * m) / float(dx * n))
        fig = plt.figure(figsize=figsize)
        gs1 = gridspec.GridSpec(m, n)    
        gs1.update(wspace = 0.0, hspace=0.0, top=0.95, bottom = 0.05)
        # Generating specified number of axis
        axis = []
        axis_ind = []
        for i in range(m*n):
            axis.append(fig.add_subplot(gs1[i]))
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
                                            prof_v7_39_38['I_Ic_mJy'].max()])]
        
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
                ax.set_ylim([intensity_limits_total[0], intensity_limits_total[1]+intensity_limits_total[1]*0.4]) # Para que el maximo en el eje y sea el mismo en todos los espectros
            
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
    
        fig.savefig(figoutdir+'/'+model_name+'_f'+'%1.1f' % factor +'_profile.'+fig_format, bbox_inches='tight', transparent=True, dpi=400, format=fig_format)
        plt.close()
            
    
    
    
    # =============================================================================
    # ## SED and tau_continuum
    # =============================================================================
        
    
    ## ISO
    model_iso = pd.read_csv(model_path+'/'+model_name+'_.iso', delim_whitespace= True)
    model_iso.columns = ['lambda_um', 'I_Ic_ergscm2um', 'Ic_ergscm2um', 'Inorm_ergscm2um', 'I_ergscm2um']
    # Conversion erg/s/cm^2/um to mJy
    model_iso['Ic_Jy'] = u_conversion.ergscmum_to_mjy(model_iso['Ic_ergscm2um'], model_iso['lambda_um'])/1000
    # Multiplying intensities by factor
    model_iso['Ic_Jy'] = model_iso['Ic_Jy']*(factor**2.)
    model_iso['Ic_ergscm2um'] = model_iso['Ic_ergscm2um']*(factor**2.)
    
    ## Spire
    model_spire = pd.read_csv(model_path+'/'+model_name+'_.spire', delim_whitespace= True)
    model_spire.columns = ['lambda_um', 'I_Ic_ergscm2um', 'Ic_ergscm2um', 'Inorm_ergscm2um', 'I_ergscm2um']
    # Conversion erg/s/cm^2/um to mJy
    model_spire['Ic_Jy'] = u_conversion.ergscmum_to_mjy(model_spire['Ic_ergscm2um'], model_spire['lambda_um'])/1000
    # Multiplying intensities by factor
    model_spire['Ic_Jy'] = model_spire['Ic_Jy']*(factor**2.)
    model_spire['Ic_ergscm2um'] = model_spire['Ic_ergscm2um']*(factor**2.)
    
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
    
    
    orig_luminosity_factor = (luminosity_total/(factor**2.))/luminosity_obs
    new_factor = np.sqrt(1./ orig_luminosity_factor)
    new_luminosity_total =(luminosity_total/(factor**2.))*(new_factor**2)# lsun
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
        fig = plt.figure()
        gs1 = gridspec.GridSpec(m, n)    
        gs1.update(wspace = 0.0, hspace=0.0, top=0.95, bottom = 0.05)
        # Generating specified number of axis
        axis = []
        for i in range(m*n):
            axis.append(fig.add_subplot(gs1[i]))
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
        
        fig.savefig(figoutdir+'/'+model_name+'_f'+'%1.1f' % factor +'_SED.'+fig_format, bbox_inches='tight', transparent=True, dpi=400, format=fig_format)
        plt.close()
            
    
    
    return luminosity_total, peak_v0_39_38, peak_v0_24_23, peak_v7_39_38, peak_v7_24_23
