#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 10:06:11 2017

@author: frico
"""

#==============================================================================
#                   For model parameters in Radiative Transfer Code
#==============================================================================


#==============================================================================
# Packages
#==============================================================================

version = '0'

import numpy as np
import os
import pandas as pd
import glob
from copy import deepcopy
from ferpy import u_conversion
from ferpy import fir_lines
import gc

def pars_wrapper(model_paths, model_names, run_ref_only, 
                 luminosity_total, peak_v0_39_38, peak_v0_24_23, peak_v7_39_38,
                 peak_v7_24_23, ratio_24_23, ratio_39_38, ratio_v7, ratio_v0,
                 factor, rad_profile_list,
                  imp, suffix):
    """
    Function to get model parameters from input model .inp files
    Need to give path to .inp models and 
    if we want to save it (imp=True) or just print it (imp=False)
    only_fort55=True only get pars from models in mapas_greg.txt
    """
    
    # Starting dir
    origdir =  os.getcwd()
    # Change to Data dir
    os.chdir(model_paths)
    
    # Constants
    k100 = 44.5             # Opacidad a 100 um cm^2 /g 
    mh = 1.67352*10**-24    # M hidrogeno g
    
    # Model list
    if run_ref_only==False:
        model_list = []
        for i,mod in enumerate(model_names):
            model_list.append(mod+'.inp')
    else:
        model_list = glob.glob(run_ref_only[1]+'.inp')*len(factor)
    
    # Model input parameters
    a_dust = []         # Abundancia polvo 
    t_dust = []         # Temperatura polvo K
    a_molec1   = []     # Abundancia molecula 1
    nh = []             # Densidad Hidrogeno cm-3
    rnube = []          # Radio de la nube cm
    vturb = []          # Velocidad turbulencia km/s8
    tgas = []           # Temperatura del gas K
    model_name = []          # Model name list
    distance = []          # Model distance in pc
    central_T = []
    central_nh = []
    
    # Col_dens lists
    Ndust = []
    Nmolec1 = []
    NH = []
    for i, modelo in enumerate(model_list):
        a_dust_mod = []         # Abundancia polvo 
        t_dust_mod = []         # Temperatura polvo K
        a_molec1_mod   = []     # Abundancia molecula 1
        nh_mod = []             # Densidad Hidrogeno cm-3
        rnube_mod = []          # Radio de la nube cm
        vturb_mod = []          # Velocidad turbulencia km/s8
        tgas_mod = []           # Temperatura del gas K
        model_name_mod = []          # Model name list
        distance_mod = []          # Model distance in pc
        
        tau_v0_24_23 = []       # Opacities for each line for
        tau_v0_39_38 = []       # the closest impact parameter
        tau_v7_24_23 = []
        tau_v7_39_38 = []
        
        
        # Reading input parameters from model.inp
        with open(model_paths+'/'+modelo, 'r') as file:
                        # read a list of lines into data
                        model = file.readlines()
        file.close()            
        rnube_tot = float(model[11].split()[1])
        for j in range(33,43):
            rc_m    = float(model[j].split()[2])
            nh_m    = float(model[j].split()[3])
            a_m     = float(model[j].split()[4])
            vkin    = float(model[j].split()[5])
            tgas_m  = float(model[j].split()[6])
            
            rnube_mod.append(rc_m)
            nh_mod.append(nh_m)
            a_molec1_mod.append(a_m)
            vturb_mod.append(vkin)
            tgas_mod.append(tgas_m)
        for j in range(53,63):
            a_dust_mod.append(float(model[j].split()[1]))
            t_dust_mod.append(float(model[j].split()[2]))
            
        h = 0
        Nd_tot  = 0
        Nmol_tot= 0
        NH_tot  = 0
        print '\nCalculating Column densities:'
        # Calculating Column densities per layer
        for r, rad in enumerate(rnube_mod):
            if h < (len(rnube_mod)-1):
                Nd_tot_i = (rnube_mod[h+1]-rnube_mod[h])*a_dust_mod[h]*nh_mod[h]
                Nd_tot = Nd_tot + Nd_tot_i
                
                Nmol_tot_i =( rnube_mod[h+1]-rnube_mod[h])*a_molec1_mod[h]*nh_mod[h]
                Nmol_tot = Nmol_tot + Nmol_tot_i
                
                NH_tot_i =  (rnube_mod[h+1]-rnube_mod[h])*nh_mod[h]
                NH_tot = NH_tot + NH_tot_i
            else:
                Nd_tot_i = (rnube_tot-rnube_mod[h])*a_dust_mod[h]*nh_mod[h]
                Nd_tot = Nd_tot + Nd_tot_i
                
                Nmol_tot_i = (rnube_tot-rnube_mod[h])*a_molec1_mod[h]*nh_mod[h]
                Nmol_tot = Nmol_tot + Nmol_tot_i
                
                NH_tot_i =  (rnube_tot-rnube_mod[h])*nh_mod[h]
                NH_tot = Nmol_tot + NH_tot_i 
            print '\t'+str(h+1)+':\tNd = '+'%.3E' % Nd_tot_i + '\tNmol = '+'%.3E' % Nmol_tot_i
            h = h+1
            
            
            
        print '\tTotal:\tNd = '+'%.3E' % Nd_tot + '\tNmol = '+'%.3E' % Nmol_tot
            
        model_name_mod.append(modelo.split('.inp')[0]) # To delete the .inp
        distance_mod.append(float(model[27].split()[0]))  
        
        rnube.append(float(model[11].split()[1]))
        nh.append(np.mean(nh_mod))
        a_molec1.append(np.mean(a_molec1_mod))
        vturb.append(float(model[40].split()[5]))
        tgas.append(np.mean(tgas_m))
        central_T.append(tgas_mod[0])
        central_nh.append(nh_mod[0])
        #tgas.append(tgas_mod[0])
        a_dust.append(np.mean(a_dust_mod))
        t_dust.append(np.mean(t_dust_mod[0]))
        model_name.append(modelo.split('.inp')[0]) # To delete the .inp
        distance.append(float(model[27].split()[0]))
        
        Ndust.append(Nd_tot)
        Nmolec1.append(Nmol_tot)
        NH.append(NH_tot)
        
        
        # Reading opacities
        #modelo = 'mhc3n_r0.0_n5.00E+06_T300.inp'
        #model_paths = '/Users/frico/Documents/Ed/modelos_hc3n/v8_tcte2/r_0.0_X1.0E-07'
        tauvel = modelo.split('.inp')[0] + '_1rec.tauvel'
        model_tau = pd.read_csv(model_paths+'/'+tauvel, delim_whitespace= True, header=None)
        model_tau.columns = ['NPar', 'ParImp/Rmax', 'Vel', 'tau_v0_24_23', 'tau_v0_39_38', 'tau_v7_24_23', 'tau_v7_39_38', 'tau_Line5', 'tau_Line6']
        
        # Selecing opacities at the peak of the line and for the closest impact parameter to the center
        parimp = model_tau['ParImp/Rmax'].min()
        tau_sel = model_tau[(model_tau['ParImp/Rmax'] == parimp) & (model_tau.Vel == 0.0)]
        
        tau_v0_24_23.append(tau_sel['tau_v0_24_23'].tolist()[0])
        tau_v0_39_38.append(tau_sel['tau_v0_39_38'].tolist()[0])
        tau_v7_24_23.append(tau_sel['tau_v7_24_23'].tolist()[0])
        tau_v7_39_38.append(tau_sel['tau_v7_39_38'].tolist()[0])
        
        # Abundance refered to total column density
        #a_molec1.append(Nmol_tot/np.mean(nh_mod)*float(model[11].split()[1]))
  
    # tau100
    tau100 = np.array(Ndust) * 2. * mh * k100
    
    #luminosity_total, peak_v0_39_38, peak_v0_24_23, peak_v7_39_38, peak_v7_24_23, 
    # Creating df
    dictio = {'NH': NH,
              'nH': nh,
              'rnube': rnube,
              'vturb': vturb,
              'Tgas': tgas,
              'tau100': tau100,
              #'tau100_all': tau100_all,
              #'tau100_tot': tau100_tot,
              'Tdust':  t_dust,
              'Ndust':  Ndust,
              'Xdust':  a_dust,
              'Nline':  Nmolec1,
              'Xline':  a_molec1,
              #'Ltot' : luminosity_total,
              #'peak_v0_24_23' : peak_v0_24_23,
              #'peak_v0_39_38' : peak_v0_39_38,
              #'peak_v7_24_23' : peak_v7_24_23,
              #'peak_v7_39_38' : peak_v7_39_38,
              #'ratio_24_23' : ratio_24_23,
              #'ratio_39_38' : ratio_39_38,
              #'ratio_v7' : ratio_v7,
              #'ratio_v0' : ratio_v0,
              'tau_v0_24_23': tau_v0_24_23*len(model_name), 
              'tau_v0_39_38': tau_v0_39_38*len(model_name),
              'tau_v7_24_23': tau_v7_24_23*len(model_name),
              'tau_v7_39_38': tau_v7_39_38*len(model_name),
              'central_nh'  : central_nh,
              'central_T'   : central_T,
              'model': model_name
                }
    N_df = pd.DataFrame.from_dict(dictio)
    # Ordering df columns
    N_df = N_df[['Tdust', 'NH', 'nH', 'rnube', 'vturb', 'Tgas', 'tau100', 'Ndust', 'Xdust', 'Nline', 'Xline',
                 'tau_v0_24_23', 'tau_v0_39_38', 'tau_v7_24_23', 'tau_v7_39_38', 'central_nh', 'central_T',
                 'model']]
    df_save = N_df.sort_values('Tdust')
    
    # Juntarlo bien por nombres de modelos de ambos df!
    df_save['Ltot'] = luminosity_total
    df_save['peak_v0_24_23'] = peak_v0_24_23
    df_save['peak_v0_39_38'] = peak_v0_39_38
    df_save['peak_v7_24_23'] = peak_v7_24_23
    df_save['peak_v7_39_38'] = peak_v7_39_38
    df_save['ratio_24_23'] = ratio_24_23
    df_save['ratio_39_38'] = ratio_39_38
    df_save['ratio_v7'] = ratio_v7
    df_save['ratio_v0'] = ratio_v0
    df_save['factor'] = factor
    #df_save['central_nh'] = central_dens_list
    #df_save['central_T'] = central_temp_list
    df_save['profile'] = rad_profile_list

    # Assuming same distance for all models
    df_save['ss_arcsec'] = u_conversion.ang_size(distance[0]/1e6, df_save['factor']*2.*rnube/100.) # Multiplicar por 2 (es el diametro el source size)

    if imp:
        # Saving df
        #decimales = pd.Series([0, 3, 0, 3,0,0,3,3,3,3,2,4,2,2,2,2,2,2,2,2,2,2], index=['Tdust', 'NH', 'nH', 'rnube', 'vturb',
        #                      'Tgas', 'tau100', 'Ndust', 'Xdust', 'Nline', 'Xline',
        #         'Ltot', 'ss_arcsec','peak_v0_24_23', 'peak_v0_39_38', 'peak_v7_24_23', 'peak_v7_39_38', 
        #         'ratio_24_23', 'ratio_39_38', 'ratio_v7', 'ratio_v0', 'factor'
        #         ])
        #rounded = df_save.round(decimales)
        rounded = df_save[['Tdust', 'NH', 'nH', 'rnube', 'vturb', 'Tgas', 'tau100', 'Ndust', 'Xdust', 'Nline', 'Xline',
                 'Ltot', 'ss_arcsec','peak_v0_24_23', 'peak_v0_39_38', 'peak_v7_24_23', 'peak_v7_39_38', 
                 'tau_v0_24_23', 'tau_v0_39_38', 'tau_v7_24_23', 'tau_v7_39_38',
                 'ratio_24_23',  'ratio_39_38', 'ratio_v7', 'ratio_v0', 'factor', 'central_nh', 'central_T', 'profile' ,'model'
                 ]]
        if run_ref_only==False:
            rounded.to_csv(model_paths+'/parameter_model_inp_res_'+suffix+'.txt',
                                               sep='\t', index=False, na_rep=-1)
        else:
            rounded.to_csv(model_paths+'/parameter_model_inp_res_'+suffix+'.txt',
                                               sep='\t', index=False, na_rep=-1)
    else:
        # Printing
        print df_save
    # Returning to original dir
    os.chdir(origdir)
    gc.collect()
    return rounded

