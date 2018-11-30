#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 16:16:07 2017

@author: frico
"""
#==============================================================================
#  Python script to run radiative transfer models
#==============================================================================

import os
from copy import deepcopy
import string
from ferpy import radtransf_modelpars_hc3n_v6 as radtransf_modelpars
from ferpy import radtransf_model_runner_hc3n_v6 as radtransf_runner
from ferpy import radtransf_HC3N_plot_and_luminosity_v6 as hc3n_plotter
import gc

def model_conv_sing(results_path, reference_model, run_ref_only,
                    nH_dens, taudust, abun_line,
                    tgas, tdust, vturb, factor):
    
    """
    Function to modify parameters from existing CII models
    Modifies: Tgas, Abun2, n(H2)
    abun_line_max = 0 --> not modifying abundances
    compmodel_only = True --> only do composite models
    """
    if run_ref_only==False:
        abun_inp = deepcopy(abun_line)
        tdust_inp = deepcopy(tdust)
        nH_dens_inp = deepcopy(nH_dens)
        tgas_inp = deepcopy(tgas)
        vturb_inp = deepcopy(vturb)
        taudust_inp = deepcopy(taudust)
        
        
        # Path to reference model
        #reference_model = '/Users/frico/Desktop/ferpy/ferpy/rad_transfer/CII_models'
        # Input model directory
        #inputdir = '/Users/frico/Documents/Ed/modelos_hc3n/modelos2'
        
        
            
        # ISM Conditions
        # High density  HII         nH = 10, 100, 1000   Tgas = 1E4       (Langer 2015, Goldsmith 2015)
        # Low density   HII (WIM)   nH = 0.01 - 0.1      Tgas = 100 - 1E4 (Langer 2015, Goldsmith 2015)
        # PDR density               nH = 1000            Tgas = 150       (Malhorta et al 2001)
        
        if not os.path.exists(results_path):
            os.makedirs(results_path)
        os.chdir(results_path)
            
        # Change to Data dir
        os.chdir(results_path)
        # List to save model parameters
    
        # Reading input parameters from model.inp
        modelo = reference_model
        mod_name = modelo.split('.inp')[0] # To delete the .inp
        with open(results_path+'/'+modelo, 'r') as file:
                        # read a list of lines into data
                        model = file.readlines()
    
        print mod_name
        print modelo
        
        # Gas Properties
        for m in range(33, 43):
            gas_propert = model[m].split()
            # Original pars
            nH_dens = float(gas_propert[3])
            tgas = float(gas_propert[6])
            vturb = float(gas_propert[5])
            a_molec1 = float(gas_propert[4])
            
        # Dust Properties
        for dd in range(53,63):
            dust_prop = model[dd].split()
            abundust = float(dust_prop[1])
            Tdust = float(dust_prop[2])
            
        # Original model properties (if they are constant in the cloud)
        
        # Optical Depth at 100 um
        k100 = 44.5             # Opacidad a 100 um (cm^2 /g) 
        mh = 1.67352*10**-24    # M hidrogeno (g)
        # R nube
        rnube = float(model[11].split()[1])
        # H density
        nhden = nH_dens
        # NH
        NH = nhden * rnube
        # Abund Dust
        a_dust = abundust
        # Ndust
        Ndust = a_dust * nhden * rnube
        # tau100
        tau100 = Ndust * 2. * mh * k100
        # Abun OH
        a_oh = a_molec1
        # Noh
        Noh = nhden * rnube * a_oh
        # Abund line
        a_line = a_molec1
        # Nline
        Nline= nhden * rnube * a_line
        # NH from tau100
        NH_tau100 = tau100/(a_dust*2.*mh*k100)
        nh_tau100 = NH_tau100/rnube
            
        # New properties
        #if taudust_inp != 'cte':
        # tau100_change = nH_change * a_dust_change
        # tau100 is input, nH is input, we just dont know the change in the abundances
        if nH_dens_inp != 0:
            tau_change  = taudust_inp/tau100
            nh_change   = nH_dens_inp/nh_tau100
            a_change    = tau_change/nh_change
            
            # New abundances:
            a_dust_new = a_change * a_dust
            a_oh_new   = a_change * a_oh
            a_line_new = a_change * a_line
            
                
            
            if abun_inp != 0:
                a_line_new = abun_inp
            # Oh
            Nmolec1_new = nH_dens_inp * rnube * a_oh_new
            # Nline
            Nmolec2_new = nH_dens_inp * rnube * a_line_new
            # Ndust
            Ndust_new   = nH_dens_inp * rnube * a_dust_new
            # tau100
            tau100_new = Ndust_new * 2. * mh * k100
            
        nh_new      = deepcopy(nH_dens_inp)
        vturb_new   = deepcopy(vturb_inp)
        tgas_new    = deepcopy(tgas_inp)
        tdust_new   = deepcopy(tdust_inp)
                
        # Changing Gas Properties
        for m in range(33, 43):
            no_oh_abun = model[m].split()
            #               N rayos            tipo capa           radio interior      dens H2 (cm-3)           Abun OH              Abun Line                   Vturb                 Tgas     
            #model[m] = '  '+no_oh_abun[0]+'   '+no_oh_abun[1]+'    '+no_oh_abun[2]+'   '+nH_dens+'   '+no_oh_abun[4]+'  '+str(float(no_oh_abun[5])*abun_factor)+'   '+no_oh_abun[6]+'  '+tgas+'\n'
            model[m] = '  '+no_oh_abun[0]+'   '+no_oh_abun[1]+'    '+no_oh_abun[2]+'   '+"%.5E" % nH_dens+'   '+"%.3E" % a_line+'   '+"%.1f" % vturb_new+'  '+"%.2E" % tgas_new+'\n'
        # Changing Dust properties
        for dd in range(53,63):
            dust_prop = model[dd].split()
            model[dd] = ' '+dust_prop[0]+' '+"%.4E" % a_dust+' '+"%.1f" % tdust_new+'\n'
        
        # Model name
        mod_name='mhc3n_2_tgas'+"%1.f" % tgas_new+'_tdust' +'%1.f' % tdust_new
        model[32] = "'"+results_path+'/'+mod_name+"_'"+'\n'
        # Writting Model input
        model_name = mod_name
        with open(results_path+'/'+model_name+'.inp', 'w') as file:
            	file.writelines(model)
        file.close()
        model_name = mod_name
        print '\t'+model_name
        radtransf_runner.model_runner(results_path+'/', model_name)
    else:
        # Running only given model
        radtransf_runner.model_runner(results_path+'/', reference_model)
        mod_name = reference_model
    
        
    
    #luminosity_total, peak_v0_39_38, peak_v0_24_23, peak_v7_39_38, peak_v7_24_23 = hc3n_plotter.luminosity_calculator(model_path=results_path,
    #                                                                                                     model_name=mod_name,
    #                                                                                                     factor=factor,
    #                                                                                                     plot_profiles=True, plot_SED=True,
    #                                                                                                     figoutdir=results_path+'/figures', fig_format='png')

    
    # Function to make a summary of the model parameters
    #df_pars = radtransf_modelpars.pars_wrapper(results_path, luminosity_total, peak_v0_39_38, peak_v0_24_23, peak_v7_39_38, peak_v7_24_23, True)
    
    gc.collect()
    return mod_name
    
def model_builder(results_path, model_name,
                    nH_dens, abun_line, abun_dust,
                    tgas, tdust, vturb):
    
    """
    Function to modify parameters from existing reference models
    """
    
    reference_model = '/Users/frico/Desktop/ferpy/molecules/HC3N/HC3N_model_builder.inp'
    
    
    # ISM Conditions
    # High density  HII         nH = 10, 100, 1000   Tgas = 1E4       (Langer 2015, Goldsmith 2015)
    # Low density   HII (WIM)   nH = 0.01 - 0.1      Tgas = 100 - 1E4 (Langer 2015, Goldsmith 2015)
        
    if not os.path.exists(results_path):
        os.makedirs(results_path)
    # Change to Data dir
    os.chdir(results_path)
            
    # Reading input parameters from model.inp
    modelo = reference_model
    mod_name = modelo.split('.inp')[0] # To delete the .inp
    with open(modelo, 'r') as file:
        # read a list of lines into data
        model = file.readlines()
        
            
    # Changing Gas Properties
    index = 0
    for m in range(33, 43):
        model[m] = string.replace(model[m]  , 'dens', "%.5E" % nH_dens[index])
        model[m] = string.replace(model[m]  , 'a_molec', "%.3E" % abun_line[index])
        model[m] = string.replace(model[m]  , 'vturb', "%.1f" % vturb[index])
        model[m] = string.replace(model[m]  , 'tgas', "%.3E" % tgas[index])
        index = index +1 
        #no_oh_abun = model[m].split()
        #               N rayos            tipo capa           radio interior      dens H2 (cm-3)           Abun OH              Abun Line                   Vturb                 Tgas     
        #model[m] = '  '+no_oh_abun[0]+'   '+no_oh_abun[1]+'    '+no_oh_abun[2]+'   '+nH_dens+'   '+no_oh_abun[4]+'  '+str(float(no_oh_abun[5])*abun_factor)+'   '+no_oh_abun[6]+'  '+tgas+'\n'
        #model[m] = '  '+no_oh_abun[0]+'   '+no_oh_abun[1]+'    '+no_oh_abun[2]+'   '+"%.5E" % nH_dens+'   '+"%.3E" % a_line+'   '+"%.1f" % vturb_new+'  '+"%.2E" % tgas_new+'\n'
    index2 = 0
    # Changing Dust properties
    for dd in range(53,63):
        model[dd] = string.replace(model[dd]  , 'a_dust', "%.3E" % abun_dust[index2])
        model[dd] = string.replace(model[dd]  , 'tdust', "%.3E" % tdust[index2])
        index2 = index2 +1 
        #dust_prop = model[dd].split()
        #model[dd] = ' '+dust_prop[0]+' '+"%.4E" % a_dust+' '+"%.1f" % tdust_new+'\n'
    
    # Out Model name
    out_mod_name = model_name.split('.inp')[0] 
    model[32] = "'"+results_path+'/'+out_mod_name+"_'"+'\n'
    file.close()
    # Writting New Model
    with open(results_path+'/'+out_mod_name+'.inp', 'w') as file:
        	file.writelines(model)
    file.close()

    print '\t'+model_name
    radtransf_runner.model_runner(results_path+'/', out_mod_name)
   
    return out_mod_name
