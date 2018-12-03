#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  3 13:27:06 2018

@author: frico
python 3.7.1
"""
version = '2'

"""
He comentado las lineas 79-82 de colis/hc3n_h2_faure16.f
"""


from ferpy import radtrasf_hc3n_runner_v6 as radtransf_hc3n
from ferpy import radtransf_modelpars_hc3n_v6 as radtransf_modelpars
from ferpy import radtransf_HC3N_plot_and_luminosity_v6 as hc3n_plotter # old one
#from ferpy import radtransf_HC3N_plot_and_luminosity_v6_noabs as hc3n_plotter # To avoid selecting ratios with abs
from ferpy import u_conversion
from ferpy import utiles
from ferpy import radtransf_model_calculations as m_calc

import scipy
import numpy as np
import os
import pandas as pd
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy import wcs
from radio_beam import Beam
from glob import glob
import resource # to control memory usage
import gc


class data:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        
modelling = True

def models_run(start_dens, profile, ref_coldust, ref_colmol, luminosity,
               gradient, Tdust_mod, vturbulence, factor_names, factors_array,
               run_ref_only, plot_profiles, plot_SED,
               model_name_list, results_path, build_from_abunds, only_pars=False
               ):
    '''
    ref_coldust/ref_abdust if build_from_abunds=False/True
    ref_colmol/ref_abmol if build_from_abunds=False/True
    '''
    
    ########## Observed values HC14 ##########
    luminosity_total_obs = 3E+09
    peak_v0_24_23_obs = 54.08
    peak_v0_39_38_obs = 70.29
    peak_v7_24_23_obs = 22.05
    peak_v7_39_38_obs = 44.59
    size_obs = 0.01995
    ratio_2423_obs = peak_v0_24_23_obs/peak_v7_24_23_obs
    ratio_3938_obs = peak_v0_39_38_obs/peak_v7_39_38_obs
    ratio_v0_obs = peak_v0_39_38_obs/peak_v0_24_23_obs
    ratio_v7_obs = peak_v7_39_38_obs/peak_v7_24_23_obs
    
    ########## Radii ##########
    radio_list = np.array([3E17, 6E17,9E17, 1.2E18, 1.5E18, 1.8E18, 2.1E18, 2.4E18, 2.7E18, 2.8E18])

    # Parameters output
    suffix = 'last_model'
    
    if only_pars == True:
        model_name_list = []
        model_list_all = glob(results_path+'/*.inp')
        for m, model in enumerate(model_list_all):
            model_name_list.append(model.split('/')[-1].split('.inp')[0])
    else:
        # Checking if name already exists
        for f, fn in enumerate(factor_names):
            if os.path.isfile(results_path+'/models_summary_'+fn+'.txt'):
                mod_short = model_name_list[0].split('.inp')[0]
                df_sum = pd.read_csv(results_path+'/models_summary_'+fn+'.txt', delim_whitespace= True, header=0, comment='#')
                if df_sum['model'].str.contains(mod_short).any():
                    modn1 = mod_short.split('_')
                    if modn1[-1].isdigit():
                        modn1[-1] = str(len(df_sum['model'])+1)
                        model_name_list[0] = '_'.join(modn1)
                    # Input other name if already exists
                    #print '\n\tModel Name already used, change it:'
                    #model_name_list[0] = raw_input('Enter new name:\t')
                    while df_sum['model'].str.contains(model_name_list[0].split('.inp')[0]).any():
                        print('\n\tModel Name already used, change it:')
                        model_name_list[0] = raw_input('Enter new name:\t')
            print('\t\t'+model_name_list[0])
        # Transforming model parameters
        #tlist = [100,150,200,250,300,350,400,450,500,550,600]
        #factor = [0.7]*len(tlist)#[0.7,0.35]
        #tlist = [1,1]
    
    
        
    
        ########## Density profile ##########
        dens_profile, dens_mean, mass_tot= m_calc.densprofile(start_dens, radio_list, profile)
        #mass_tot, dens_mean = m_calc.Mass_from_dens(dens_profile, radio_list)
    
        ########## Generating abunds ##########
        if build_from_abunds == False:
            abunds_densgrad = m_calc.abunds_capa([ref_coldust/10]*len(radio_list), [ref_colmol/10]*len(radio_list), dens_profile, radio_list)
            abun_dust = abunds_densgrad[0]
            abun_line = abunds_densgrad[1]
        else:
            abun_dust = [ref_coldust]*len(radio_list)
            abun_line = [ref_colmol]*len(radio_list)
            
        ########## Making all Col density the same ##########
        Ncols_tots, Ncols_lists = m_calc.Col_dens(abun_dust, abun_line, dens_profile, radio_list)
    
        ########## Generating Temperature Gradient ################
        tdust_gradient = []
        print('\nModel Parameters:')
        for i, rad in enumerate(radio_list):
            if gradient == 'manual' or gradient == 'manual_dens':
                tdust_gradient.append(np.round(Tdust_mod[i],3))
                ttt_st = tdust_gradient[0]
            else:
                td = m_calc.tdust_devicente(luminosity, rad, gradient, Tdust_mod)
                tdust_gradient.append(np.round(td,3))
                #print td
                #print tdust_gradient
                ttt_st = tdust_gradient[i]
            
            #print dens_profile[i]
            #print abun_dust[i]
            #print Ncols_lists[0][i]
            #print abun_line[i]
            #print Ncols_lists[1][i]
            #print tdust_gradient[i]
            
            print(str(i+1)+':\t nH = '+ '%1.2E' % dens_profile[i]+'\tXd = '+ '%1.2E ' % abun_dust[i] +'\tNd = '+ '%1.2E ' % Ncols_lists[0][i] +'\t Xm = '+ '%1.2E ' % abun_line[i]+'\tNm = '+ '%1.2E ' % Ncols_lists[1][i] +'\t T= '+'%1.2f ' % ttt_st)
        print('Mean\tnH = '+'%1.2E' % dens_mean + ' \tRef Mean\tnH = '+'%1.2E' % dens_mean_or)
        print('\tNd_tot = '+'%1.2E' % Ncols_tots[0] + ' \tNdmtot = '+'%1.2E' % Ncols_tots[1])
    
        ########## Turbulence Velocity ################
        vturb = [vturbulence]*len(radio_list)


    luminosity_total_list = []
    peak_v0_39_38_list = []
    peak_v0_24_23_list = []
    peak_v7_39_38_list = []
    peak_v7_24_23_list = []
    ratio_24_23_list = []
    ratio_39_38_list = []
    ratio_v7_list = []
    ratio_v0_list = []
    model_names = []
    factor_list = []
    central_dens_list = []
    rad_profile_list = []
    central_temp_list = []
    print('\nRunning model:')
    for i, m_name in enumerate(model_name_list):
        if only_pars == False:
            mod_name = radtransf_hc3n.model_builder(results_path=results_path, model_name=m_name,
                        nH_dens=dens_profile, abun_line=abun_line, abun_dust=abun_dust,
                        tgas=tdust_gradient, tdust=tdust_gradient, vturb=vturb)
        else:
            mod_name = m_name.split('.inp')[0] 
        print('\t'+mod_name)
        
        
        #mod_name = radtransf_hc3n.model_conv_sing(results_path=results_path, reference_model=reference_model, run_ref_only=run_ref_only,
        #                           nH_dens=0, taudust=0, abun_line=0,
        #                           tgas=temp, tdust=temp, vturb=20, factor=factor[i])
        for f, factor in enumerate(factors_array):
            luminosity_total, peak_v0_39_38, peak_v0_24_23, peak_v7_39_38, peak_v7_24_23 = hc3n_plotter.luminosity_calculator(model_path=results_path,
                                                                                                             model_name=mod_name,
                                                                                                             factor=factor,
                                                                                                             plot_profiles=plot_profiles, plot_SED=plot_SED,
                                                                                                             figoutdir=results_path+'/figures', fig_format='png')
        
            
            print('\n Comparing values HC14 - Factor'+factor_names[f]+':')
            print('Modelled:\tL='+ '%1.2E' % luminosity_total+'\tv0_24-23='+ '%1.2f ' % peak_v0_24_23 +'\tv0_39-38='+ '%1.2f ' % peak_v0_39_38+'\tv7_24-23='+'%1.2f ' % peak_v7_24_23+'\tv7_39-38='+'%1.2f' % peak_v7_39_38)
            print('Observed:\tL='+ '%1.2E' % luminosity_total_obs+'\tv0_24-23='+ '%1.2f ' % peak_v0_24_23_obs +'\tv0_39-38='+ '%1.2f ' % peak_v0_39_38_obs+'\tv7_24-23='+'%1.2f ' % peak_v7_24_23_obs+'\tv7_39-38='+'%1.2f ' % peak_v7_39_38_obs)
            print('\n')
            ratio_2423 = peak_v0_24_23/peak_v7_24_23
            ratio_3938 = peak_v0_39_38/peak_v7_39_38
            ratio_v0 = peak_v0_39_38/peak_v0_24_23
            ratio_v7 = peak_v7_39_38/peak_v7_24_23
            print('Modelled:\tR_24-23='+ '%1.2f' % ratio_2423+'\tR_39-38='+ '%1.2f ' % ratio_3938 +'\tR_v0='+ '%1.2f ' % ratio_v0+'\tR_v7='+'%1.2f ' % ratio_v7)
            print('Observed:\tR_24-23='+ '%1.2f' % ratio_2423_obs+'\tR_39-38='+ '%1.2f ' % ratio_3938_obs +'\tR_v0='+ '%1.2f ' % ratio_v0_obs+'\tR_v7='+'%1.2f ' % ratio_v7_obs)
        
            model_names.append(mod_name)
            luminosity_total_list.append(luminosity_total)
            peak_v0_39_38_list.append(peak_v0_39_38)
            peak_v0_24_23_list.append(peak_v0_24_23)
            peak_v7_39_38_list.append(peak_v7_39_38)
            peak_v7_24_23_list.append(peak_v7_24_23)
            ratio_24_23_list.append(peak_v0_24_23/peak_v7_24_23)
            ratio_39_38_list.append(peak_v0_39_38/peak_v7_39_38)
            ratio_v7_list.append(peak_v7_39_38/peak_v7_24_23)
            ratio_v0_list.append(peak_v0_39_38/peak_v0_24_23)
            factor_list.append(factor)
            #central_dens_list.append(start_dens)
            rad_profile_list.append(profile)
            #central_temp_list.append(tdust_gradient[0])
    # Function to make a summary of the model parameters
    df_pars = radtransf_modelpars.pars_wrapper(results_path, model_names, run_ref_only, luminosity_total_list, peak_v0_39_38_list,
                                               peak_v0_24_23_list, peak_v7_39_38_list, peak_v7_24_23_list,
                                               ratio_24_23_list, ratio_39_38_list, ratio_v7_list, ratio_v0_list,
                                               factor_list, rad_profile_list,
                                               True, suffix)
    
    for f, factor in enumerate(factors_array):
        if f == 0:
            ind3 = 1
        elif f == 1:
            ind3 = 2
        df_factor = df_pars[df_pars['factor'] == factor]
        # Appending to models result with factor 0.7
        if only_pars == True:
            if not os.path.isfile(results_path+'/models_summary_'+factor_names[f]+'.txt'):
                df_factor.to_csv(results_path+'/models_summary_'+factor_names[f]+'.txt', header=True, sep='\t', index=False, na_rep=-1)
            else:
                df_factor.to_csv(results_path+'/models_summary_'+factor_names[f]+'.txt', mode='a', header=False, sep='\t', index=False, na_rep=-1)
        else:
            if not os.path.isfile(results_path+'/models_summary_'+factor_names[f]+'.txt'):
                df_factor.to_csv(results_path+'/models_summary_'+factor_names[f]+'.txt', header=True, sep='\t', index=False, na_rep=-1)
            else:
                df_factor.to_csv(results_path+'/models_summary_'+factor_names[f]+'.txt', mode='a', header=False, sep='\t', index=False, na_rep=-1)
        
    gc.collect()
    return df_pars


        



## Reference Model
#reference_model = 'mhc3n_2_densandtgrad_12' #without .inp

## Results path
#results_path = '/Users/frico/Documents/Ed/modelos_hc3n/v0/n_cte_t150'
#run_ref_only_path = results_path+'/'+reference_model
## Model name (if building from 0)
#model_name_list = ['mhc3n_ncte_t350_1.inp'] # With the .inp

########## Radii ##########
radio_list = np.array([3E17, 6E17, 9E17, 1.2E18, 1.5E18, 1.8E18, 2.1E18, 2.4E18, 2.7E18, 2.8E18])


########## Original pars from first model (Ed) ##########
ref_coldust_or = 3.086E22
ref_colmol_or  = 1.0E17
nh_or = 3.2e6
original_dens_profile = [3.2e6]*len(radio_list)
mass_tot_or, dens_mean_or = m_calc.Mass_from_dens(original_dens_profile, radio_list)


"""
model_pars
ref_colmol/abunds
ref_coldust/abunds

"""    
    
# =============================================================================
#  H2 Density Parameters
# =============================================================================

# Density of the smallest shell
start_dens_list =   [2e7, 4e7,7e7, 8e7] 
# Density profile exponent
profile_list    =   [0.0]    

# =============================================================================
#  Temperatures
# =============================================================================
"""
 Luminosity for HC if size = 0.1
 """
luminosity_list  = [5e8, 7e8, 9e8, 1e9, 1.5e9, 2.0, 2.9e9]    
gradient = 'manual_dens' #'tgrad' #'manual' #False # Gradient temperature
Tdust_mod = 150  # Temperature gradient if Tdust_mod = False or manual
Tdust_mod_list = [[300]*len(radio_list)]
"""
Tlist = [[100]*len(radio_list), [150]*len(radio_list),
 [200]*len(radio_list), [250]*len(radio_list),
 [300]*len(radio_list), [350]*len(radio_list),
 [400]*len(radio_list), [450]*len(radio_list),
 [500]*len(radio_list), [550]*len(radio_list),
 [600]*len(radio_list)]
"""
# =============================================================================
#  Abundances
# =============================================================================
# Desired Dust abundances:
abund_dust_list = [1E-2]

# Desired HC3N abundances:
low_abun = [1E-7,1E-8, 1E-9, 1E-10]
high_abun = [5E-6, 5E-7, 1E-7, 5E-8]
abund_hc3n_list = [5E-6, 5E-7, 1E-7, 5E-8, 1E-8, 1E-9, 1E-10]


# =============================================================================
#  Cloud Size
# =============================================================================
# Cloud total radius (cm)
rnube = 3.086E18


# =============================================================================
#  Column densities
# =============================================================================
# Column densities to keep reasonable HC3N and dust abundances
Nhc3n = {}
Ndust= {}
Xhc3n = {}
Xdust = {}
# Proper Column densities to get abundances from
obs_col = [1.00E16, 3.16E16] #[3.16E14, 1.00E15, 3.16E15, 1.00E16, 3.16E16, 1.00E17, 3.16E17, 1E18]#[1E17, 1E16, 1E15, 1E14]
dust_col = [1E22]#[1E19, 5E19, 1E20, 5E20, 1E21, 5E21, 1E22, 5E22, 8E22, 1E23]#[4.0E20, 9E20, 6.172E21, 1E22, 1E23] #[3.086E18, 3.086E19, 3.086E20, 6.172E20, 3.086E21]
for hh, nh in enumerate(start_dens_list):
    Nhc3n['%1.1E' %nh] = {}
    Ndust['%1.1E' %nh] = {}
    Xhc3n['%1.1E' %nh] = {}
    for aaa, abun in enumerate(abund_hc3n_list):
        Nhc3n['%1.1E' %nh]['%1.1E' % abun] = (rnube*abun*nh)
    for aad, abun_dust in enumerate(abund_dust_list):
        Ndust['%1.1E' %nh]['%1.1E' % abun_dust]  = (rnube*abun_dust*nh)   
    for o, col in enumerate(obs_col):
        Xhc3n['%1.1E' %nh]['%1.2E' % col] = col/(nh*rnube)
    


"""
Getting only pars
"""    
get_pars_only = False
"""
Ed model parameters:
    Td=Tgas
    nH2
    Ndust
    Nmolec
    Tamaño
"""
# =============================================================================
#  Starting Model runner
# =============================================================================
############ Model Runner ########################
if modelling == True:
    
    version = '14_a'        # Version to name the models folder
    build_from_abun = True  # Building models from abunds or from col_dens
    vturbulence = 20.0      # Turbulence velocity
    abun_from_col = True    # Abundances from column densities
    
    # Plotting options
    plot_profiles, plot_SED = False, False
    # Runing only reference model
    run_ref_only = False #[True,reference_model] # else: False
    
    factor_names = ['07', '035'] # Factors to get luminosity or line intensities #### MAAL
    factors_array = [0.7, 0.35]
    
    if gradient == True:
        # Look hc3n_modelling previous versions
        print('Nothing to do')
    elif gradient == 'manual_dens':
        """
        Specifying density in each shell
        """
        version = version + '_tcte_T'
        if abun_from_col == True and build_from_abun==True:
            # Calculates the abundances of the molecule and dust given a certain 
            # reference column for the molecule and dust
            decimals = 1
            ref_colmol_list = obs_col
            ref_coldust_list = dust_col
            string = 'X'
        else:
            string = 'N'
        for r, profile in enumerate(profile_list):
            # Loop through density profiles ## Creo que solo vale con 0.0
            for tt, Tdust_mod in enumerate(Tdust_mod_list):
                # Loop through dust temperatures (Tdust = Tkin)
                td = Tdust_mod[0]
                # Writting Results path
                results_path = '/Users/frico/Documents/Ed/modelos_hc3n/v'+version+'/r_'+'%1.1f' % profile +'/r_'+'%1.1f' % profile +'_T'+'%1.f'  % td
                for ro, start_dens in enumerate(start_dens_list):
                    # Loop through densities ### MIRAR si vale con perfil de densidad
                        for cd, ref_coldust in enumerate(ref_coldust_list):
                            
                            # Getting dust abundances given certain column density and density
                            abun_d = ref_coldust/(start_dens*rnube)
                            # rRunding exponential numbers
                            abund_dust = utiles.rounding_exp(abun_d, decimals)
                            
                            # Removing previous model summary if only updating pars summary
                            if os.path.isdir(results_path):
                                if get_pars_only == True:
                                    for f, factor in enumerate(factors_array):
                                        os.rename(results_path+'/models_summary_'+factor_names[f]+'.txt', results_path+'/models_summary_'+factor_names[f]+'_old.txt')
                           
                            
                            for cm, ref_colmol in enumerate(ref_colmol_list):
                                # Loop through desired molecule densities
                                # Getting molecule abundances given certain column density and density
                                abun = ref_colmol/(start_dens*rnube)
                                abund_hc3n = utiles.rounding_exp(abun, decimals)
                                print('\n++++++++++')
                                print(abund_hc3n)
                                print('+++++++++++++')
                                luminosity = 1e9 # Luminosity without Tgrad is not relevant
                                model_name_list = ['mhc3n_r'+'%1.1f' % profile +'_Nd'+'%1.2E' % ref_coldust+'_X'+'%1.2E' % abund_hc3n +'_nH'+'%1.2E' % start_dens+'.inp'] # With the .inp
                                print('\n\tStart Density: '+'%1.3E'  % start_dens)
                                print('\tProfile: '+'%1.3E'  % profile)
                                print('\tNd: '+'%1.3E'  % abund_dust+'\tNHC3N: '+'%1.3E'  % abund_hc3n)
                                Ndddd = rnube*abund_dust*start_dens
                                Nmolll = rnube*abund_hc3n*start_dens
                                print('\tXd: '+'%1.3E'  % abund_dust+'\tXHC3N: '+'%1.3E'  % abund_hc3n)
                                print('\tNd: '+'%1.3E'  % Ndddd+'\tNHC3N: '+'%1.3E'  % Nmolll)
                                print('\tLuminosity: '+'%1.3E'  % luminosity)
                                a = models_run(start_dens, profile, abund_dust, abund_hc3n, luminosity,
                                               gradient, Tdust_mod, vturbulence, factor_names, factors_array,
                                               run_ref_only, plot_profiles, plot_SED,
                                               model_name_list, results_path, build_from_abun, only_pars=get_pars_only
                                               )
                                a = 0
                                gc.collect()
                                print('\tMemory Usage:'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
                                    
                                    
                                    
# =============================================================================
# Observed values MADCUBA
# =============================================================================

# Observed data workdir
dworkdir_cubes = '/Users/frico/Documents/data/NGC253_H3O+'
dworkdir_spec = dworkdir_cubes+'/Hotcores_v4_all'

# Out dir
out_dir = dworkdir_spec+'/Results_v'+version+'/'
if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
# HC3N properties
hc3nv7 = data(lmin='hc3nv7',
              freq=[218.8608*u.GHz, 219.173757*u.GHz, 364.676275*u.GHz, 365.195232*u.GHz],
              lgint=[-1.76, -1.759, -1.419, -1.418],
              elo=[306.961, 307.081, 460.252, 460.590]
              )

# NGC253 properties
D = 3.5 # Mpc to NGC253
vsis = 258 #km/s



# Observation data
hc3n_prop = pd.read_csv(dworkdir_spec+'/tabla_hc3n_v12.txt', delim_whitespace= True, header=0, comment='#')
#hc3n_prop.columns = ['Ind', 'Ind_2', 'RA',	'Dec', 'VLSR', 'VLSR_err', 'v', 'v_err',
#                     'N(HC3N)', 'N(HC3N)_err', 'N(HC3N)_LR','N(HC3N)_LR_err', 'N(rot)','N(rot)_err',
#                     'Tvib', 'Tvib_err', 'Tvib_LR', 'Tvib_LR_err', 'Trot', 'Trot_err',
#                     'Size', 'hc3nv0_peak_mJy/beam', 'hc3nv7_peak_mJy/beam', 'hc3nv6_peak_mJy/beam',
#                     'tau_v0', 'tau_v0_err', 'tau_v7', 'tau_v7_err', 'tau_v6', 'tau_v6_err']

hc3n_prop['Ind_ok'] = ['14', '13', '10', '11', '12', '8', '9', '5', '4', '3', '2', '1']

#hc3n_prop = hc3n_prop.replace('-', np.nan)
#hc3n_prop['HC3Nv6_219'].replace('-', np.nan)
#hc3n_prop['HC3Nv6_355'].replace('-', np.nan)
#hc3n_prop['HC3Nv6_364'].replace('-', np.nan)
#hc3n_prop['Size_m'] = u_conversion.lin_size(D,hc3n_prop['Size'])
#hc3n_prop['Radio_m'] = hc3n_prop['Size_m']/2.

header_cont = fits.getheader(dworkdir_cubes+'/MAD_CUB_NGC253_TE_219_220_0.19X0.29_briggs_v7.pbcor.fits')#'/MAD_CUB_NGC253_TE_cont_218_HR_briggs.pbcor.fits') 
beam_cont = Beam.from_fits_header(header_cont) 
bmin = beam_cont.minor.to(u.arcsec)
bmaj = beam_cont.major.to(u.arcsec)


# Brightness temperature
T_B = u_conversion.Jybeam_to_T(0.046, 218.373, bmin.value, bmaj.value)
                        

#==============================================================================
# Luminosidades
#==============================================================================
    
# Simulated peaks HC3Nv7=1 219GHz
#hc3n_peak = []
#hc3nv6_peak = []
#for i in hot_cores:
#    hc3n_peak.append(local_max_y_list[i][0][-1])
#    hc3nv6_peak.append(local_max_y_list[i][1][-1])
    
# HC3N v=0
hc3n_prop['HC3Nv0_218_JyBeam'] = hc3n_prop['HC3Nv0_218']/1000.
hc3n_prop['HC3Nv0_354_JyBeam'] = hc3n_prop['HC3Nv0_354']/1000.
freq_v0_218 =   218.324723 # GHz
freq_v0_354 =   354.697463 # GHz

# HC3N v7=1
hc3n_prop['HC3Nv7_219_JyBeam'] = hc3n_prop['HC3Nv7_219']/1000. 
hc3n_prop['HC3Nv7_355_JyBeam'] = hc3n_prop['HC3Nv7_355']/1000. 
freq_v7_219 =   219.173757 # GHz
freq_v7_355 =   355.566254 # GHz

# HC3N v6=1
hc3n_prop['HC3Nv6_219_JyBeam'] = hc3n_prop['HC3Nv6_219']/1000. 
hc3n_prop['HC3Nv6_355_JyBeam'] = hc3n_prop['HC3Nv6_355']/1000. 
hc3n_prop['HC3Nv6_364_JyBeam'] = hc3n_prop['HC3Nv6_364']/1000. 
freq_v6_219 =   218.682561 # GHz
freq_v6_355 =   355.277594 # GHz
freq_v6_364 =   364.380302 # GHz





# Brightness temperature
hc3n_prop['T_B_v0_218'] = u_conversion.Jybeam_to_T(hc3n_prop['HC3Nv0_218_JyBeam'], freq_v0_218, bmin.value, bmaj.value)
hc3n_prop['T_B_v0_354'] = u_conversion.Jybeam_to_T(hc3n_prop['HC3Nv0_354_JyBeam'], freq_v0_354, bmin.value, bmaj.value)

hc3n_prop['T_B_v7_219'] = u_conversion.Jybeam_to_T(hc3n_prop['HC3Nv7_219_JyBeam'], freq_v7_219, bmin.value, bmaj.value)
hc3n_prop['T_B_v7_355'] = u_conversion.Jybeam_to_T(hc3n_prop['HC3Nv7_355_JyBeam'], freq_v7_355, bmin.value, bmaj.value)

hc3n_prop['T_B_v6_219'] = u_conversion.Jybeam_to_T(hc3n_prop['HC3Nv6_219_JyBeam'], freq_v6_219, bmin.value, bmaj.value)
hc3n_prop['T_B_v6_355'] = u_conversion.Jybeam_to_T(hc3n_prop['HC3Nv6_355_JyBeam'], freq_v6_355, bmin.value, bmaj.value)
hc3n_prop['T_B_v6_364'] = u_conversion.Jybeam_to_T(hc3n_prop['HC3Nv6_364_JyBeam'], freq_v6_364, bmin.value, bmaj.value)


# Source size from LTE Tvib (diameter)
hc3n_prop['Source_size_v0_218'] = utiles.sourcesize_from_fit(hc3n_prop['T_B_v0_218'], hc3n_prop['Tvib_HC3N_J24-23'], bmin.value*bmaj.value)
hc3n_prop['Source_size_m_v0_218'] = u_conversion.lin_size(D,hc3n_prop['Source_size_v0_218'])
hc3n_prop['Source_size_v0_354'] = utiles.sourcesize_from_fit(hc3n_prop['T_B_v0_354'], hc3n_prop['Tvib_HC3N_J39-38'], bmin.value*bmaj.value)
hc3n_prop['Source_size_m_v0_354'] = u_conversion.lin_size(D,hc3n_prop['Source_size_v0_354'])

hc3n_prop['Source_size_v7_219'] = utiles.sourcesize_from_fit(hc3n_prop['T_B_v7_219'], hc3n_prop['Tvib_HC3N_J24-23'], bmin.value*bmaj.value)
hc3n_prop['Source_size_m_v7_219'] = u_conversion.lin_size(D,hc3n_prop['Source_size_v7_219'])
hc3n_prop['Source_size_v7_355'] = utiles.sourcesize_from_fit(hc3n_prop['T_B_v7_355'], hc3n_prop['Tvib_HC3N_J39-38'], bmin.value*bmaj.value)
hc3n_prop['Source_size_m_v7_355'] = u_conversion.lin_size(D,hc3n_prop['Source_size_v7_355'])

hc3n_prop['Source_size_v6_219'] = utiles.sourcesize_from_fit(hc3n_prop['T_B_v6_219'], hc3n_prop['Tvib_HC3N_J24-23'], bmin.value*bmaj.value)
hc3n_prop['Source_size_m_v6_219'] = u_conversion.lin_size(D,hc3n_prop['Source_size_v6_219'])
hc3n_prop['Source_size_v6_355'] = utiles.sourcesize_from_fit(hc3n_prop['T_B_v6_355'], hc3n_prop['Tvib_HC3N_J39-38'], bmin.value*bmaj.value)
hc3n_prop['Source_size_m_v6_355'] = u_conversion.lin_size(D,hc3n_prop['Source_size_v6_355'])
hc3n_prop['Source_size_v6_364'] = utiles.sourcesize_from_fit(hc3n_prop['T_B_v7_355'], hc3n_prop['Tvib_HC3N_J39-38'], bmin.value*bmaj.value)
hc3n_prop['Source_size_m_v6_364'] = u_conversion.lin_size(D,hc3n_prop['Source_size_v6_364'])

# Upper limit for source size 0.1"
hc3n_prop['Source_size_m_uplim_v7_219'] = u_conversion.lin_size(D,0.1)

# Line Luminosities
hc3n_prop['L_Watt_v0_218'] = u_conversion.stef_boltz(hc3n_prop['Source_size_m_v0_218']/2., hc3n_prop['Tvib_HC3N_J24-23'])
hc3n_prop['L_Watt_v0_218_err'] = u_conversion.stef_boltz_error(hc3n_prop['Source_size_m_v0_218']/2., 0, hc3n_prop['Tvib_HC3N_J24-23'], hc3n_prop['Tvib_HC3N_J24-23_err'])
hc3n_prop['L_Lsun_v0_218'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v0_218'])
hc3n_prop['L_Lsun_v0_218_err'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v0_218_err'])
hc3n_prop['L_Watt_v0_354'] = u_conversion.stef_boltz(hc3n_prop['Source_size_m_v0_354']/2., hc3n_prop['Tvib_HC3N_J39-38'])
hc3n_prop['L_Watt_v0_354_err'] = u_conversion.stef_boltz_error(hc3n_prop['Source_size_m_v0_354']/2., 0, hc3n_prop['Tvib_HC3N_J39-38'], hc3n_prop['Tvib_HC3N_J39-38_err'])
hc3n_prop['L_Lsun_v0_354'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v0_354'])
hc3n_prop['L_Lsun_v0_354_err'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v0_354_err'])

hc3n_prop['L_Watt_v7_219'] = u_conversion.stef_boltz(hc3n_prop['Source_size_m_v7_219']/2., hc3n_prop['Tvib_HC3N_J24-23'])
hc3n_prop['L_Watt_v7_219_err'] = u_conversion.stef_boltz_error(hc3n_prop['Source_size_m_v7_219']/2., 0, hc3n_prop['Tvib_HC3N_J24-23'], hc3n_prop['Tvib_HC3N_J24-23_err'])
hc3n_prop['L_Lsun_v7_219'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v7_219'])
hc3n_prop['L_Lsun_v7_219_err'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v7_219_err'])
hc3n_prop['L_Watt_v7_355'] = u_conversion.stef_boltz(hc3n_prop['Source_size_m_v7_355']/2., hc3n_prop['Tvib_HC3N_J39-38'])
hc3n_prop['L_Watt_v7_355_err'] = u_conversion.stef_boltz_error(hc3n_prop['Source_size_m_v7_355']/2., 0, hc3n_prop['Tvib_HC3N_J39-38'], hc3n_prop['Tvib_HC3N_J39-38_err'])
hc3n_prop['L_Lsun_v7_355'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v7_355'])
hc3n_prop['L_Lsun_v7_355_err'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v7_355_err'])

hc3n_prop['L_Watt_v6_219'] = u_conversion.stef_boltz(hc3n_prop['Source_size_m_v6_219']/2., hc3n_prop['Tvib_HC3N_J24-23'])
hc3n_prop['L_Watt_v6_219_err'] = u_conversion.stef_boltz_error(hc3n_prop['Source_size_m_v6_219']/2., 0, hc3n_prop['Tvib_HC3N_J24-23'], hc3n_prop['Tvib_HC3N_J24-23_err'])
hc3n_prop['L_Lsun_v6_219'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v6_219'])
hc3n_prop['L_Lsun_v6_219_err'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v6_219_err'])
hc3n_prop['L_Watt_v6_355'] = u_conversion.stef_boltz(hc3n_prop['Source_size_m_v6_355']/2., hc3n_prop['Tvib_HC3N_J39-38'])
hc3n_prop['L_Watt_v6_355_err'] = u_conversion.stef_boltz_error(hc3n_prop['Source_size_m_v6_355']/2., 0, hc3n_prop['Tvib_HC3N_J39-38'], hc3n_prop['Tvib_HC3N_J39-38_err'])
hc3n_prop['L_Lsun_v6_355'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v6_355'])
hc3n_prop['L_Lsun_v6_355_err'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v6_355_err'])
hc3n_prop['L_Watt_v6_364'] = u_conversion.stef_boltz(hc3n_prop['Source_size_m_v6_364']/2., hc3n_prop['Tvib_HC3N_J39-38'])
hc3n_prop['L_Watt_v6_364_err'] = u_conversion.stef_boltz_error(hc3n_prop['Source_size_m_v6_364']/2., 0, hc3n_prop['Tvib_HC3N_J39-38'], hc3n_prop['Tvib_HC3N_J39-38_err'])
hc3n_prop['L_Lsun_v6_364'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v6_364'])
hc3n_prop['L_Lsun_v6_364_err'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v6_364_err'])

# Upper limit for luminosity if ss is 0.1"
hc3n_prop['L_Watt_v7_uplim_219'] = u_conversion.stef_boltz(hc3n_prop['Source_size_m_uplim_v7_219']/2., hc3n_prop['Tvib_HC3N_J24-23'])
hc3n_prop['L_Watt_v7_uplim_219_err'] = u_conversion.stef_boltz_error(hc3n_prop['Source_size_m_uplim_v7_219']/2., 0, hc3n_prop['Tvib_HC3N_J24-23'], hc3n_prop['Tvib_HC3N_J24-23_err'])
hc3n_prop['L_Lsun_v7_uplim_219'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v7_uplim_219'])
hc3n_prop['L_Lsun_v7_uplim_219_err'] = u_conversion.watt_to_lsun(hc3n_prop['L_Watt_v7_uplim_219_err'])


# Size in parsecs (diameter)
hc3n_prop['Source_size_pc_v7_219'] = hc3n_prop['Source_size_m_v7_219'] * 3.2407793E-17

# Virial Mass
hc3n_prop['M_vir_v7_219'] = utiles.virial_mass(hc3n_prop['Source_size_m_v7_219'], hc3n_prop['vlsr'])
                        
# Ratios
hc3n_prop['ratio_24_23'] =hc3n_prop['HC3Nv0_218_JyBeam']/hc3n_prop['HC3Nv7_219_JyBeam']
hc3n_prop['ratio_39_38'] = hc3n_prop['HC3Nv0_354_JyBeam']/hc3n_prop['HC3Nv7_355_JyBeam']
hc3n_prop['ratio_v7'] = hc3n_prop['HC3Nv7_355_JyBeam']/hc3n_prop['HC3Nv7_219_JyBeam']
hc3n_prop['ratio_v0'] = hc3n_prop['HC3Nv0_354_JyBeam']/hc3n_prop['HC3Nv0_218_JyBeam']
    

#hc3n_prop[['ind2', 'Source_size_m_v7_219', 'L_Lsun_v7_219', 'ratio_24_23', 'ratio_39_38', 'ratio_v7', 'ratio_v0']]
                 
model_radius_arcsec = u_conversion.ang_size(D, 3.086e16)


# HC3N properties
hc3nv0_ok = data(lmin='hc3nv0',
              freq=[218.324723*u.GHz, 354.697463*u.GHz],
              elo_cm=[83.7551, 224.8392],
              j_up = [24, 39]
              )

hc3nv0_ok.elo_eV = np.array(hc3nv0_ok.elo_cm) *u.cm**-1 * 1.2398e-4 * u.eV / u.cm**-1
hc3nv0_ok.elo_K = hc3nv0_ok.elo_eV.to(u.K, equivalencies=u.temperature_energy())

hc3nv7_ok = data(lmin='hc3nv7',
              freq=[219.173757*u.GHz, 355.566254*u.GHz],
              elo_cm=[307.0811, 448.391],
              j_up = [24, 39]
              )

hc3nv7_ok.elo_eV = np.array(hc3nv7_ok.elo_cm) *u.cm**-1 * 1.2398e-4 * u.eV / u.cm**-1
hc3nv7_ok.elo_K = hc3nv7_ok.elo_eV.to(u.K, equivalencies=u.temperature_energy())

hc3nv7_2_ok = data(lmin='hc3nv7_2',
              freq=[219.741866*u.GHz, 219.707349*u.GHz, 219.675114*u.GHz, 292.989758*u.GHz, 292.909165*u.GHz, 292.831602*u.GHz],
              elo_cm=[532.566, 532.5596, 530.2792, 599.7559, 599.7355, 597.4417],
              j_up = [24, 24, 24, 32, 32, 32]
              )

hc3nv7_2_ok.elo_eV = np.array(hc3nv7_2_ok.elo_cm) *u.cm**-1 * 1.2398e-4 * u.eV / u.cm**-1
hc3nv7_2_ok.elo_K = hc3nv7_2_ok.elo_eV.to(u.K, equivalencies=u.temperature_energy())


hc3nv6_ok = data(lmin='hc3nv6',
              freq=[218.682561*u.GHz, 355.277594*u.GHz],
              elo_cm=[582.7042, 724.0193],
              j_up = [24, 39]
              )
hc3nv6_ok.elo_eV = np.array(hc3nv6_ok.elo_cm) *u.cm**-1 * 1.2398e-4 * u.eV / u.cm**-1
hc3nv6_ok.elo_K = hc3nv6_ok.elo_eV.to(u.K, equivalencies=u.temperature_energy())


# T_ex = Ts /((Ts/Tkin)+np.log(1+A/C))
# ncrit = A/C
Ts = 838.
Tkin = 450.
Tex = 88.
ncrit = np.exp(Ts*((1./Tex)-(1./Tkin)))-1.

#table_out1 = ['Ind_2','vlsr', 'vlsr_err', 'fwhm', 'fwhm_err',
#             'logN_vib_rdiag_HC3N_J24-23', 'logN_vib_rdiag_HC3N_J24-23', 'Tvib', 'Tvib_err', 'Trot', 'Trot_err',
#             'Source_size_v7', 'hc3nv7_peak_mJy/beam', 'hc3nv6_peak_mJy/beam',
#             'L_Lsun_v7', 'L_Lsun_v7_err']

#hc3n_prop.to_latex(buf=None, columns=table_out1)
