#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 10:34:04 2018

@author: frico
"""

import numpy as np
import astropy.units as u

def Col_dens(a_dust_mod, a_molec1_mod, nh, rnube):
    rnube_tot = (1*u.pc).to(u.cm).value #3.0856776E18 #1 pc in cm
    Nd_tot = 0
    Nmol_tot = 0
    NH_tot = 0
    h = 0
    Nd_list = []
    Nmol_list = []
    NH_list = []
    
    #rnube = np.array([3E17, 6E17,9E17, 1.2E18, 1.5E18, 1.8E18, 2.1E18, 2.4E18, 2.7E18, 2.8E18])
    #nh = [ 28734760.40771883,   7183690.10192971,   3192751.1564132 ,
    #     1795922.52548243,   1149390.41630875,    798187.7891033 ,
    #      586423.68179018,    448980.63137061,    354750.12849036,
    #      329863.32100698]
    #a_dust_mod = [1e-2]*len(rnube)
    #a_molec1_mod = [1e-8]*len(rnube)
    for r, rad in enumerate(rnube):
        if h < (len(rnube)-1):
            Nd_tot_i = (rnube[h+1]-rnube[h])*a_dust_mod[h]*nh[h]
            Nd_tot = Nd_tot + Nd_tot_i
            Nd_list.append(Nd_tot_i)
            
            Nmol_tot_i =( rnube[h+1]-rnube[h])*a_molec1_mod[h]*nh[h]
            Nmol_tot = Nmol_tot + Nmol_tot_i
            Nmol_list.append(Nmol_tot_i)
            
            NH_tot_i =  (rnube[h+1]-rnube[h])*nh[h]
            NH_tot = NH_tot + NH_tot_i
            NH_list.append(NH_tot_i)
        else:
            Nd_tot_i = (rnube_tot-rnube[h])*a_dust_mod[h]*nh[h]
            Nd_tot = Nd_tot + Nd_tot_i
            Nd_list.append(Nd_tot_i)
            
            Nmol_tot_i = (rnube_tot-rnube[h])*a_molec1_mod[h]*nh[h]
            Nmol_tot = Nmol_tot + Nmol_tot_i
            Nmol_list.append(Nmol_tot_i)
            
            NH_tot_i =  (rnube_tot-rnube[h])*nh[h]
            NH_tot = Nmol_tot + NH_tot_i 
            NH_list.append(NH_tot_i)
        h = h +1
        

    print '\n\t NH\t= '+ '%1.3E Ndust' % NH_tot
    print '\t Ndust\t= '+ '%1.3E Ndust' % Nd_tot
    print '\t Nmol\t= '+ '%1.3E Nmol' % Nmol_tot
    return [NH_tot, Nmol_tot, Nd_tot], [Nd_list, Nmol_list, NH_list]

def Mass_from_dens(nh, rnube):
    rnube_tot = 3.086E18
    mass_tot = 0
    h = 0
    mass_list = []
    for r, rad in enumerate(rnube):
        if h < (len(rnube)-1):
            mass_tot_i = (nh[h]*4.*np.pi*(rnube[h+1]-rnube[h])**3)/3.
        else:
            mass_tot_i = (nh[h]*4.*np.pi*(rnube_tot-rnube[h])**3)/3.
            
        mass_tot = mass_tot + mass_tot_i
        mass_list.append(mass_tot_i)
        h = h +1
    
    vol_tot = (4./3)*np.pi*(rnube_tot)**3
    dens_mean = mass_tot/vol_tot
    #print '\n\t Ndust\t= '+ '%1.3E Ndust' % Ndust
    #print '\t Nmol\t= '+ '%1.3E Nmol' % Nmol
    return mass_tot, dens_mean

def abunds_capa(Nd_tot_i, Nmol_tot_i, nh,rnube):
    ### Al dividir el Ntot/10 entre las 10 capas no es real...
    rnube_tot = 3.086E18
    h = 0
    a_dust_modlist = []
    amolec_modlist = []
    amolec_modlist_cte = []
    #rnube = np.array([3E17, 6E17,9E17, 1.2E18, 1.5E18, 1.8E18, 2.1E18, 2.4E18, 2.7E18, 2.8E18])
    #r_percen = rnube/rnube_tot
    #Nd_tot_i = [3E22/10]*len(rnube)
    #Nmol_tot_i = r_percen*1.5E15
    #Nmol_tot_i = [1.5E15/10]*len(rnube)
    #nh = [5E6]*len(rnube)
    for r, rad in enumerate(rnube):
        
        if h < (len(rnube)-1):
            print '\t'+str(r)+'\t'+'%1.2E' % rad
            a_dust_mod = Nd_tot_i[h]/((rnube[h+1]-rnube[h])*nh[h])
            a_dust_modlist.append(a_dust_mod)
            
            a_molec1_mod =Nmol_tot_i[h]/(( rnube[h+1]-rnube[h])*nh[h])
            amolec_modlist.append(a_molec1_mod)
        else:
           print '\tl'+str(r)+'\t'+'%1.2E' % rad
           a_dust_mod = Nd_tot_i[h]/((rnube_tot-rnube[h])*nh[h])
           a_dust_modlist.append(a_dust_mod)
            
           a_molec1_mod = Nmol_tot_i[h]/((rnube_tot-rnube[h])*nh[h])
           amolec_modlist.append(a_molec1_mod)
        amolec_modlist_cte.append(Nmol_tot_i[0]/(( rnube[0+1]-rnube[0])*nh[0]))
          
        
        #print str(h+1)+':\t Xdust\t= '+ '%1.3E ' % a_dust_mod +':\t Xmol\t= '+ '%1.3E ' % a_molec1_mod
        h = h +1
    #amolec_modlist = [1E-8]*len(rnube)
    #[NH_tot, Nmol_tot, Nd_tot], [Nd_list, Nmol_list, NH_list] = Col_dens(a_dust_modlist, amolec_modlist_cte, nh, rnube)
    
    return a_dust_modlist, amolec_modlist


def abunds(Ndust, Nmol, nh,rnube):
    adust = Ndust/(nh*rnube)
    amol = Nmol/(nh*rnube)
    #print '\n\t Xdust\t= '+ '%1.3E adust' % adust
    #print '\t Xmol\t= '+ '%1.3E amol' % amol
    return [adust, amol]

#def densprofile(dens_0, radios, profile):
#    dens = []
#    h = 0 
#    for i, rad in enumerate(radios):
#        if i == 0:
#            dens.append(dens_0)
#        else:
#            dens.append(dens_0 * (radios[0]/rad)**profile)
#        dd = dens[i]
#        #print '\t nH\t= '+ '%1.3E' % dd
#        h = h +1
#    return dens

def densprofile(dens_0, radios, profile):
    """
    dens_0 is the mean final density given a certain profile
    """
    dens = [] # Density profile (1/r)**profile
    rnube_tot = 3.086E18
    #radios = [3E17, 6E17,9E17, 1.2E18, 1.5E18, 1.8E18, 2.1E18, 2.4E18, 2.7E18, 2.8E18]
    vols = [] # Volume of each shell
    v_tot = (4.*np.pi*(rnube_tot**3)/3.)
    for h, rad in enumerate(radios):
        if h < (len(radios)-1):
            dens.append(dens_0 * (radios[0]/radios[h])**(profile))
            vols.append(4.*np.pi*((radios[h+1]**3) - (radios[h]**3))/3.)
        else:
            dens.append(dens_0 * (radios[0]/radios[h])**(profile))
            vols.append(4.*np.pi*((rnube_tot**3) - (radios[h]**3))/3.)
                
    # Mean density (weighting each shell by its volume)
    dens_mean = np.sum(np.array(dens)*np.array(vols))/v_tot
    # Factor to convert current mean dens to desired dens
    ff = dens_0/dens_mean
    # Final densities to give the desired mean density
    ff_dens = np.array(dens)*ff
    # Cheching mean density
    ff_dens_mean = np.sum(np.array(ff_dens)*np.array(vols))/v_tot
    mass_tot = ff_dens_mean * v_tot
    if profile == 0: 
        # Just to keep dens of each shell to be dens_0
        # but actually mean dens is  slightlty <dens_0 
        # (dens dep with r is huge (v=r^-3))
        ff_dens = dens
    return ff_dens, ff_dens_mean, mass_tot

def abunds_capa2(Nd_tot_i, Nmol_tot_i, nh,rnube, profile):
    ### Al dividir el Ntot/10 entre las 10 capas no es real...
    rnube_tot = 3.086E18
    h = 0
    a_dust_modlist = []
    amolec_modlist = []
    amolec_modlist_cte = []
    #rnube = np.array([3E17, 6E17,9E17, 1.2E18, 1.5E18, 1.8E18, 2.1E18, 2.4E18, 2.7E18, 2.8E18])
    #r_percen = rnube/rnube_tot
    #Nd_tot_i = [3E22/10]*len(rnube)
    #Nmol_tot_i = r_percen*1.5E15
    #Nmol_tot_i = [1.5E15/10]*len(rnube)
    #nh = [5E6]*len(rnube)
    for r, rad in enumerate(rnube):
        
        if h < (len(rnube)-1):
            print '\t'+str(r)+'\t'+'%1.2E' % rad
            a_dust_mod = Nd_tot_i[h]/((rnube[h+1]-rnube[h])*nh[h])
            a_dust_modlist.append(a_dust_mod)
            
            a_molec1_mod =Nmol_tot_i[h]/(( rnube[h+1]-rnube[h])*nh[h])
            amolec_modlist.append(a_molec1_mod)
        else:
           print '\tl'+str(r)+'\t'+'%1.2E' % rad
           a_dust_mod = Nd_tot_i[h]/((rnube_tot-rnube[h])*nh[h])
           a_dust_modlist.append(a_dust_mod)
            
           a_molec1_mod = Nmol_tot_i[h]/((rnube_tot-rnube[h])*nh[h])
           amolec_modlist.append(a_molec1_mod)
        amolec_modlist_cte.append(Nmol_tot_i[0]/(( rnube[0+1]-rnube[0])*nh[0]))
          
        
        #print str(h+1)+':\t Xdust\t= '+ '%1.3E ' % a_dust_mod +':\t Xmol\t= '+ '%1.3E ' % a_molec1_mod
        h = h +1
    #amolec_modlist = [1E-8]*len(rnube)
    #[NH_tot, Nmol_tot, Nd_tot], [Nd_list, Nmol_list, NH_list] = Col_dens(a_dust_modlist, amolec_modlist_cte, nh, rnube)
    
    return a_dust_modlist, amolec_modlist

def denstot(dens_0, radios):
    dens = []
    h = 0
    for i, rad in enumerate(radios):
        if i == 0:
            dens.append(dens_0)
        else:
            dens.append(dens_0 * (radios[0]/rad)**2)
        dd = dens[i]
        #print '\t nH\t= '+ '%1.3E' % dd
        h = h +1
    return dens

def tdust_devicente(luminosity, radio, gradient, td):
    """
    Returns Tdust given luminosity in Lsun and radio in cm
    Eq.1 From de Vicente et al 2000
    """
    if gradient == True:
        Tdust = 15.222 *(luminosity * (1e16/radio)**2.)**(1/4.)  
    else:
        Tdust = td
    return Tdust