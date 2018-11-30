#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 11:54:29 2018

@author: frico
"""

import matplotlib.pyplot as plt
from ferpy import utiles
import numpy as np
import os


# Data workdir
dworkdir_cubes = '/Users/frico/Documents/data/NGC253_H3O+'
dworkdir_spec = dworkdir_cubes+'/Hotcores_v4_all'


# Out dir
out_dir = dworkdir_spec+'/Results_v8/'
if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
# Out figura 1
out_fig1_dir = out_dir+'/figura2/'
if not os.path.exists(out_fig1_dir):
        os.makedirs(out_fig1_dir)

v0_ground = 0
v7_ground = 223.3042
v6_ground = 499.1157
v72_ground = 445.9987

v0_info = {'24-23': {'frec': 218.32, 'ELO':120.50},
           '39-38': {'frec': 354.70, 'ELO':323.48}
           }

v7_1_info = {'24-23': {'frec': 219.17, 'ELO':441.81},
             '39-38': {'frec': 355.57, 'ELO':645.11}
             }

v6_1_info = {'24-23': {'frec': 218.68, 'ELO':838.35},
             '39-38': {'frec': 355.28, 'ELO':1041.67}
             }

v7_2_info = {'24-23a': {'frec': 219.74, 'ELO':766.22},
             '24-23b': {'frec': 219.71, 'ELO':766.21},
             '24-23c': {'frec': 219.68, 'ELO':762.93},
             '32-31a': {'frec': 292.99, 'ELO':862.89},
             '32-31b': {'frec': 292.91, 'ELO':862.86},
             '32-31c': {'frec': 292.93, 'ELO':859.56}
             }




fig, ax = plt.subplots()

# Hide the right bottom and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)

# Only show ticks on the left spines
ax.yaxis.set_ticks_position('left')
# Hide Xaxis ticks
ax.get_xaxis().set_visible(False)
# Hide Yaxis ticks
#ax.get_yaxis().set_visible(False)

label_sep = 50
# Plotting levels
posv0 = [3,4]

# Ground State vib level
ax.plot([posv0[0]-0.5, posv0[1]], [v0_ground,v0_ground],color='k', linewidth=2, zorder=1)
ax.text(np.mean([posv0[0]-0.5, posv0[1]]), v0_ground-label_sep, r'$v=0$', ha='center', va='center', backgroundcolor='none', fontsize=8 , color='k', fontweight='bold')
# v=0
for key in v0_info:
    labels = key.split('-')
    # lower level
    ax.plot(posv0, [v0_info[key]['ELO'],v0_info[key]['ELO']],color='k', linewidth=1)
    ax.text(np.mean(posv0)-0.6, v0_info[key]['ELO']-10, r'$'+labels[1]+'$', ha='right', va='center', backgroundcolor='none', fontsize=6, color='k', fontweight='bold')

    # upper level
    freq_E = utiles.trans_energy(v0_info[key]['frec'])+30
    ax.plot(posv0, [v0_info[key]['ELO']+freq_E,v0_info[key]['ELO']+freq_E],color='k', linewidth=1)
    ax.text(np.mean(posv0)-0.6, v0_info[key]['ELO']+freq_E+10, r'$'+labels[0]+'$', ha='right', va='center', backgroundcolor='none', fontsize=6 , color='k', fontweight='bold')
    # v=0 -> v=0
    #ax.annotate("", xy=(np.mean([posv0[0]-0.5, posv0[1]]), v0_info[key]['ELO']+freq_E+5),
    #            xytext=(np.mean([posv0[0]-0.5, posv0[1]]), v0_ground),
    #            arrowprops=dict(arrowstyle="->, head_length=0.2, head_width=0.1", fc='red', ec='red', shrinkA = 0, shrinkB = 0, linewidth=0.7), zorder=2) 
    # Rotational
    ax.annotate("", xy=(np.mean(posv0)+0.0, v0_info[key]['ELO']),
                    xytext=(np.mean(posv0)+0.0, v0_info[key]['ELO']+freq_E),
                    arrowprops=dict(arrowstyle="->, head_length=0.2, head_width=0.1",
                                #connectionstyle="arc3,rad=-4",
                                fc='green', ec='green', shrinkA = 0, shrinkB = 0, linewidth=0.7)
                                , zorder=2) 
# v7=1
posv7 = [5,6]

ax.plot([posv7[0]-0.5, posv7[1]], [v7_ground,v7_ground],color='k', linewidth=2, zorder=1)
ax.text(np.mean([posv7[0]-0.5, posv7[1]]), v7_ground-label_sep, r'$v7=1$', ha='center', va='center', backgroundcolor='none', fontsize=8 , color='k', fontweight='bold')
# v=0 -> v7=1
ax.annotate("", xy=(np.mean([posv7[0]-0.5, posv7[1]]), v7_ground),
                xytext=(np.mean([posv0[0]-0.5, posv0[1]]), v0_ground),
                arrowprops=dict(arrowstyle="->, head_length=0.2, head_width=0.1", fc='red', ec='red', shrinkA = 0, shrinkB = 0, linewidth=0.7), zorder=2)
for key in v7_1_info:
    labels = key.split('-')
    # lower level
    ax.plot(posv7, [v7_1_info[key]['ELO'],v7_1_info[key]['ELO']],color='k', linewidth=1)
    ax.text(np.mean(posv7)-0.6, v7_1_info[key]['ELO']-10, r'$'+labels[1]+'$', ha='right', va='center', backgroundcolor='none', fontsize=6, color='k', fontweight='bold')

    # upper level
    freq_E = utiles.trans_energy(v7_1_info[key]['frec'])+30
    ax.plot(posv7, [v7_1_info[key]['ELO']+freq_E,v7_1_info[key]['ELO']+freq_E],color='k', linewidth=1)
    ax.text(np.mean(posv7)-0.6, v7_1_info[key]['ELO']+freq_E+10, r'$'+labels[0]+'$', ha='right', va='center', backgroundcolor='none', fontsize=6 , color='k', fontweight='bold')
    # v=0 -> v7=1
    #ax.annotate("", xy=(np.mean([posv7[0]-0.5, posv7[1]]), v7_1_info[key]['ELO']+freq_E+5),
    #            xytext=(np.mean([posv0[0]-0.5, posv0[1]]), v0_ground),
    #            arrowprops=dict(arrowstyle="->, head_length=0.2, head_width=0.1", fc='red', ec='red', shrinkA = 0, shrinkB = 0, linewidth=0.7), zorder=2) 
    
    # Rotational
    ax.annotate("", xy=(np.mean(posv7)+0.0, v7_1_info[key]['ELO']),
                    xytext=(np.mean(posv7)+0.0, v7_1_info[key]['ELO']+freq_E),
                    arrowprops=dict(arrowstyle="->, head_length=0.2, head_width=0.1",
                                #connectionstyle="arc3,rad=-4",
                                fc='green', ec='green', shrinkA = 0, shrinkB = 0, linewidth=0.7)
                                , zorder=2) 
    
    
    
#v6=1
posv6 = [1,2]
# v=0 -> v6=1
ax.annotate("", xy=(np.mean([posv6[0]-0.5, posv6[1]]), v6_ground),
                xytext=(np.mean([posv0[0]-0.5, posv0[1]]), v0_ground),
                arrowprops=dict(arrowstyle="->, head_length=0.2, head_width=0.1", fc='red', ec='red', shrinkA = 0, shrinkB = 0, linewidth=0.7), zorder=2) 
ax.plot([posv6[0]-0.5, posv6[1]], [v6_ground,v6_ground],color='k', linewidth=2, zorder=1)
ax.text(np.mean([posv6[0]-0.5, posv6[1]]), v6_ground-label_sep, r'$v6=1$', ha='center', va='center', backgroundcolor='none', fontsize=8 , color='k', fontweight='bold')
for key in v6_1_info:
    labels = key.split('-')
    # lower level
    ax.plot(posv6, [v6_1_info[key]['ELO'],v6_1_info[key]['ELO']],color='k', linewidth=1)
    ax.text(np.mean(posv6)-0.6, v6_1_info[key]['ELO']-10, r'$'+labels[1]+'$', ha='right', va='center', backgroundcolor='none', fontsize=6, color='k', fontweight='bold')

    # upper level
    freq_E = utiles.trans_energy(v6_1_info[key]['frec'])+30
    ax.plot(posv6, [v6_1_info[key]['ELO']+freq_E,v6_1_info[key]['ELO']+freq_E],color='k', linewidth=1)
    ax.text(np.mean(posv6)-0.6, v6_1_info[key]['ELO']+freq_E+10, r'$'+labels[0]+'$', ha='right', va='center', backgroundcolor='none', fontsize=6 , color='k', fontweight='bold')
    
    # Rotational
    ax.annotate("", xy=(np.mean(posv6)+0.0, v6_1_info[key]['ELO']),
                    xytext=(np.mean(posv6)+0.0, v6_1_info[key]['ELO']+freq_E),
                    arrowprops=dict(arrowstyle="->, head_length=0.2, head_width=0.1",
                                #connectionstyle="arc3,rad=-4",
                                fc='green', ec='green', shrinkA = 0, shrinkB = 0, linewidth=0.7)
                                , zorder=2) 
#v7=2
posv72 = [7,8]
# v7=1 -> v7=2
ax.annotate("", xy=(np.mean([posv72[0]-0.5, posv72[1]]), v72_ground),
                xytext=(np.mean([posv7[0], posv7[1]]), v7_ground),
                arrowprops=dict(arrowstyle="->, head_length=0.2, head_width=0.1", fc='red', ec='red', shrinkA = 0, shrinkB = 0, linewidth=0.7), zorder=2)
# v7=2 low
ax.plot([posv72[0]-0.5, posv72[1]], [v72_ground,v72_ground],color='k', linewidth=2, zorder=1)
ax.text(np.mean([posv72[0]-0.5, posv72[1]]), v72_ground-label_sep, r'$v7=2$', ha='center', va='center', backgroundcolor='none', fontsize=8 , color='k', fontweight='bold')
for key in v7_2_info:
    if key in ['24-23a', '32-31a']:
        labels = key.split('-')
        # lower level
        ax.plot(posv72, [v7_2_info[key]['ELO'],v7_2_info[key]['ELO']],color='k', linewidth=1)
        ax.text(np.mean(posv72)-0.6, v7_2_info[key]['ELO']-10, r'$'+labels[1][:-1]+'$', ha='right', va='center', backgroundcolor='none', fontsize=6, color='k', fontweight='bold')
        # upper level
        freq_E = utiles.trans_energy(v7_2_info[key]['frec'])+30
        ax.plot(posv72, [v7_2_info[key]['ELO']+freq_E,v7_2_info[key]['ELO']+freq_E],color='k', linewidth=1)
        ax.text(np.mean(posv72)-0.6, v7_2_info[key]['ELO']+freq_E+10, r'$'+labels[0]+'$', ha='right', va='center', backgroundcolor='none', fontsize=6 , color='k', fontweight='bold')
         

        # Rotational
        ax.annotate("", xy=(np.mean(posv72)+0.0, v7_2_info[key]['ELO']),
                    xytext=(np.mean(posv72)+0.0, v7_2_info[key]['ELO']+freq_E),
                    arrowprops=dict(arrowstyle="->, head_length=0.2, head_width=0.1",
                                #connectionstyle="arc3,rad=-4",
                                fc='green', ec='green', shrinkA = 0, shrinkB = 0, linewidth=0.7)
                                , zorder=2)

#ax.set_yscale('symlog')    
ax.set_xlim(0,9)
ax.set_ylim(-10,1100)
ax.set_ylabel('Energy (K)')#, fontsize=labelsize, fontname=labelfont, fontweight='bold')



fig.savefig(out_fig1_dir+'/hc3n_Ediag_v2.eps', bbox_inches='tight', transparent=True, dpi=400)
plt.close()