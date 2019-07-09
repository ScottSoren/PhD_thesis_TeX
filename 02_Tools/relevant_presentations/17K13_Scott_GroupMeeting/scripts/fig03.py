#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 23:16:46 2017

This has been edited by Scott 17K13 to plot calibrated signals

@author: cinf
"""

import os
#import numpy as np
from matplotlib import pyplot as plt
from EC_MS import import_data, sync_metadata, synchronize, plot_experiment
from EC_MS import plot_signal, plot_vs_potential, select_cycles, smooth_data
from EC_MS import set_figparams, is_time, Molecule, plot_flux

plt.close('all')
smoothpoints = 20 #moving-average smoothing for EC data. Set to 1 for no smoothing.

set_figparams(figwidth=9)

importrawdata = True
if importrawdata:
    data_dir = os.path.expanduser('~/Dropbox/Sniffer_Experiments/03_Pt_Sputtered/Data/17I20_sniff2data')
    tag = '19'
    lslist = os.listdir(data_dir)
    EC_files = [f for f in lslist if f[0:2]==tag and '.mpt' in f]
    MS_file = 'QMS_2.txt'
    MS_data_0 = import_data(data_dir + os.sep + MS_file, data_type='MS')
    EC_datas_0 = []
    for EC_file in EC_files:
        EC_datas_0 += [import_data(data_dir + os.sep + EC_file, data_type='EC')]

#syncrhonize in two steps so that we can get the default t-span right.
EC_dataset = synchronize(EC_datas_0, t_zero='first', override=True)
dataset = synchronize([MS_data_0, EC_dataset], t_zero='start')

V_str, J_str = sync_metadata(dataset, RE_vs_RHE=0.687, A_el=0.196)

normalize_to_H2O = False
if normalize_to_H2O:
    normalizer = dataset['M18-y'][19500] / dataset['M18-y']  # normilization value taken just before onset of HER in sniff2fig03
    #print(normalizer)
    #plt.plot(np.arange(len(dataset['M18-y'])),dataset['M18-y'])
    
    for col in dataset['data_cols']:
        if col[-2:]=='-y':
            #print(col)
            #print(dataset[col])
            dataset[col] = dataset[col] * normalizer

dataset = smooth_data(dataset, points=smoothpoints, cols=[V_str, J_str]) 
'''
F_cal_O2
Out[23]: 4.7135336141326336

F_cal_CO2
Out[24]: 5.2447534530273812

F_cal_H2
Out[25]: 7.1336703492931335
'''
mols_1 = ['He', 'CO']

O2 = Molecule('O2')
CO2 = Molecule('CO2')
H2 = Molecule('H2')
O2.F_cal = 4.7135336141326336 # from fig05.py
CO2.F_cal = 5.2447534530273812
H2.F_cal = 7.1336703492931335

He = Molecule('He')
CO = Molecule('CO')



mols_2 = [H2, O2, CO2]

plotfiga = True
if plotfiga:

    tspan = [0, 1800]
    
    axa = plot_experiment(dataset, mols=mols_1, tspan=tspan, logplot=False, unit='pmol/s')
    
    axa += [axa[0].twinx()]
    plot_flux(dataset, mols=mols_2, ax=axa[-1], tspan=tspan, logplot=False, unit='pmol/s')
     
    #axa[0].set_ylim([-1, 20])
    #axa[0].set_ylabel('MS signal / [nA]\n')
    axa[1].set_ylim([-0.1, 1.9])
    axa[2].set_ylim([-0.5, 1])
    #axa[3].set_ylim([-0.1, 2.0])
    axa[0].set_ylim([-450, 4500])  
    axa[3].set_ylim([-45, 450])      
    axa[0].set_ylabel('cal. signal / [pmol s$^{-1}$]\n')
    axa[3].set_ylabel('cal. signal / [pmol s$^{-1}$]')
    
    max_yticks = 5
    yloc = plt.MaxNLocator(max_yticks)
    axa[1].set_yticks([0, 0.5, 1, 1.5])
    axa[2].yaxis.set_major_locator(yloc)
    
    max_xticks = 6
    xloc = plt.MaxNLocator(max_xticks)
    axa[0].xaxis.set_major_locator(xloc)
    axa[2].xaxis.set_major_locator(xloc)
    axa[1].set_xlim(axa[0].get_xlim())    
    axa[2].set_ylabel('J / [mA cm$^{-2}$]')
     
    plt.savefig('fig03a.png')
    


plotfigb = True
if plotfigb:
    cycle_He = select_cycles(dataset, 2, cycle_str='cycle number') 
    cycle_CO = select_cycles(dataset, 3, cycle_str='cycle number') 
    
    axb = plot_vs_potential(cycle_He, mols=mols_2, unit='pmol/s', 
                            logplot=False, leg=False)
    
    plot_vs_potential(cycle_CO, mols=mols_2, spec={'linestyle':'--','dashes':(2, 1)}, 
                      unit='pmol/s', ax=axb, logplot=False, leg=False)
    
    max_yticks = 5
    yloc = plt.MaxNLocator(max_yticks)
    axb[1].yaxis.set_major_locator(yloc)

    #axb[1].set_ylim([-0.1,1.9])
    #axb[2].set_ylim([-0.5, 1])
 
    #axb[0].set_ylim([-0.1, 2.0])
    axb[1].set_ylim([-0.75, 1.499])
    
    max_xticks = 8
    xloc = plt.MaxNLocator(max_xticks)
    axb[0].xaxis.set_major_locator(xloc)
    axb[1].xaxis.set_major_locator(xloc)
    
    #axb[0].set_ylabel('MS signal / [nA]\n')
    
    axb[1].set_ylabel('J / [mA cm$^{-2}$]')
    axb[0].set_ylabel('cal. signal / [pmol s$^{-1}$]')
    
    axb[0].set_xlim([-0.05, 1.85])
    axb[1].set_xlim([-0.05, 1.85])

    plt.savefig('fig03b.png')
    
    
    
    
    
    
    
plot_sub_panels = True
if plot_sub_panels:
#    os.chdir('./fig03_sub_panels_Daniel')
    
    #sub panel a    
    masses = ['M2',
              #'M32',
              #'M44',
              ]
    
#    for col in dataset['data_cols']:
#        if is_time(col):
#            dataset[col] -= 62
              
    tspan = [62, 362]
    
    axa = plot_experiment(dataset, masses=masses, tspan=tspan, logplot=False, unit='pmol/s')
       
    axa[0].set_ylim([-0.5, 5.5])
    axa[1].set_ylim([-0.05,0.55])
    axa[2].set_ylim([-0.5, 0.5])
    axa[2].set_ylabel('J / [mA cm$^{-2}$]')
    
    axa[0].set_xlim([62, 312])
    axa[1].set_xlim([62, 312])
    
    max_yticks = 5
    yloc = plt.MaxNLocator(max_yticks)
    axa[1].yaxis.set_major_locator(yloc)

    max_yticks = 6
    yloc = plt.MaxNLocator(max_yticks)
    axa[2].yaxis.set_major_locator(yloc)
     
    #plt.savefig('fig03a.png')
    
#    for col in dataset['data_cols']:
#        if is_time(col):
#            dataset[col] += 62

    #sub panel b    
    masses = [#'M2',
              'M32',
              #'M44',
              ]
    
#    for col in dataset['data_cols']:
#        if is_time(col):
#            dataset[col] -= 57 + 250
              
    tspan = [307, 607]
    
    axb = plot_experiment(dataset, masses=masses, tspan=tspan, logplot=False, unit='pmol/s')
       
    axb[0].set_ylim([-0.5, 3.5])
    axb[1].set_ylim([0.1,2.1])
    axb[2].set_ylim([-0.5, 1.0])
    axb[2].set_ylabel('J / [mA cm$^{-2}$]')
    
    axb[0].set_xlim([307, 557])
    axb[1].set_xlim([307, 557])
    
    max_yticks = 5
    yloc = plt.MaxNLocator(max_yticks)
    axb[0].yaxis.set_major_locator(yloc)
    
    max_yticks = 5
    yloc = plt.MaxNLocator(max_yticks)
    axb[1].yaxis.set_major_locator(yloc)

    max_yticks = 5
    yloc = plt.MaxNLocator(max_yticks)
    axb[2].yaxis.set_major_locator(yloc)
     
    #plt.savefig('fig03b.png')
    
#    for col in dataset['data_cols']:
#        if is_time(col):
#            dataset[col] += 57 + 250

    #sub panel c    
    masses_1 = ['M4', 'M28']
    masses_2 = ['M2', 'M32', 'M44']
    
    tspan = [500, 1300]
    
    axc = plot_experiment(dataset, mols=mols_1, tspan=tspan, logplot=False, unit='pmol/s')
    
    axc += [axc[0].twinx()]
    plot_flux(dataset, mols=mols_2, ax=axc[-1], tspan=tspan, logplot=False, unit='pmol/s')
     
#    axc[0].set_ylim([-2.5, 82.5])
    #axc[0].set_ylabel('MS signal / [nA]\n')
    axc[1].set_ylim([-0.1,1.9])
    axc[2].set_ylim([-0.5, 1])
    axc[3].set_ylim([-35, 350])
    
    max_yticks = 6
    yloc = plt.MaxNLocator(max_yticks)
    axc[0].yaxis.set_major_locator(yloc)
    
    max_yticks = 5
    yloc = plt.MaxNLocator(max_yticks)
    axc[2].yaxis.set_major_locator(yloc)
    
    max_xticks = 6
    xloc = plt.MaxNLocator(max_xticks)
    axc[0].xaxis.set_major_locator(xloc)
    axc[2].xaxis.set_major_locator(xloc)
    
    axc[2].set_ylabel('J / [mA cm$^{-2}$]')
    axc[1].set_xlim(axc[0].get_xlim())     
    plt.savefig('fig03c.png')
    
   #sub panel d   
    masses = [#'M2',
              #'M32',
              'M44',
              ]
    
#    for col in dataset['data_cols']:
#        if is_time(col):
#            dataset[col] -= 40 + 250 + 250 + 800
              
    tspan = [1340, 1640]
    
    axd = plot_experiment(dataset, mols=CO2, tspan=tspan, logplot=False, unit='pmol/s')
       
    axd[0].set_ylim([-0.5, 8.5])
    axd[1].set_ylim([0.3,1.2])
    axd[2].set_ylim([-0.5, 1.0])
    axd[2].set_ylabel('J / [mA cm$^{-2}$]')
    
    axd[0].set_xlim([1340, 1590])
    axd[1].set_xlim([1340, 1590])
    
    max_yticks = 5
    yloc = plt.MaxNLocator(max_yticks)
    axd[1].yaxis.set_major_locator(yloc)

    max_yticks = 5
    yloc = plt.MaxNLocator(max_yticks)
    axd[2].yaxis.set_major_locator(yloc)
     
    plt.savefig('fig03d.png')
    
#    for col in dataset['data_cols']:
#        if is_time(col):
#            dataset[col] += 40 + 250 + 250 + 800