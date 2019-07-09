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
from EC_MS import plot_flux, plot_vs_potential, select_cycles, smooth_data, set_figparams

#from EC_MS import seconds_to_timestamp, timestamp_to_seconds 
from EC_MS import Molecule, ML_strip_cal
#I need to use these because I started the EC file after midnight

plt.close('all')
smoothpoints = 1 #moving-average smoothing for EC data. Set to 1 for no smoothing.

set_figparams(figwidth=9)

importrawdata = True
if importrawdata:
    data_dir = os.path.expanduser('~/Dropbox/Sniffer_Experiments/03_Pt_Sputtered/Data/17I20_sniff2data')
    tag = '09' #just change this to plot something else.
    lslist = os.listdir(data_dir)
    EC_files = [f for f in lslist if f[0:2]==tag and '.mpt' in f]
    MS_file = 'QMS_1.txt'
    MS_data_0 = import_data(data_dir + os.sep + MS_file, data_type='MS')
    EC_datas_0 = []
    for EC_file in EC_files:
        EC_datas_0 += [import_data(data_dir + os.sep + EC_file, data_type='EC')]

EC_datas = [d.copy() for d in EC_datas_0]
# So that I can, later, just select cycles after the pause.  

#syncrhonize in two steps so that we can get the default t-span right.
EC_dataset = synchronize(EC_datas, t_zero='first', override=True)
#EC_dataset['timestamp'] = seconds_to_timestamp(timestamp_to_seconds(EC_dataset['timestamp']) + 24*60*60)
#EC_dataset = dayshift(EC_dataset, days=1) #fixes timestamp, since I started after midnight.
#add a day to the timestamp of the EC_dataset

dataset = synchronize([MS_data_0, EC_dataset], t_zero='start')
dataset['selector'] = dataset['cycle number'] + 10*dataset['file_number']

V_str, J_str = sync_metadata(dataset, RE_vs_RHE=0.687, A_el=0.196)

dataset = smooth_data(dataset, points=smoothpoints, cols=[V_str, J_str]) 

mols_1 = ['He', 'CO']

O2 = Molecule('O2')
CO2 = Molecule('CO2')
H2 = Molecule('H2')
O2.F_cal = 4.7135336141326336 # from fig05.py
CO2.F_cal = 5.2447534530273812
H2.F_cal = 7.1336703492931335

mols_2 = [H2, CO2]

plotfiga = True
dataset['M2-y'] = dataset['M2-y'] - min(dataset['M2-y']) - 3.5e-11
dataset['M44-y'] = dataset['M44-y'] - min(dataset['M44-y']) - 2.5e-12

if plotfiga:

    tspan = [0, 567]
    
    axa = plot_experiment(dataset, mols=mols_1, tspan=tspan, logplot=False, 
                          unit='pmol/s', removebackground=False)
    
    axa += [axa[0].twinx()]
    plot_flux(dataset, mols=mols_2, ax=axa[-1], tspan=tspan, logplot=False, 
              unit='pmol/s', removebackground=False)
    
#    axa[0].set_ylim([-1, 20])
#    axa[0].set_ylabel('MS signal / [nA]\n')
    axa[2].set_ylim([-0.12, 0.14])
    axa[0].set_ylim([-450, 4500])
    axa[3].set_ylim([-1.2, 12])
    axa[1].set_xlim(axa[0].get_xlim())
    axa[1].set_yticks([0, 0.5, 1])
    max_yticks = 4
    yloc = plt.MaxNLocator(max_yticks)
    axa[1].yaxis.set_major_locator(yloc)

    #axa[1].set_xlabel('Testing [mA/cm]')    
    axa[2].set_ylabel('J / [mA cm$^{-2}$]')
    axa[0].set_ylabel('cal. signal / [pmol s$^{-1}$]\n')
    axa[3].set_ylabel('cal. signal / [pmol s$^{-1}$]')

    
    plt.savefig('fig04a.png')
    

plotfigb = True
if plotfigb:
    cycle_CO = select_cycles(dataset, 11, cycle_str='selector') 
    cycle_He = select_cycles(dataset, 12, cycle_str='selector') 
    
    axb = plot_vs_potential(cycle_CO, mols=mols_2, unit='pmol/s', 
                            logplot=False, leg=False, removebackground=False)
    
    plot_vs_potential(cycle_He, mols=mols_2, spec={'linestyle':'--','dashes':(2, 1)}, 
                      unit='pmol/s', ax=axb, logplot=False, leg=False, removebackground=False)
    
    axb[1].set_ylim([-0.09,0.15])
    axb[1].set_yticks([])
    
    max_yticks = 5
    yloc = plt.MaxNLocator(max_yticks)
    axb[1].yaxis.set_major_locator(yloc)
    axb[0].set_ylabel('cal. signal / [pmol s$^{-1}$]')
    
    axb[1].set_ylabel('J / [mA cm$^{-2}$]')

    plt.savefig('fig04b.png')
    
    
plotfigc = True #gives me a problem.
if plotfigc:
    CO2 = Molecule('CO2')
    vspan = [0.65, 1.0]
    calibration, ax = ML_strip_cal(dataset, cycles=[11,12], t_int=100,
                                   cycle_str='selector',
                                 mol=CO2, mass='M44', n_el=2, 
                                 Vspan=vspan, redox='ox',
                                 ax='two', title='fig_4c', 
                                 verbose=True)
    ax[0].set_ylim([-0.10,0.25])
    
    #plt.savefig('fig04c.png',dpi=600)
    
    
    
    