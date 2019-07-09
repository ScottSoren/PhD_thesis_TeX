#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 15:22:42 2017

@author: scott@fysik.dtu.dk
EC_MS package available at https://github.com/ScottSoren/EC_MS
"""

import os
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from EC_MS import import_folder, sync_metadata, plot_experiment, plot_flux
from EC_MS import Molecule, cut, get_signal, get_flux, Chem, is_time
from EC_MS import standard_colors, select_cycles, plot_vs_potential
from EC_MS import clip_cycles, stagnant_operator, Chip

plt.close('all')

A_el = 0.196
RE_vs_RHE = 0.747
t_str = 'time/s'
J_str = 'J /[mA/cm^2]'
V_str = 'U vs RHE / [V]'
        
importrawdata = True
if importrawdata:
    data_directory = os.path.expanduser('~/Dropbox/Sniffer_Experiments/02_NiFeNPs/Data')
    folder = '17J05_Pt_isotope_exchange'
    datasets = import_folder(data_directory + os.sep + folder, tags='all')

plotrawdata = True
if plotrawdata:
    for tag, data in datasets.items():
        sync_metadata(data, RE_vs_RHE=RE_vs_RHE, A_el=A_el)
        data['J / [uA/cm^2]'] = data[J_str].copy() * 1e3  # for plotting uA/cm^2 # the copy() is necessary
        #data['J_str'] = J_str # plot uA/cm^2 by default
        t_str = 'time/s'
        ax = plot_experiment(data, masses=['M2','M4','M28','M32','M34','M36','M44','M46','M48'])
        ax[1].set_title(tag)
        plt.savefig('COox_17J05_' + tag + '.png')

calibrateit = True
if calibrateit:
    data = datasets['03']
    mols = {}
    
    tspan_CO2 = [2150, 2200]
    t, j = cut(data[t_str], data[J_str], tspan_CO2)
    n_CO2 = np.trapz(j, t) * A_el * 1e-3 / (2 * Chem.Far) 
    Y = {}
    for mass in ['M44', 'M46', 'M48']:
        x, y = get_signal(data, mass, tspan=tspan_CO2, unit='mol/s')
        Y[mass] = np.trapz(y, x)
    Y_CO2 = sum(Y.values())
    F_cal_CO2 = Y_CO2 / n_CO2
    for mass in ['M44', 'M46', 'M48']:
        m = Molecule('CO2')
        m.primary = mass
        m.F_cal = F_cal_CO2
        mols[mass] = m
        
    tspan_O2 = [1320, 1380]
    t, j = cut(data[t_str], data[J_str], tspan_O2)
    n_O2 = np.trapz(j, t) * A_el * 1e-3 / (4 * Chem.Far) 
    Y = {}
    for mass in ['M32', 'M34', 'M36']:
        x, y = get_signal(data, mass, tspan=tspan_O2, unit='mol/s')
        Y[mass] = np.trapz(y, x)
    Y_O2 = sum(Y.values())
    F_cal_O2 = Y_O2 / n_O2
    for mass in ['M32', 'M34', 'M36']:
        m = Molecule('O2')
        m.primary = mass
        m.F_cal = F_cal_O2
        mols[mass] = m

    #Gives me 4.2 C/mol and 5.5 C/mol, respectively, within 10% of what's in the Snfifer 2 paper. Nice!

fitL = True
if fitL:
    data = datasets['03']
    tspan_O2tail = [1813, 1913]
    x,y = cut(data['M34-x'], data['M34-y'], tspan=tspan_O2tail)
    x = x-x[0]
    y = y - y[-1]
    y = y/max(y) #normalize
    chip = Chip('SI-3iv1-1-C5', l_cap=3.33e-3)
    q0 = chip.capillary_flow(gas='CO')/Chem.NA  
    L = 160e-6
    results = stagnant_operator(j_el=1, tpulse=40, tspan=[0, 70],
                                #startstate = 'steady is not working!
                                L=L, q0=q0, normalize=True,
                                mol=mols['M32'], n_el=4, )  
    j = results['j']/max(results['j'])
    t = results['t']
    fig, ax9 = plt.subplots()
    ax9.plot(x, y,'k-')
    ax9.plot(t, j, 'k--')


modelCV = True
if modelCV:

    data = datasets['02']
    # the cycles I want are close to 14 and 17, but the cycle switch is at an unfortunate place
    cycles_CO = select_cycles(data, cycles=[13, 14])
    cycles_He = select_cycles(data, cycles=[16, 17])
    cycle_CO = clip_cycles(cycles_CO, cycles=1, V_clip=0.45, redox=1, verbose=True)
    cycle_He = clip_cycles(cycles_He, cycles=1, V_clip=0.45, redox=1, verbose=True)
    
    plotit = True
    if plotit:
        tspan_CVs = [10400, 11650]

        ax0 = plot_experiment(data, mols=['He','CO'], tspan=tspan_CVs,  #J_str='cycle number',
                              unit='pmol/s',logplot=False,)
        ax0 += [ax0[0].twinx()]
        plot_flux(data, ax=ax0[3], mols=list(mols.values()), tspan=tspan_CVs, 
                  unit='pmol/s', logplot=False)
        ax0[1].set_xlim(ax0[0].get_xlim())
        plt.savefig('COox_CVs.png')

        ax1 = plot_vs_potential(cycle_CO, mols=list(mols.values()), logplot=False)
        plot_vs_potential(cycle_He, mols=list(mols.values()), ax=ax1, spec={'linestyle': '--'}, logplot=False)
        plt.savefig('COox_CVs_vs_potential.png')    
    
    predictsignal = True
    if predictsignal:
        N_points = 2700 #min([len(cycle_CO[t_str]), len(cycle_He[t_str])])
        t_zero = cycle_CO[t_str][0]
        t = cycle_CO[t_str][:N_points] - t_zero
        j_diff = cycle_CO[J_str][:N_points] - cycle_He[J_str][:N_points]

        #lets make at dataset dictionary to store what we need for this model:
        diff = {}
        diff[J_str] = j_diff
        diff[t_str] = t
        for mass in ['M44','M46','M48']:
            diff[mass + '-x'] = cycle_CO[mass + '-x'] - t_zero
            diff[mass + '-y'] = cycle_CO[mass + '-y']
        # and a whole bunch of bullshit it needs to run properly:
        diff['J_str'] = J_str
        diff['V_str'] = V_str
        diff['title'] = 'title'
        diff['A_el'] = A_el
        tspan_model = [0, 350]
        ax3 = plot_experiment(diff, mols=[mols['M44'], mols['M46'], mols['M48']],
                              tspan=tspan_model, plotpotential=False, logplot=False) 
        ax3[1].set_xlim(ax3[0].get_xlim())
        plt.savefig('COox_CVs_predicted1.png')
        
        runsimulation = True
        if runsimulation:
            chip = Chip('SI-3iv1-1-C5', l_cap=3.33e-3)
            q0 = chip.capillary_flow(gas='CO')/Chem.NA  
            j_el = j_diff * 10 # so that it's in A/m^2
            #need to fit L (above) to get a good model. Turns out we were way indented! L=160um
            results = stagnant_operator(tj = [t, j_el], tspan=tspan_model, L=L,
                                        mol=mols['M44'], n_el=2, q0=q0)
        ax3[0].plot(results['t'], results['j']/3, 'k--')
        plt.savefig('COox_CVs_predicted2.png')
        
        ax3[0].clear()  
        for mass in ['M44','M46','M48']:
            x, y = diff[mass + '-x'], diff[mass + '-y']
            y_norm = (y-min(y))/(max(y)-min(y)) 
            ax3[0].plot(x, y_norm, standard_colors[mass])
      
        ax3[0].plot(results['t'], results['j']/max(results['j']), 'k--')
        ax3[0].set_xlim(ax3[1].get_xlim())
        ax3[0].set_ylabel('normalized flux')
        plt.savefig('COox_CVs_predicted_normalized.png')    
    
definefunctions = True
if definefunctions:  #to be used later, when fitting the carbonic rate constant
    def carbonic_ode(S, t, pars):
        '''
        Equations from Scott's notebook 17K11
        '''
        k = pars[0] # rate constant / s^-1
        g = pars[1] # H2(16)O / H2(18)O ratio
        
        S44, S46, S48 = S[0], S[1], S[2]
        
        dS44 = k * (-2/3*(1-g)*S44 + 1/3*g*S46)
        dS46 = k * (2/3*(1-g)*S44 - 1/3*S46 + 2/3*g*S48)
        dS48 = k * (1/3*(1-g)*S46 - 2/3*g*S48)
        
        return [dS44, dS46, dS48]
    
    def solve_carbonic_burst(k=0.026, g=0.27, tspan=[0,60]):
        '''
        Returns the partial concentrations at M44, M46, and M48 following an
        initial burst of CO(16) oxidation given:
            g = the H2(18)O/H2(16)O ratio
            k = the rate constant for H2O + CO2 --> H2CO3 in s^-1
        '''
        print('k = ' + str(k))
        print('g = ' + str(g))
        S0 = np.array([g, 1-g, 0])
        pars = [k, g]
        if len(tspan) == 2:
            tspan = np.linspace(tspan[0], tspan[-1], 200)
        SS = odeint(carbonic_ode, S0, tspan, args=(pars,))    
        return SS

calculatecarbonicrate = True
if calculatecarbonicrate:
    data = datasets['02'].copy()
    for col in data['data_cols']:
        if is_time(col):
            data[col] = data[col].copy() - 18621
    data[t_str] = data[t_str] + 0.8 #so the datasets look more aligned. We need triggers!
    plotit = True
    if plotit: 
        ax1 = plot_experiment(data, mols=['He','CO'], unit='pmol/s', tspan=[-600, 200], logplot=False)
        ax1 += [ax1[0].twinx()]
        plot_flux(data, mols=list(mols.values()), ax=ax1[3], unit='pmol/s', tspan=[-600, 200], logplot=False)
        ax1[2].set_ylim([-0.040, 0.4])
        plt.savefig('COox_faststrip.png')
        ax2 = plot_experiment(data, mols=[mols['M44'], mols['M46'], mols['M48']], tspan=[-2, 60], 
                              logplot=False, J_str='J / [uA/cm^2]')
        ax2[2].set_ylim([-20, 75])
        plt.savefig('COox_faststrip_zoom.png')
    
    #get water isotope ratio from O2 signals
    getgamma = True
    if getgamma:
        tspan_O2 = [50, 100]
        for mass in ['M32', 'M34', 'M36']:
            x, y = get_signal(data, mass, tspan=tspan_O2, unit='nA', removebackground=True)
            y = y/F_cal_O2 #in nmol/s
            Y[mass] = np.trapz(y, x)
        y_hat = np.array([Y['M32'], Y['M34'], Y['M36']])
        y_hat = y_hat / np.sum(y_hat)
        
        # g is the H2(16)O / H2(18)O ratio, called gamma elsewhere
        def sqerror(g):
            return (g**2 - y_hat[0])**2 + (2*g*(1-g) - y_hat[1])**2 + ((1-g)**2 - y_hat[2])**2
        def testr(g):
            return np.array([g**2, 2*g*(1-g), (1-g)**2])
        res = minimize(sqerror, 0.5)
        g = res.x[0]
    
    tspan_CO2 = [-10, 100]
    x = np.arange(3, 60, 0.5)
    ys = []
    for i, mass in [(0,'M44'), (1,'M46'),(2,'M48')]:
        x0, y0 = get_signal(data, mass, tspan=tspan_CO2, unit='pA', removebackground=True)
        y0 = y0/F_cal_CO2 # in pmol/s
        f = interp1d(x0, y0, kind='linear')
        ys += [f(x)]
    ys = np.array(ys)
    ysum = np.sum(ys, axis=0)
    yhats = ys / ysum
    
    fig3, ax3 = plt.subplots()
    ax3.plot(x, yhats[0], standard_colors['M44'])
    ax3.plot(x, yhats[1], standard_colors['M46'])
    ax3.plot(x, yhats[2], standard_colors['M48'])
    #ax3.set_yscale('log')
    ax3.set_ylim([0,1])
    ax3.set_xlim([0, 60])
    ax3.set_ylabel('partial signal')
    ax3.set_xlabel('time / [s]')
    plt.savefig('COox_partial_signals.png')
    
    #-------------- solve the model ------------
    x_model = np.arange(0, 60, 0.5)
    
    literaturek = False
    if literaturek:
        SS = solve_carbonic_burst(k=0.026, g=g, tspan=x_model)
        ax3.plot(x_model, SS[:,0], ':', color=standard_colors['M44'])
        ax3.plot(x_model, SS[:,1], ':', color=standard_colors['M46'])                
        ax3.plot(x_model, SS[:,2], ':', color=standard_colors['M48']) 
        plt.savefig('COox_literature_k.png')
    
    fitk = True
    if fitk:
        SS = solve_carbonic_burst(k=0.08, g=g, tspan=x_model)
        ax3.plot(x_model, SS[:,0], '--', color=standard_colors['M44'])
        ax3.plot(x_model, SS[:,1], '--', color=standard_colors['M46'])                
        ax3.plot(x_model, SS[:,2], '--', color=standard_colors['M48'])    
        plt.savefig('COox_fit_k.png')        
        


isotopex = True
if isotopex:
    data = datasets['02']
    
    tspan_try1 = [900, 1800]
    ax1a = plot_experiment(data, mols=['He','CO',mols['M32']], tspan=tspan_try1, 
                           logplot=False, unit='pmol/s')
    ax1a += [ax1a[0].twinx()]
    plot_flux(data, mols=[mols['M44'], mols['M46'], mols['M48']], ax=ax1a[3], 
              tspan=tspan_try1, logplot=False, unit='pmol/s')
    ax1a[0].set_xlim(ax1a[1].get_xlim())
    ax1a[3].set_ylim([0, 5])
    plt.savefig('COox_try1_a')
    t_switch1 = 1584
    data1 = data.copy()
    for col in data1['data_cols']:
        if is_time(col):
            data1[col] = data1[col] - t_switch1
    ax1b = plot_experiment(data1, mols=[mols['M44'], mols['M46'], mols['M48']], 
                           tspan=[-100, 300], logplot=False, J_str='J / [uA/cm^2]')
    ax1b[1].set_ylim([0.6, 1.0])
    plt.savefig('COox_try1_b')


    tspan_try2 = [4100, 5200]
    ax2a = plot_experiment(data, mols=['He','CO',mols['M32']], tspan=tspan_try2, 
                           logplot=False, unit='pmol/s')
    ax2a += [ax2a[0].twinx()]
    plot_flux(data, mols=[mols['M44'], mols['M46'], mols['M48']], ax=ax2a[3], 
              tspan=tspan_try2, logplot=False, unit='pmol/s')
    ax2a[0].set_xlim(ax2a[1].get_xlim())
    ax2a[3].set_ylim([0, 5])
    plt.savefig('COox_try2_a')
    t_switch2 = 5043
    data2 = data.copy()
    for col in data2['data_cols']:
        if is_time(col):
            data2[col] = data2[col] - t_switch2
    ax2b = plot_experiment(data2, mols=[mols['M44'], mols['M46'], mols['M48']], 
                           tspan=[-100, 300], logplot=False, J_str='J / [uA/cm^2]')
    ax2b[1].set_ylim([0.6, 1.0])
    plt.savefig('COox_try2_b')


    tspan_try3 = [12000, 13000]
    ax3a = plot_experiment(data, mols=['He','CO',mols['M32']], tspan=tspan_try3, 
                           logplot=False, unit='pmol/s')
    ax3a += [ax3a[0].twinx()]
    plot_flux(data, mols=[mols['M44'], mols['M46'], mols['M48']], ax=ax3a[3], 
              tspan=tspan_try3, logplot=False, unit='pmol/s')
    ax3a[0].set_xlim(ax3a[1].get_xlim())
    ax3a[3].set_ylim([0, 5])
    plt.savefig('COox_try3_a')
    t_switch3 = 12825
    data3 = data.copy()
    for col in data3['data_cols']:
        if is_time(col):
            data3[col] = data3[col] - t_switch3
    ax3b = plot_experiment(data3, mols=[mols['M44'], mols['M46'], mols['M48']], 
                           tspan=[-100, 300], logplot=False, J_str='J / [uA/cm^2]')
    ax3b[1].set_ylim([0, 0.2])
    plt.savefig('COox_try3_b')


    #not bothering to change all the variable names again.
    tspan_try3 = [15500, 17000]
    ax3a = plot_experiment(data, mols=['He','CO',mols['M32']], tspan=tspan_try3, 
                           logplot=False, unit='pmol/s')
    ax3a += [ax3a[0].twinx()]
    plot_flux(data, mols=[mols['M44'], mols['M46'], mols['M48']], ax=ax3a[3], 
              tspan=tspan_try3, logplot=False, unit='pmol/s')
    ax3a[0].set_xlim(ax3a[1].get_xlim())
    ax3a[3].set_ylim([0, 5])
    plt.savefig('COox_try4_a')
    t_switch3 = 16825
    data3 = data.copy()
    for col in data3['data_cols']:
        if is_time(col):
            data3[col] = data3[col] - t_switch3
    ax3b = plot_experiment(data3, mols=[mols['M44'], mols['M46'], mols['M48']], 
                           tspan=[-100, 300], logplot=False, J_str='J / [uA/cm^2]')
    ax3b[1].set_ylim([0, 0.2])
    plt.savefig('COox_try4_b')

        