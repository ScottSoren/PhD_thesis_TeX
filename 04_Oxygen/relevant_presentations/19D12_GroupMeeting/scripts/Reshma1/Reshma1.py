#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 18:12:20 2018

@author: scott
"""

import os
import pickle
import numpy as np
from matplotlib import pyplot as plt
from EC_MS import trigger_cal, sync_metadata, cut_dataset, correct_ohmic_drop
from EC_MS import plot_experiment, Molecule, plot_vs_potential, make_selector
from EC_MS import activity_steps, calibration_curve, Chem, select_cycles
from EC_MS import get_current, get_potential, colorax


plt.close('all')

# constants
A_el = 0.196
data_directory = os.path.expanduser('~/Dropbox/DTU-MIT RuO2/Data/ECMS/')


with open('calibration.pkl', 'rb') as f:
    calibration = pickle.load(f)



measurements = {
                'Reshma1Ai':('18H28_Reshma1A_16O', '01', 'O16', 0.713, 500,
                             'Jazz4_CA', [1700, 5000],),
                'Reshma1Bi':('18H31_Reshma1B_18O', '05', 'O18', 0.730, 0,
                             'Emil3', [1300, 4600],),
                }


mols = []
for mass in ['M32','M34','M36']:
    mol = Molecule('O2')
    mol.primary = mass
    mol.name = mass
    mols += [mol]

He = Molecule('He')
mols += [He]
def T(M):
    return M**(-1/2)
ratio = He.get_RSF(transmission_function=T) / mols[0].get_RSF(transmission_function=T)

if True: # load data, make old plots
    datasets = {}
    for name, (folder, tag, isotope, RE_vs_RHE, R, cal, tspan) in measurements.items():

        F_cal = calibration[cal]['F_cal']
        for mol in mols:
            mol.F_cal = F_cal
        He.F_cal = F_cal * ratio

        F = data_directory + folder
        pkl_file_name = [f for f in  os.listdir(F) if f[0:2]==tag and f[-4:]=='.pkl'][0]
        with open(F + os.sep + pkl_file_name, 'rb') as pkl:
            dataset = pickle.load(pkl)
        V_str, J_str = sync_metadata(dataset, RE_vs_RHE=RE_vs_RHE, A_el=A_el)
        trigger_cal(dataset)
        V_str = correct_ohmic_drop(dataset, R_ohm=R)
        make_selector(dataset)
        datasets[name] = dataset

        data = cut_dataset(dataset, tspan=tspan, t_zero='start')
        ax = plot_experiment(data, mols=mols, removebackground=False)
       #ax[0].set_ylim([1e-1, 3e2])
        plt.savefig('fig2_' + name + '.png')

    if False: # zoom-in at 90 mV overpotential point
        x, y = mols[2].get_flux(datasets['Reshma1Bi'], tspan=[2940, 3130], unit='pmol/s')
        fig, ax = plt.subplots()
        ax.plot(x-2965, y-np.mean(y[0:5]), 'g')
        ax.set_xlabel('time / s')
        ax.set_ylabel('flux / [pmol/s]')
        fig.set_figwidth(2)
        fig.set_figheight(1.5)
        fig.savefig('fig2b_inset')

if True: # additional plots for Reshma1
    dataset = datasets['Reshma1Ai']

    plot_experiment(dataset, J_str='selector')

    cycles = [28,
              29, 30, 31, 32, 33, 34, 35, 36, 37, 38
              ]


    O2 = calibration_curve(dataset, mol='O2', cycles=cycles,
                             t_bg=[6100, 6200], n_el=4,
                             find_max=True, t_max_buffer=5, V_max_buffer=10,
                             mode='average', t_int=10,
                           )
    O2.name = 'M32'
    results = activity_steps(dataset, mols=O2, cycles=cycles,
                             t_bg=[6100, 6200], unit='pmol/s',
                             find_max=True, t_max_buffer=5, V_max_buffer=10,
                             mode='average', t_int=30)
    for ax in results['ax']:
        ax.set_xlim([3300, 5500])
    plt.savefig('quant_experiment.png')

    I = results['Qs'] * 1e6
    n_I = results['Qs'] / (4*Chem.Far) * 1e12
    n = results['ns']['M32']

    fig, ax = plt.subplots()

    ax.plot(I, n, '.', color='k', markersize=10)
    ax.plot(np.append(I, 0), np.append(n_I, 0), '--', color='r')
    ax.plot(0, 0, 'r.')

    ax.set_xlabel('final current / [$\mu$A]', fontsize=12)
    ax.set_ylabel('final cal. signal / [pmol/s]', fontsize=12)

    plt.savefig('activity_steps.svg')

    c = select_cycles(dataset, cycles=[33])

    fig, ax = plt.subplots()
    ax2 = ax.twinx()

    x, y = O2.get_flux(c, unit='pmol/s', background='constant')
    I_x = y * 1e-12 * 4 * Chem.Far
    t, I = get_current(c, unit='A')
    t_v, V = get_potential(c)

    t = t - 1700 # fixing t_zero
    t_v = t_v - 1700 # fixing t_zero
    x = x - 1700 # fixing t_zero
    ax.plot(x, I_x*1e6, 'k-')
    ax.plot(t, I*1e6, 'r--')
    ax2.plot(t_v, V, 'b-')

    ax.set_xlabel('time/s', fontsize=12)
    ax.set_ylabel('(partial) current / [$\mu$A]', fontsize=12)
    ax2.set_ylabel(V_str[:-1], fontsize=12)
    ax2.set_ylim([1.1, 2.2])
    ax2.set_yticks([1.2, 1.3, 1.4])
    ax.set_ylim([-12, 10])
    colorax(ax2, 'b')

    plt.savefig('partial_current.svg')




