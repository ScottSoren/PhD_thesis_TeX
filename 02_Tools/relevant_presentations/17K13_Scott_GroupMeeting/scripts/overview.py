#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 21:14:26 2017

@author: scott@fysik.dtu.dk
EC_MS package available at https://github.com/ScottSoren/EC_MS
"""

'''
Ahhhh.... that really happened. And now my life is about getting back to 
Berkeley. I know I can go back, probably go back often, maybe even go back soon. 
But I need to focus. Focusing on this is a prerequisite for getting back to the 
reality of the girl who I'd like to daydream about now. Here we go.
'''

'''
experimental NOTES:
    17J03, https://cinfelog.fysik.dtu.dk/EC-MS/47
    17J04, https://cinfelog.fysik.dtu.dk/EC-MS/48
    17J05, https://cinfelog.fysik.dtu.dk/EC-MS/49
'''

import os
import numpy as np
from matplotlib import pyplot as plt
from EC_MS import import_folder, sync_metadata, plot_experiment

plt.close('all')

importrawdata = True
if importrawdata:
    data_directory = os.path.expanduser('~/Dropbox/Sniffer_Experiments/02_NiFeNPs/Data')
    folders = {'17J03':'17J03_Pt_calibration',
               '17J04':'17J04_Pt_isotope_exchange',
               '17J05':'17J05_Pt_isotope_exchange'
               }
    
    tagtype = None # 'all' to import separated by tags
    all_data = {}
    for key, f in folders.items():
        all_data[key] = import_folder(data_directory + os.sep + f, tags=tagtype)

plotrawdata = True
if plotrawdata:
    for key, data in all_data.items():
        sync_metadata(data, RE_vs_RHE=0.747, A_el=0.196)
        plot_experiment(data, masses=['M2','M4','M28','M32','M34','M36','M44','M46','M48'])
        plt.savefig('overview_' + key + '.png')


        
    
    
    

