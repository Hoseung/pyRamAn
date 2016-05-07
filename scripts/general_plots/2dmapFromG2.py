# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 00:11:20 2015

@author: hoseung
"""

#from os import path
import pickle

#base = '/home/hoseung/Work/data'
#sim_dir = ['C29195', 'C01605_hi/', 'C04466/', 'G2/kisti'][3]
#base_dir = path.join(base, '', sim_dir, '')
base_dir = './'

#ptype = 'dm'
#ptype = 'star'
#nout_ini = 1
#nout_fi = 12

for nout in range(120, 140, 2):
    snout = str(nout).zfill(3)
    fin = base_dir + snout + 'map.pickle'
    try:
        with open(fin, 'rb') as f:
            img = pickle.load(f)

        ptimg = img.ptden2d
        fout = base_dir + snout + "dmmap_" + ptimg.proj + ".png"
        img.ptden2d.plot_2d_den(save=fout, show=False, vmin=1e7, vmax=3e10, dpi=400)
    #    fout = base_dir + snout + "gasmap_" + img.hydro.proj + ".png"
    #    img.hydro.plot_2d_den(save=fout, show=True)#, vmin=1e7, vmax=3e10, dpi=500)
    except:
        pass
    #%%
