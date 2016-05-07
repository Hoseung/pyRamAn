# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:28:43 2015

@author: hoseung
"""

import numpy as np
import load
import time
#base = '/home/hoseung/Work/data/group3826/'
base = '/home/hoseung/Work/data/Lmax_136/'
#base = '/home/hoseung/Work/data/AGN2/'
nout=136
# add a method to automatically set ranges
import utils.sampling as smp

region = smp.set_region(xc=0.205, yc=0.78, zc=0.395, radius=0.03)
#print(region)
xrs = region["ranges"]
#xrs = [[0.21, 0.22], [0.77, 0.78], [0.39, 0.40]]
#%%
s = load.sim.Sim(nout, base, ranges=xrs)

#%%
s.add_amr()
s.add_hydro()

#%%
s.hydro.amr2cell(lmax=14)
s.hydro.cell.x