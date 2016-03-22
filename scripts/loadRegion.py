# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 17:30:29 2015

@author: hoseung
"""


import load
import utils.sampling as smp

rr = load.runparam.RefineParam()
base_dir = "/home/hoseung/Work/data/C04466/"
rr.loadRegion(base_dir + 'refine_params.txt')

#print(rr.aexp)

nn = load.runparam.Nml()
nn.loadNml(base_dir + 'cosmo_200.nml')

nout = 100
aexp = nn.aout[nout]

from utils import match
i_aexp = match.closest(aexp, rr.aexp)

x_refine = rr.x_refine[i_aexp]
y_refine = rr.y_refine[i_aexp]
z_refine = rr.z_refine[i_aexp]
r_refine = rr.r_refine[i_aexp] * 0.5

region = smp.set_region(xc = x_refine, yc = y_refine, zc = z_refine, radius = r_refine)

#%%