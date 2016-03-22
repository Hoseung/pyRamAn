# -*- coding: utf-8 -*-
"""
Plots merger history

First used for SAMI project.
This is also useful in showing merging cluster's merger history.

Created on Thu Jun  4 22:58:03 2015
@author: hoseung
"""
# Load a tree 
from astropy.io import fits
from astropy.table import Table

import load
import tree
import numpy as np
import utils.sampling as smp
from utils import util
from tree import tmtree
import tree.draw_merger_tree as dmt
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

wdir = '/home/hoseung/Work/data/05427/'
#data = fits.getdata(wdir + "halo/TMtree.fits", 1)
data = fits.getdata(wdir + "GalaxyMaker/TMtree.fits", 1)
tt = Table(data)
tmtree.fix_nout(tt, 10, 187)

#%% 
nout_ini = 37
nout_fi = 187
nstep = nout_fi - nout_ini - 100

ttt = tt[tt["NOUT"] == 187] # final galaxies
for galid in ttt["HALNUM"][6:7]:
    fig = plt.figure(figsize=[6,6])
    plt.ioff()
    ax = fig.add_subplot(111)
    #galid = 6033
    sidgal = str(galid).zfill(5)
    
    nouts = np.unique(range(nstep + 2))
    zreds = np.unique(tt["Z"])[:nstep]
    zreds = ["%.2f" % i for i in zreds]
    #print(zreds)
    idx0 = dmt.get_idx(tt, galid, 187)
    dmt.recursive_tree(idx0, tt, nstep, ax, 0, 0, 1, mass_unit=1e9)
    
    # y axis label (redshift)
    ax.set_ylabel("Redshift")
    ax.set_xlim([-2.1,1.5])
    ax.set_ylim([-6, nstep + 2])
    plt.yticks(nouts[1:nstep:10], zreds[1:nstep:10])
    ax.set_title(sidgal)
    #fig.show()
    plt.savefig(wdir + "mergertrees/" + sidgal + '.png')
    plt.close()

