# -*- coding: utf-8 -*-
"""
Plots merger trees of (seemingly) merging clusters.
List of merging clusters is from visual inspection. 
TreeMaker gives not satisfying result. 
Try again with Consistent trees.

Created on Wed Jun 24 22:31:53 2015

@author: hoseung
"""
import load
from tree import halomodule
import utils.sampling as smp
import numpy as np
import tree.draw_merger_tree as dmt
import utils.match as mtc
import matplotlib.pyplot as plt

work_dir = '/home/hoseung/Work/data/DMO/'

m_threshold = 5e13  # in solar mass (physical).
nout_fi = 80
nout_ini= 10
info = load.info.Info(nout=nout_fi, base=work_dir)
info.read_info()

# tree
#%%
from astropy.io import fits
from astropy.table import Table
data = fits.getdata(work_dir + "/tree/tree_eccen_v2.fits", 1)
tt = Table(data)


#%%
# halo
hall = halomodule.Halo(base=work_dir, nout=nout_fi, info=info)
hall.set_halofinder('HaloMaker')
hall.load_hm()

# filter by id
ids = [2697, 4822, 9801, 17893, 25284, 29176, 32967,
       36413, 49368, 49919, 58214, 58803, 80181]
hind = mtc.match_list_ind(h.data.id, ids)

h = halomodule.Halo()
h.derive_from(hall, hind)

#%%
for cluster in h.data.id:
    print(cluster)
    fig = plt.figure(figsize=[6,6])
    plt.ioff()
    ax = fig.add_subplot(111)
    #galid = 6033
    scluster = str(cluster).zfill(5)
    
    nouts = np.unique(range(nout_fi - nout_ini + 2))
    zreds = np.unique(tt["Z"])[:len(nouts)]
    zreds = ["%.2f" % i for i in zreds]
    #print(zreds)
    idx0 = dmt.get_idx(tt, cluster, 0)
    dmt.recursive_tree(idx0, tt, 123, ax, 0, 0, 1)
    
    # y axis label (redshift)
    ax.set_ylabel("Redshift")
    ax.set_xlim([-1.1,2])
    ax.set_ylim([0,122])
    plt.yticks(nouts[1:122:10], zreds[1:122:10])
    ax.set_title(scluster)
    #fig.show()
    plt.savefig(work_dir + "mergertrees/" + scluster + '.png')
