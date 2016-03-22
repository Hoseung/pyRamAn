# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 00:10:54 2015

@author: hoseung
"""
def all_tree_branches(tree, target_gal, nout_now,
                      r_flyby=100,
                      nout_ini_tree=None,
                      nout_fi_tree=None):
    """
    Searches for all fly-by trees around 'target_gal'.
    All trees that are
    
    parameters
    ----------
    tree :
        dd
    target_gal :
        ID of the target galaxy at the given nout "nout_now"
    nout_now :
        
    r_flyby : float [in kpc] (optional)
        All halos inside this radius at "nout_now" are tracked.
    
    """  
    import utils.sampling as smp
    # TreeMaker.
    if nout_ini_tree is None:
        nout_ini_tree = max(tree["NOUT"])
    if nout_fi_tree is None:
        nout_fi_tree = min(tree["NOUT"])
        
    i_flyby = smp.extract_halos_within(hh.data,
                                i_hal_now,
                                scale=0.001/hh.data['rvir'][i_hal_now],
                                Mcut=1e10)
                                
#%%
def get_idx(hnu, nout, tree):
    """
        Returns unique tree ID of a given halo.
    """
    i_now = np.where(tree["NOUT"] == nout)[0]
    return tree["IDX"][np.where(tree["HALNUM"][i_now] == hnu)[0]]

def full_tree(idx, tree):
    """
        Returns full main prg list of a galaxy.
    """
    idx_org = idx
    nouts = []
    idxs =[]
    while idx > 0 :
        nouts.append(tree[idx]["NOUT"])
        idxs.append(idx)  
        idx = tree[idx]["SONID"]
            
    # remove the first element 
    # because it will be added once more by the loop below.
    nouts.pop(0)
    idxs.pop(0)
    
    idx = idx_org
    while idx > 0:
        nouts.insert(0, tree[idx]["NOUT"])
        idxs.insert(0, idx)
        idx = tree["TREE"][idx][0]

    return nouts, idxs

#%%
def _is_new_tree(trees, this_tree):
    """
        True if it's a new tree.
        
        Parameters
        ----------
        trees: 
    """
    # check for length
    if trees == []:
        return True

    ll = [len(i) for i in trees]
#    print(len(this_tree))

    l_this = len(this_tree)
    same_length=[]
    for i, l in enumerate(ll):
        if l == l_this:
            same_length.append(i)

#    same_length = np.where(ll == len(this_tree))[0]
#    print(same_length, len(same_length))

    if len(same_length) == 0:
        return True
    #same_length = [same_length]
    import collections
    for i_tree in same_length:
#        print(i_tree)
# list[indices] does not work.
# So, there is no way iterate over a sub sample of a list
# except using a much simpler 'slicing'. 
# Instead, iterate over the indices.
        test_tree = trees[i_tree]
#        print(i_tree)
        if collections.Counter(test_tree) == collections.Counter(this_tree):
            return False
    return True

def compile_trees(trees, tree, idx_nearby, l_min=10):
    """
    
    """
    # Only NEW trees are appended.
    for idx_this in idx_nearby:
        tree_nout, this_tree = full_tree(idx_this, tree)
#        print(len(this_tree) > l_min)
        if (_is_new_tree(trees, this_tree) and len(this_tree) > l_min):
            trees.append(this_tree)
    

#%%
# normalize Tree XP first. 
# Halo position in physical size.
# Need pboxsize of all nouts...! 
# tree_c = tree.copy()
                                
# How about radius or mass??

#import tree.halomodule as hmo
import numpy as np
#import utils.util
from tree import tmtree
#import utils.sampling as smp
work_dir = '/home/hoseung/Work/data/05427/'
tree = tmtree.load(work_dir=work_dir, filename="halo/TMtree.fits")
pboxsizes =[]
nouts=[]
import load

#l_ind = 0
# Positions, virial radii, mass of all halos are needed. 
ll = len(tree["XP"][:,0])
x = np.zeros(ll)
y = np.zeros(ll)
z = np.zeros(ll)
rvir = np.zeros(ll)
mass = np.zeros(ll)

nout_fi = 187
nout_ini = 150

# Store needed values and normalize them to code unit (0 ~ 1)
# Normalize
for nout in range(7, nout_fi+1):
    info=load.info.Info(nout=nout, base=work_dir, load=True)
    pboxsizes.append(info.pboxsize)
    nouts.append(nout)
    ind_now = np.where(tree["NOUT"] == nout_fi - nout)[0]

    x[ind_now] = tree["XP"][ind_now][:,0] / info.pboxsize + 0.5
    y[ind_now] = tree["XP"][ind_now][:,1] / info.pboxsize + 0.5
    z[ind_now] = tree["XP"][ind_now][:,2] / info.pboxsize + 0.5
    rvir[ind_now] = tree["RVIR"][ind_now] / info.pboxsize + 0.5
    mass[ind_now] = tree["M"][ind_now]# / info.pboxsize + 0.5
# Table data는 변경이불가능한듯??
#%%    
# Nearby galaxies 
target_gal = 1618
#nout_now = 187

ind_prg, target_prg = tmtree.get_main_prg(tree, target_gal)
# ind_prg, target_prg are 2D lists.
#%%

xx=[]
yy=[]
zz=[]

trees=[]
for nout_now in range(nout_fi, nout_ini -1, -1):
    nout_tree = 187 - nout_now
    ind_now = np.where(tree["NOUT"] == nout_tree)[0]
    assert (len(ind_now) > 0), "No halos at the current nout"
    
    i_this_gal = np.where(tree["HALNUM"][ind_now] == target_prg[0][nout_tree])[0]

    xc = x[ind_now[i_this_gal]]
    yc = y[ind_now[i_this_gal]]
    zc = z[ind_now[i_this_gal]]

    xx.append(xc[0])
    yy.append(yc[0])
    zz.append(zc[0])

    idx_nearby = np.where((np.sqrt([(xc - x[ind_now])**2 +
                                   (yc - y[ind_now])**2 +
                                   (zc - z[ind_now])**2])[0] < 0.001) &
                                   (mass[ind_now] > 1e9))[0]
    # only branches with length > 10 are stored.
    compile_trees(trees, tree, idx_nearby, l_min=10)


#%%
# get x,y,z positions from idx.
x_trace=[]
y_trace=[]
z_trace=[]
for this_tree in trees:
    x_tmp = [x[i] for i in this_tree]
    y_tmp = [y[i] for i in this_tree]
    z_tmp = [z[i] for i in this_tree]
    x_trace.append(x_tmp)
    y_trace.append(y_tmp)
    z_trace.append(z_tmp)

#%%
# Plot trees

# Mark the main halo
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
from itertools import cycle
lines = ["-","--","-.",":"]
linecycler = cycle(lines)
for i in range(100):
    ax.plot(x_trace[i],y_trace[i], next(linecycler))

plt.show()

#%%
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range():
    ax.plot(x_trace[i],y_trace[i],z_trace[i])

ax.plot(xx,yy,zz, 'g:', linewidth=4)

# 3D 인데 2D인줄 아는듯..? 
plt.show()
