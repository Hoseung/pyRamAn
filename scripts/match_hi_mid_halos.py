# -*- coding: utf-8 -*-
"""
Matching halos in different resolution simulation is too difficult... 

I will get back to it later.. 

Created on Sat Jun 20 14:06:58 2015

@author: hoseung
"""
import tree
import numpy as np
import load.info

wdir_m = '/home/hoseung/Work/data/C01605_mid/'
wdir_h = '/home/hoseung/Work/data/C01605_hi/'
fnml_m =  wdir_m + 'cosmo_150.nml'
fnml_h =  wdir_h + 'cosmo_370.nml'
nml_m = load.info.Nml(fnml_m)
nml_m.load()
nml_h = load.info.Nml(fnml_h)
nml_h.load()

#%%
def convertnout(nout1, nml1, nml2):   
    """
    Returns nout of sim2 corresponding to nout1 in sim1. 
    Note that nout starts from 1. axep[0] = first nout.
    """
    aexp = nml1.aout[nout1 -1]
    print(aexp)
    nout2 = np.where(nml2.aout == aexp)[0]+1
    return int(nout2)

nout_m = 130
nout_h = convertnout(nout_m, nml_m, nml_h)
print(nml_h.aout[nout_h - 1])
print(nout_h)
# So far so good.
# But, need odd nunmber outputs. 
# Lets get back to this later. 
# - 2015. 06. 20



#%%
hh = tree.halomodule.Halo(nout=convertnout(nout_m, nml_m, nml_h), halofinder="RS", base=wdir_h)
hh.load()
hm = tree.halomodule.Halo(nout=nout_m, halofinder="RS", base=wdir_m)
hm.load()
#%%
# First, crop larger array so that the two array cover roughly the same range of
# values. 
# This case, mass is importance.

ind_hh = np.where(hh.data['mvir'] > 5 * min(hm.data['mvir']))[0]
ind_hm = np.where(hm.data['mvir'] > 5 * min(hm.data['mvir']))[0]
print("original length:", len(hh.data.id))
hh2= tree.halomodule.Halo()
hh2.derive_from(hh, ind_hh)

hm2= tree.halomodule.Halo()
hm2.derive_from(hm, ind_hm)

print("Cropped :",len(hh2.data.id))

indmm = np.argsort(hm2.data['mvir'])
indmm = indmm[-2:]
indhh = np.argsort(hh2.data['mvir'])
indhh = indhh[-20:]

hhm= tree.halomodule.Halo()
hhm.derive_from(hm2, indmm)
hhh= tree.halomodule.Halo()
hhh.derive_from(hh2, indhh)
#%%
# Compensate for cluster offset



#%%
#hhm = hhh
import utils
utils.util.reimport(zt)
#mi = zt.match_diff_arr(np.log10(hhh.data['mvir']) * 1e-4, np.log10(hhm.data['mvir'])* 1e-4,
mi = zt.match_diff_arr(hhh.data['x'], hhm.data['x'],
                       hhh.data['y'], hhm.data['y'],
                       z1=hhh.data['z'], z2=hhm.data['z'],
                        tolerance=0.0095, window=0.9, reverse=True)

i_h = np.where(mi > -1)
i_m = mi[i_h]
print(len(i_m))
#%%
for i in range(len(i_m)):
    print(hhh.data['x'][i_h][i], hhm.data['x'][i_m][i])
print(" ")
for i in range(len(i_m)):
    print(hhh.data['y'][i_h][i], hhm.data['y'][i_m][i])
print(" ")
for i in range(len(i_m)):
    print(hhh.data['z'][i_h][i], hhm.data['z'][i_m][i])
print(" ")    
for i in range(len(i_m)):
    print(np.log10(hhh.data['mvir'][i_h][i]), np.log10(hhm.data['mvir'][i_m][i]))
#%%
import zoom.utils as zt
# mi = "matched" index of second data set
mi = zt.match_diff_arr(hh.data['x'], hm.data['x'],
                       hh.data['y'], hm.data['y'],
                       z1=hh.data['z'], z2=hm.data['z'])
#%%

                       
                    
#%%                       
i_h = np.where(mi > 0)
i_m = mi[i_h]
#%%
print(len(i_ok))
#print(hh.data.id[mi[mi > 0]])
import matplotlib.pylab as plt

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
#ax.scatter(hh.data['x'][i_h], hh.data['y'][i_h])
#h1 = ax.hist2d(hh.data['x'][i_h], hh.data['y'][i_h], bins=50)
h2 = ax.hist2d(hm.data['x'][i_m], hm.data['y'][i_m], bins=50)
#ax.scatter(hh.data['x'][i_m], hm.data['y'][i_m])

