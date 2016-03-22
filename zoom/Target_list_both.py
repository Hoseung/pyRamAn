# -*- coding: utf-8 -*-
"""
Created on Tue May 12 21:39:32 2015

@author: hoseung
"""
import pickle

def dump_list(data, fname):
    """outputs zoom-in target list.
        Assuming no normalization of units
    """
    ids = cat.data.id
    x = cat.data.x
    y = cat.data.y
    z = cat.data.z
    vx = cat.data.vx
    vy = cat.data.vy
    vz = cat.data.vz
    mvir = cat.data.mvir * 1e-13
    rvir = cat.data.rvir * cat.info.pboxsize
    nsub = cat.data.nsub

    with open(fname, 'w') as f:
        f.write("#      ID        x          y         z       vx      vy     vz   ")
        f.write("    Rvir(Mpc)   Mvir    #sub \n")
        for i in range(len(ids)):
            f.write("{:<4}   {:<4}   {:.5f}  {:.5f}  {:.5f}".format(i, ids[i],x[i],y[i],z[i]))
            f.write("  {:.3f}  {:.3f}  {:.3f}".format(vx[i],vy[i],vz[i]))
            f.write("   {:.3f}   {:.3f}   {:<3}     \n".format(rvir[i],mvir[i],nsub[i]))


f = open("/home/hoseung/Work/data/DMO_000/halo_sub.pickle", 'rb')
cat = pickle.load(f)
cat.data.x = cat.data.x - 0.5
cat.data.y = cat.data.y - 0.5
cat.data.z = cat.data.z - 0.5

f = open("/home/hoseung/Work/data/DMO/halo_sub.pickle", 'rb')
cat_org = pickle.load(f)
#%%
# Match two halos
import zoom.utils as zt
# mi = "matched" index of second data set
mi = zt.match_diff_arr(cat_org.data.x, cat.data.x, cat_org.data.y, cat.data.y, tolerance=500)
print(cat.data.id[mi])

#%%

dump_list(cat.data, './zoom_target_offset.txt')
dump_list(cat_org.data, './zoom_target_org.txt')            
