# -*- coding: utf-8 -*-
"""
dump halo list into a CSV file.
Created on Wed Jun 24 16:05:27 2015

@author: hoseung
"""
import tree
import load
import utils.sampling as smp
import numpy as np
work_dir = '/home/hoseung/Work/data/DMO/'
info = load.info.Info(nout=80, base=work_dir)
info.read_info()
h = tree.halomodule.Halo(base=work_dir, nout=80, info=info)
h.set_halofinder('HaloMaker')
h.load_hm()


m_threshold = 5e13
region = smp.set_region(centers=[0.5]*3, radius=0.5)

ind = np.where((h.data.x > region["xr"][0]) & (h.data.x < region["xr"][1]) &
                (h.data.y > region["yr"][0]) & (h.data.y < region["yr"][1]) &
                (h.data.z > region["zr"][0]) & (h.data.z < region["zr"][1]) &
                (h.data.mvir > m_threshold))[0]

with open(work_dir + 'halo_list.txt', 'w') as f:
    x = h.data['x'][ind]# * info.pboxsize
    y = h.data['y'][ind]# * info.pboxsize
    z = h.data['z'][ind]# * info.pboxsize
#    vx = h.data['vx'][ind]# * info.kms
#    vy = h.data['vy'][ind]# * info.kms
#    vz = h.data['vz'][ind]# * info.kms
    r = h.data['rvir'][ind] * info.pboxsize
    m = h.data['mvir'][ind]# * 1e11
    m2 = h.data['m'][ind]
    ids = [int(i) for i in h.data['id'][ind]]

    for i in range(len(ids)):
        f.write("{:<4}   {:.5f}  {:.5f}  {:.5f}".format(ids[i],x[i],y[i],z[i]))
#        f.write("  {:.3f}  {:.3f}  {:.3f}".format(vx[i],vy[i],vz[i]))
        f.write("  {:.6f}  {:.0f}  {:.0f}     \n".format(r[i],m[i],m2[i]))
        
  
#%%