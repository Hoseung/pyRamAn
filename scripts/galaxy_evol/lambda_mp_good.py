# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 22:06:42 2015

@author: hoseung
"""

import numpy as np

a = np.zeros(10)

b = [0,1,4,7]

c= a[b]

print(c)
c[2] = 1.2

print(c)

print(a)

#%%

x = np.array([(1.5, 4), (1.0, 2), (3.0, 4)], dtype=[('x', float), ('y', int)])
ind = np.where(x['x'] < 2)
b = x[ind]


#%%
from tree import tmtree
import tree.halomodule as hmo 
import utils.sampling as smp
wdir =  '/home/hoseung/Work/data/01605/'
tt = tmtree.load(work_dir=wdir, filename="halo/TMtree.fits")


m_halo_min = 2e10
nout_fi = 187

hh = hmo.Halo(base=wdir, nout=nout_fi, halofinder='HM', info=info, load=True)
    
#halo = hmu.load_data(nout_fi, work_dir=work_dir, normalize=True)
i_center = np.where(hh.data['np'] == max(hh.data['np']))
i_satellites = smp.extract_halos_within(hh.data, i_center, scale=r_cluster_scale)
print("Total {0} halos \n{1} halos are selected".format(
      len(i_satellites),sum(i_satellites)))

# halos found inside the cluster and has tree back to nout_ini
large_enugh = hh.data['mvir'] > m_halo_min
halo_list = hh.data['id'][i_satellites * large_enugh]
    
h_ind_ok, halo_ok = tmtree.check_tree_complete(tt, 0, nout_fi - nout_ini0, halo_list)
print(len(halo_ok), "halos left")
final_gal = halo_ok[:,0]
ngals = len(final_gal)

#%%
import matplotlib.pyplot as plt
plt.plot(np.log10(hh.data['mvir'][large_enugh]))
plt.show()

#%%
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]

#%%
print(list(chunks(range(10),10)))

#%%

l = 21
n = 5

arr=[]
[arr.append([]) for i in range(5)]
for i in range(l):
    j = i % n
    arr[j].append(i)

#%%


def chunks(l, n):
    n = max(1, n)
    return [l[i:i + n] for i in range(0, len(l), n)]
    
print(chunks(np.arange(21),4))