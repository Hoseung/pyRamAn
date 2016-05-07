# -*- coding: utf-8 -*-
"""
Match halos in two different simulations

Created on Thu Jun 25 19:32:19 2015

@author: hoseung
"""


import tree
import zoom.utils as zt



wdir = '/home/hoseung/Work/data/DMO/'
nout = 80

hrs_all = tree.halomodule.Halo(nout=nout, halofinder="RS", base=wdir)
htm_all = tree.halomodule.Halo(nout=nout, halofinder="TM", base=wdir)
hrs_all.load()
htm_all.load()
#%% Histogram first.
import matplotlib.pylab as plt

fig = plt.figure()
ax = fig.add_subplot(121)
h1 = ax.hist2d(hrs.data['x'], hrs.data['y'], bins=100)
ax = fig.add_subplot(122)
h2 = ax.hist2d(htm.data['x'], htm.data['y'], bins=100)
# Looks reasonably good.


#%% only large enough halos

ind_hrs = np.where(hrs_all.data['mvir'] > 3e13)[0]
ind_htm = np.where(htm_all.data['mvir'] > 3e13)[0]

hrs= tree.halomodule.Halo()
hrs.derive_from(hrs_all, ind_hrs)

htm= tree.halomodule.Halo()
htm.derive_from(htm_all, ind_htm)

print(len(htm.data.id))

#%%
i1 = np.argmax(hrs.data.mvir)
i2 = np.argmax(htm.data.mvir)
print(hrs.data['x'][i1], hrs.data['y'][i1], hrs.data['z'][i1])
print(htm.data['x'][i2], htm.data['y'][i2], htm.data['z'][i2])

#%%
import utils
utils.util.reimport(zt)
mi = zt.match_diff_arr(hrs.data['x'], htm.data['x'],
                       hrs.data['y'], htm.data['y'],
                        tolerance = 0.0005)#,
#                       z1=hrs.data['z'], z2=htm.data['z'],
#                        tolerance = 0.0001)
print("Matched halos: {} out of {}".format(sum(mi > -1), len(hrs.data.id)))
#%%
import numpy as np
i_r = np.where(mi > -1)
i_t = mi[i_r][0:100]
for i in range(len(i_t)):
    print(hrs.data['x'][i_r][i], htm.data['x'][i_t][i])
print(" ")
for i in range(len(i_t)):
    print(hrs.data['y'][i_r][i], htm.data['y'][i_t][i])
print(" ")
for i in range(len(i_t)):
    print(hrs.data['z'][i_r][i], htm.data['z'][i_t][i])
print(" ")    
for i in range(len(i_t)):
    print(np.log10(hrs.data['mvir'][i_r][i]), np.log10(htm.data['mvir'][i_t][i]))
#%%
for i in range(len(i_t)):
    print(hrs.data['id'][i_r][i], htm.data['id'][i_t][i])
#%%
# Manual labor..

def close_halos(h1, h2, id1, tol=0.0001):
    """
    returns indices of the counterpart candidates.
    """
    import numpy as np
    ind = np.where(h1.data['id'] == id1)[0]
    
    x = h1.data['x'][ind]
    y = h1.data['y'][ind]
    z = h1.data['z'][ind]

    xx = h2.data['x']
    yy = h2.data['y']
    zz = h2.data['z']

    return np.where(np.sqrt((x-xx)**2 + (y-yy)**2 + (z-zz)**2 < tol))[0]

ids = [32967,9782,5420,49919,36413,35663,32675,29176,24962,17891,10002,5427,
       4822,74010,49096,46814,41608,29195,29172,28928,28927,24954,14172,6098,
       49464,49367,41092,36370,36297,35861,35582,35401,28432,18212,18206,15233,
       5743,14168,28930,11402,26649,17545] # 43 target clusters
for tm_id in ids:
    ind = close_halos(htm, hrs, tm_id) # 11402 = 39990
    print(hrs.data['id'][ind], hrs.data['x'][ind], hrs.data['y'][ind],
          hrs.data['z'][ind], hrs.data['mvir'][ind])
    