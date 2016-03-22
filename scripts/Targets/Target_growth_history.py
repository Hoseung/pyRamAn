# -*- coding: utf-8 -*-
"""
Mass growth history of the target clusters

Created on Thu Jun 25 15:43:16 2015

@author: hoseung
"""

import load
import tree
import utils.match as mtc
#%%
work_dir = '/home/hoseung/Work/data/DMO/'
nout=80

ids = [32967,9782,5420,49919,36413,35663,32675,29176,24962,17891,10002,5427,
       4822,74010,49096,46814,41608,29195,29172,28928,28927,24954,14172,6098,
       49464,49367,41092,36370,36297,35861,35582,35401,28432,18212,18206,15233,
       5743,14168,28930,11402,26649,17545] # 43 target clusters

#%%
# load halo
hrs_all = tree.halomodule.Halo(nout=nout, halofinder="RS", base=wdir)
htm_all = tree.halomodule.Halo(nout=nout, halofinder="TM", base=wdir)
hrs_all.load()
htm_all.load()

ind_hrs = np.where(hrs_all.data['mvir'] > 3e13)[0]
ind_htm = np.where(htm_all.data['mvir'] > 3e13)[0]

hrs= tree.halomodule.Halo()
hrs.derive_from(hrs_all, ind_hrs)

htm= tree.halomodule.Halo()
htm.derive_from(htm_all, ind_htm)

#%%
def close_halos(h1, h2, id1, tol=0.000005):
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


for tm_id in ids:
    ind = close_halos(htm, hrs, tm_id) # 11402 = 39990
    print(hrs.data['id'][ind], hrs.data['x'][ind], hrs.data['y'][ind],
          hrs.data['z'][ind], hrs.data['mvir'][ind])
          
# What if there are more than one candidates?           
          
          
          
          
          
#%%
# load DMO 
#s = load.sim.Sim(nout=nout, base=work_dir, dmo=True)
#s.add_part(['dm id pos mass'])
#s.set_ranges() # whole volume
#s.part.load()
#s.part.dm['m'] *= s.info.msun
