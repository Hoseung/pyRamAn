
# coding: utf-8

# In[1]:
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from utils.catalogs import *
import tree.ctutils as ctu
import pickle
import tree.halomodule as hmo
import load


# In[2]:

def clustocentric_dist(cat, hal, info):
    cluster = hal.data[hal.data['np'].argmax()]
    cxc = cluster['x']*info.pboxsize
    cyc = cluster['y']*info.pboxsize
    czc = cluster['z']*info.pboxsize
    dist = np.sqrt(np.square(cxc - cat['xc'])
                 + np.square(cyc - cat['yc'])
                 + np.square(czc - cat['zc']))
    clu_rvir = cluster['rvir'] * info.pboxsize

    return clu_rvir, dist


# In[3]:
clusters = ['35663', '17891', '01605', '05427', '36413', '39990', '36415',
            '07206', '06098', '10002', '28928']


# In[7]:

# single nout.

cdir = 'catalog_GM/'
cmap = plt.cm.get_cmap('jet', 5)

verbose=True
nout = 67
ngals_tot = 0

cuts = [0,2,100]
lambdas = [[]] * 4


colors1 = ['red', 'orangered', 'orange']
colors2 = ['blue', 'royalblue', 'skyblue']
titles  = ['z = 0.0', 'z = 0.5', 'z = 1.0', 'z = 2.0']

rvir = []
fig, ax = plt.subplots(2,2, sharex=False, sharey=False)  # [187, 121, 87, 54, 37]
ax = ax.ravel()
for inout, nout in enumerate([187, 121, 87, 54]):
    cluster = clusters[0]
    wdir = './' + cluster + '/'
    # main galaxy list
    print(cluster)
    cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout) + '.pickle', 'rb'))
    info = load.info.Info(base=wdir, nout=nout, load=True)
    hal = hmo.Halo(base=wdir, nout=nout, is_gal=False, 
                   halofinder='HM', return_id=False,
                   load=True, info=info)
    rr, dist = clustocentric_dist(cat, hal, info)
    if nout == 187:
        rvir.append(rr)
    near = cat['lambda_r'][dist < rvir[0]]
    far = cat['lambda_r'][dist > rvir[0]]
#    print(near)
    for iclu in range(1,len(clusters)):
        cluster = clusters[iclu]
        print(cluster)
        wdir = './' + cluster + '/'
        # main galaxy list

        cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout) + '.pickle', 'rb'))
        info = load.info.Info(base=wdir, nout=nout, load=True)
        hal = hmo.Halo(base=wdir, nout=nout, is_gal=False, 
                       halofinder='HM', return_id=False,
                       load=True, info=info)
        aa, bb = clustocentric_dist(cat, hal, info)
        if nout == 187:
            rvir.append(aa)
        near = np.concatenate((near,cat['lambda_r'][bb < rvir[iclu]]))
        far = np.concatenate((far,cat['lambda_r'][bb > rvir[iclu]]))
        print(rvir[iclu], len(near), len(far))
    
    ax[3-inout].hist(near, alpha=0.9, color=colors1[0],
                     histtype='step', label='near',
                     lw=2, bins=12, range=[0.0, 0.6])
    ax[3-inout].hist(far, alpha=0.9, color=colors2[0],
                     histtype='step', label='far',
                     lw=2, bins=12, range=[0.0, 0.6])#, stacked=True)

    ax[3-inout].set_title(titles[inout])
    
plt.legend()
plt.savefig("environment_dependence.png")

