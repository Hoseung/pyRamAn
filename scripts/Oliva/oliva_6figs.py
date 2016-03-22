
# coding: utf-8

# In[1]:

import pickle
import numpy as np
import load

import matplotlib.pyplot as plt
import tree.halomodule as hmo
import utils.match as mtc


# In[8]:

def bcg_prop(cat, verbose=False):
    sort_mstar = cat['mstar'].argsort()
    bcg = cat[-1]
    sbcg = cat[-2]
    dominance = bcg.mstar - sbcg.mstar
    if verbose:
        pass
        #print(bcg.mstar, sbcg.mstar)
    
    return bcg, dominance


def assign_recarr(recarr, ind, data, drop=None):
    """
        If dtype varies among clusters (different version of lambda_gal),
        available 
    """
    names = list(data.dtype.names)
    if drop is not None:
        for dd in drop:
            names.remove(dd)
    for fieldname in names:
        recarr[fieldname][ind] = data[fieldname]
        
        
def halo_of_gal(halo, catalog, galid, dim6=False):
    gal = catalog[catalog.id == galid]
    
    center = np.array([gal.xc, gal.yc, gal.zc, gal.rhalo, gal.vx, gal.vy, gal.vz])
    
    if dim6 :
        norm = np.sqrt(np.square(center[0] - halo.x) + 
                   np.square(center[1] - halo.y) + 
                   np.square(center[2] - halo.z) +
                   np.square(center[3] - halo.rvir) +
                   np.square(center[4] - halo.vx) + 
                   np.square(center[5] - halo.vy) + 
                   np.square(center[6] - halo.vz))
    else:
        norm = np.sqrt(np.square(center[0] - halo.x) + 
                   np.square(center[1] - halo.y) + 
                   np.square(center[2] - halo.z) +
                   np.square(center[3] - halo.rvir))

    i_match = norm.argmin()
    
    return halo[i_match]


# In[3]:

nout = 187
clusters = ['05427', '01605', '29172', '28928']
cdir = 'catalog_GM/'


# In[10]:

# check if the clusters have relevant data
check_file=False
if check_file:
    from glob import glob
    for i, cluster in enumerate(clusters):
        wdir = '/home/hoseung/Work/data/' + cluster + '/'
        cat_list = glob("")
        for file in glob(wdir + cdir + 'catalog' + str(nout) + '.pickle'):
            print(file)
        for file in glob(wdir + 'halo/DM/tree_bricks' + str(nout)):
            print(file)


# In[10]:

bcgs = np.zeros(len(clusters), 
                dtype=[('index', '<i8'), ('boxtokpc', '<f8'), ('id', '<i8'),
                       ('idx', '<i8'), ('lambda_r', '<f8'),
                       ('mgas', '<f8'), ('mstar', '<f8'), ('nstar', '<i8'),
                       ('rgal', '<f8'), ('rhalo', '<f8'),
                       ('vx', '<f8'), ('vy', '<f8'), ('vz', '<f8'),
                       ('xc', '<f8'), ('yc', '<f8'), ('zc', '<f8'),
                       ('mhalo', '<f8'), ('dominance', '<f8'), ('cluster', '<i8')])
# lambda_arr is removed from bcgs and also will be ignored in assign_recarr
# by drop=['lambda_arr'] option.

dominance =[]
rgal = []
mhalo = []
rhalo =[]
lambdar = []
mstar = []

for i, cluster in enumerate(clusters):
    wdir = '/home/hoseung/Work/data/' + cluster + '/' #'05427/'
    
    cat = pickle.load(open(wdir + cdir + 'catalog' + str(nout) + '.pickle', 'rb'))
    bcg, dom = bcg_prop(cat, verbose=True)
#    plot_gal_merger_history(cluster, bcg)
    # exclude lambda_arr 
    
    assign_recarr(bcgs, i, bcg, drop=['lambda_arr'])
    
    info = load.info.Info(base = wdir, nout=nout, load=True)
    hh = hmo.Halo(nout=nout, base=wdir, info=info, halofinder="HM", is_gal=False, load=True)
    hh.data.x *= info.pboxsize
    hh.data.y *= info.pboxsize
    hh.data.z *= info.pboxsize
    hh.data.rvir *= info.pboxsize
    
    halo = halo_of_gal(hh.data, cat, bcg.id) # halo data
    halo = hh.data[hh.data.np.argmax()]
    print("{:.4f} {:.4f} {:.4f} {:.4e}".format(halo['x'], halo['y'], halo['z'], halo.mvir))
    print("{:.4f} {:.4f} {:.4f} {:.4e}".format(bcg['xc'], bcg['yc'], bcg['zc'], bcg.mstar))
    #print(halo.mvir, hh.data.mvir[hh.data.np.argmax()])
    rgal.append(np.log10(bcg['rgal'])) # in kpc
    rhalo.append(np.log10(bcg['rhalo'] * info.pboxsize * 1000)) # in kpc  (/h?)
    mhalo.append(np.log10(halo['mvir']))
    lambdar.append(bcg['lambda_r'])
    mstar.append(np.log10(bcg['mstar']))
    dominance.append(dom)
    bcgs[i]['mhalo'] = halo['mvir']
    bcgs[i]['dominance'] = dom
    bcgs[i]['cluster'] = cluster

    


# In[17]:

#np.savetxt("ss.txt", bcgs)
np.save("Oliva_data.npy", bcgs)


# # confirm halos matching
# fig, axs = plt.subplots(2)
# axs[0].plot(halos.id, cat.id)
# axs[0].set_title("id vs id")
# axs[1].plot(halos.rvir, cat.rhalo)
# axs[1].set_title("rvir vs rvir")
# 
# plt.show()

# In[16]:

#samples = bcgs

fig, axs = plt.subplots(3,3)
axs = axs.ravel()

#rgal = np.log10(samples['rgal']) # in kpc
#rhalo = np.log10(samples['rhalo'] * info.pboxsize * 1000) # in kpc  (/h?)
#mhalo = np.log10(samples['mvir'])
#lambdar = samples['lambda_r']
#mstar = np.log10(samples['mstar'])

axs[0].scatter(mstar, lambdar, c = lambdar)
axs[0].set_title("rotation vs Mstar, fig3")

axs[1].scatter(rgal, lambdar, c = lambdar)
axs[1].set_title("rotation vs Rgal, fig7")

axs[2].scatter(mhalo, lambdar, c = lambdar)
axs[2].set_title("rotation vs Mhalo, fig8")

axs[3].scatter(dominance, lambdar, c = lambdar)
axs[3].set_title("rotation vs dominance, fig9")

axs[4].scatter(mstar, rgal, c = lambdar)
axs[4].set_title("Rgal vs Mstar, fig10")

axs[5].scatter(mhalo, rhalo, c = lambdar)
axs[5].set_title("Mhalo vs Mstar, fig11")

plt.suptitle("nout = {}, z= {:.3f}".format(str(nout), info.zred))
plt.tight_layout()
#plt.show()
plt.savefig('Oliva_fig.png', dpi=200)


# In[ ]:



