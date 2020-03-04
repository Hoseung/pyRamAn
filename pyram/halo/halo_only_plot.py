
# coding: utf-8

# In[40]:

wdir = '/home/hoseung/Work/data/05427/'
cluster_name = wdir.split('/')[-2]

frefine= 'refine_params.txt'
fnml = 'cosmo_200.nml'
nout_ini=37
nout_fi=187
nouts = range(nout_ini,nout_fi+1)

scale = 5
npix = 800

m_halo_min = 5e9

# data loading parameters
ptype=["star pos mass"]
refine_params = True
dmo=False
draw=True
draw_halos=True
draw_part = True
draw_hydro = False
if draw_hydro:
    lmax=19
    
import load
import utils.sampling as smp
import tree.halomodule as hmo
import utils.match as mtc
import draw
import pickle
import numpy as np
from tree import tmtree

info = load.info.Info(nout_fi, wdir, load=True)
tt = tmtree.load(work_dir=wdir, filename="GalaxyMaker/TMtree.fits")

# Galaxy ID list
hh = hmo.Halo(base=wdir, nout=nout_fi, halofinder='HM', info=info, load=True, is_gal=True)
i_center = np.where(hh.data['np'] == max(hh.data['np']))
i_satellites = smp.extract_halos_within(hh.data, i_center, scale=scale)
print("Total {0} halos \n{1} halos are selected".format(len(i_satellites),sum(i_satellites)))

large_enugh = hh.data['mvir'] > m_halo_min
halo_list = hh.data['id'][i_satellites * large_enugh]
h_ind_ok, gal_ok = tmtree.check_tree_complete(tt, 0, nout_fi - nout_ini, halo_list)
final_gal = gal_ok[:,0]
ngals = len(final_gal)
print("Ngal", final_gal)


# In[41]:

import matplotlib.pyplot as plt
import tree
from draw import pp
import os

for igal, gal in enumerate(final_gal[3:6]):
    gdir = wdir+"gal_"+str(gal).zfill(5) + '/'
    if not os.path.isdir(gdir):
        os.mkdir(gdir)
    for inout, nout in enumerate(nouts[::-1]):
        snout = str(nout).zfill(3)

        # Only selected galaxy
        hh = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info, load=True, is_gal=True)
        hind = np.where(hh.data['id'] == gal_ok[igal,inout])[0]
        h = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info)
        h.derive_from(hh, hind)
        print("Final target galaxy: ")
        print(" ID {} \n POS {} {} {} \n Rvir {}".format(h.data['id'],
                                                        h.data['x'], h.data['y'], h.data['z'],
                                                        h.data['rvir']*200000))

        region = smp.set_region(xc=h.data['x'],
                        yc=h.data['y'],
                        zc=h.data['z'],
                        radius = h.data['rvir'] * 5)
        s = load.sim.Sim(nout, wdir, dmo=dmo, ranges=region["ranges"], setup=True)
        imgs = draw.img_obj.MapSet(info=s.info, region=region)

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        ind = np.where(hh.data.mvir > 1e9)
        h_sub = tree.halomodule.Halo()
        h_sub.derive_from(hh, ind) 

        ax1.set_xlim(right=npix)
        ax1.set_ylim(top=npix)
        pp.pp_halo(h_sub, npix, region=region, axes=ax1, rscale=1, name=True)

        plt.savefig(gdir + snout + "_gals.png")
        plt.close()


# In[39]:

plt.close()


# In[ ]:



