# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 03:20:07 2015

@author: hoseung
"""

import matplotlib.pyplot as plt
import numpy as np
import load
from galaxymodule import galaxy

rscale = 0.8
npix=800
wdir = './'

#frefine= 'refine_params.txt'
#fnml = 'cosmo_200.nml'

ptypes=["star id pos mass vel"]

npix = 20
dtype_galevol = [('id', int), ('mtot', '<f8'), ('mgas', '<f8'),
                 ('mstar', '<f8'), ('mdm', '<f8'), ('mhal', '<f8'),
                 ('nstar', int), ('ndm', int), ('nsink', int),
                 ('pos', '<f8', 3), ('vel', '<f8', 3), ('lx', '<f8', 3),
                 ('dcluster', '<f8'), ('b2t', float), ('mag', float, 5),
                 ('sfr', '<f8'), ('lambdar', '<f8'), ('lambdar_arr', '<f8', npix), 
                 ("nfit",float), ("morph_vi",str),
                 ("morph_b2t",str), ("add2",float), ("add3",float)]

#%%
# Tree
import pickle
import tree
wdir = '/home/hoseung/Work/data/AGN2/'
f_tree = wdir + "rhalo/tree.pickle"
with open(f_tree, "rb") as ft:
    rstree = pickle.load(ft)

# normalize temporarily. 
rstree['x'] *= 200/199.632011
rstree['y'] *= 200/199.632011
rstree['z'] *= 200/199.632011
rstree['nout'] += 51
# normalize RS halo included in halomodule.

#%%
nouts=np.arange(30, 133)

import utils.sampling as smp
from utils import util
import draw
#%%
#util.reimport(tree)
# Select target galaxies
nout_final = 132
#i_final = tree['nout'] == nout_final
tt = rstree[rstree['nout'] == nout_final]


#%%
print(sum(tt['mvir'] > 1e14))
print(sum(tt['mvir'] > 1e13))

#idgals = [5232, 5495, 5543, 6061, 5491, 6191]

wdir = '/home/hoseung/Work/data/AGN2/'
info = load.info.Info(nout=132, base=wdir)
hall = tree.halomodule.Halo(nout=132, base=wdir, halofinder="RS", info=info)

hall.load()
# convert to code unit. - done by default
#hall.normalize()

# subset of halos ONLY inside zoom-in region
i_center = np.where(hall.data['np'] == max(hall.data['np']))[0]
h_ind = smp.extract_halos_within(hall, i_center, scale=1.0)
h = tree.halomodule.Halo()
#h.derive_from(hall, h_ind[5:35])
h.derive_from(hall, h_ind)#,  [4921, 5281, 5343, 5365, 5375, 5412, 5415], 5442, 5639, 5665, 6095])
import matplotlib.pyplot as plt
plt.close(0)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(np.log10(tt['mvir']), range=[9,15])
ax.hist(np.log10(h.data.mvir), range=[9,15])

idgals = [1932, 1955, 1969]

#%%

for idgal in idgals:
#    prg_treex, prg_tree = tree.tmtree.get_main_prg(tt, idgal, nout_ini=122, nout_fi=0)
    
    gal_evol = np.zeros(len(nouts), dtype=dtype_galevol)
    sidgal = str(idgal)

#    for inout, nout in enumerate(nouts[::-1]):
    for inout, nout in enumerate(range(132, 133)):
        print(" ID: {}, nout: {}".format(sidgal, nout))
        print("\n")
        if nout > 132:
            continue
        idgal_now = prg_tree[0][inout]
        snout = str(nout).zfill(3)
        # Load halo
    
        info = load.info.Info(nout=nout, base=wdir)

        hall = tree.halomodule.Halo(nout=nout, base=wdir, halofinder="HM", info=info)
        hall.load()

        i_center = np.where(hall.data['np'] == max(hall.data['np']))[0]
        x_clu = hall.data['x'][i_center]
        y_clu = hall.data['y'][i_center]
        z_clu = hall.data['z'][i_center]
        # single halo 
        h = tree.halomodule.Halo()

        h.derive_from(hall, np.where(idgal_now == hall.data.id)[0])
        print("Halo property:")
        print( h.data['x'], h.data['y'],h.data['z'], h.data.id, '\n')
        
        # Load simulation data
        region = smp.set_region(xc=h.data.x, yc=h.data.y, zc=h.data.z, radius = 0.002)#400kpc
        print(region['centers'])
        s = load.sim.Sim(nout, wdir)
        
        s.set_ranges(region["ranges"])
        s.show_cpus()
        
        s.add_part(ptypes)
        s.part.load()
        # convert to km/s
        s.part.star['vx'] *= s.info.kms
        s.part.star['vy'] *= s.info.kms
        s.part.star['vz'] *= s.info.kms
        s.part.star["m"]  *= s.info.msun

#%%
        xc_tmp = h.data['x']
        yc_tmp = h.data['y']
        zc_tmp = h.data['z']
        rr_tmp = h.data['rvir']
        
        x = s.part.star['x']
        y = s.part.star['y'] 
        z = s.part.star['z'] 
        m = s.part.star['m']

        ind = np.where((x - xc_tmp)**2 + (y - yc_tmp)**2 
                            + (z - zc_tmp)**2 < rr_tmp**2)[0]


        data = draw.pp.den2d(x[ind], y[ind], z[ind], m[ind], 400, s.info, region=None, cic=True)

        import utils.prettyplot as ptt
        
        vmax = 2e18
        log = True
        
        vmax = data.max()
        vmin = vmax * 1e-6
        
        fig = plt.figure()
        axes = fig.add_subplot(111)
        
        if log:
            from matplotlib.colors import LogNorm
            p = axes.imshow(data, cmap=plt.get_cmap('brg'),
                          norm=LogNorm(vmin=vmin, vmax=vmax))
        else:
            p = axes.imshow(data, cmap=plt.get_cmap('brg'))
        
        fig.savefig(wdir + "tmp.png")
#%%        
        util.reimport(galaxy)
        util.reimport(draw.pp)
        gal = galaxy.Galaxy(h.data, rscale=0.5, radius_method='simple', info=s.info)
        gal.mk_gal(part=s.part, cell=None, rscale=0.6) # First guess, stars inside 0.5Rvir.
        gal.cal_lambda_r(npix=npix, method=1, rscale=1.0) # lambdar_r out to 1 * Reff
#        gal_evol[inout]['lambdar_arr'][0:len(gal.lambda_r)] = gal.lambda_r# array length?
        gal.plot_gal(base=wdir + snout + 'prg_' + sidgal)

#%%        
        gal_evol[inout]['mstar'] = gal.mstar
        gal_evol[inout]['id'] = gal.id
        gal_evol[inout]['pos'][0] = gal.xc
        gal_evol[inout]['pos'][1] = gal.yc
        gal_evol[inout]['pos'][2] = gal.zc
        gal_evol[inout]['lambdar_arr'][0:len(gal.lambda_r)] = gal.lambda_r
        gal_evol[inout]['lambdar'] = 0.5*(gal.lambda_r[0.5*len(gal.lambda_r) -1] 
                                        + gal.lambda_r[len(gal.lambda_r) -1])
        gal_evol[inout]['dcluster'] = np.sqrt((gal.xc - x_clu)**2 + 
                                        (gal.yc - y_clu)**2 + 
                                        (gal.zc - z_clu)**2)/hall.data['rvir'][i_center]

#%%
#data = gal_evol['lambdar_arr']
#fig= plt.figure()
"""
def animate(nframe):
    plt.cla()
    plt.plot(data[nframe])
    plt.ylim(0,1)
    plt.title('Exp: %.2f'%(1.*nframe/10.))
"""
#anim = animation.FuncAnimation(fig, animate, frames=len(nouts))
# You need ffmpeg, which is not included in Ubuntu standard repository.
# add a PPA and then install ffmpeg.
#anim.save('animation.mp4', fps=30, 
#          extra_args=['-vcodec', 'h264', 
#                      '-pix_fmt', 'yuv420p'])