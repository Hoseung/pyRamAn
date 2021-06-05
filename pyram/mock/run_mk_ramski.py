import os
import numpy as np
import pickle
import pyram
from pyram.tree import halomodule as hmo
from pyram.mock import make_ramski
from pyram.utils import match
from pyram.utils.cosmology import Timeconvert as TC

fn_template = './template_zubko.ski'
nouts =[870,871,873, 874, 876]#nouts = [371, 446, 599, 719, 925] #371 : z=1.3
nouts = [841,842,843,845,846,848,849,851,852]
#nouts=[599, 719]
for nout in [925]:
    out_dir = f'gal_{nout:05d}/'
    if not os.path.isdir(out_dir): os.mkdir(out_dir)

    sim_dir = './'
    s=pyram.load.sim.Sim(nout, base=sim_dir, sim_type='nh')

    #gcat = hmo.Halo(base=sim_dir, nout=nout, is_gal=True, double=True)
    pure = np.genfromtxt(f"contam/100PerCentPurity/list_gal_{nout:05d}.dat_nocontam", 
        dtype=[("id", int), ("level", int), ("m", float), ("x", float), ("y", float), ("z", float), ("r", float)])
    #good_ids = pure['id']

    #radius = 20 #20kpc fixed radius
    radius=None

    #ind_good = match.match_list_ind(gcat.data['id'], good_ids, allow_swap=False)
    #gcats = gcat.data[ind_good[1:10]]#[:200]
    gcats = pure[np.argsort(pure['m'])]
    #gcats['r'] *=10 # double the radius only for z=4 galaxies
    with open(out_dir+f"centers_{nout}_2.txt", 'w') as fsave:
        tc = TC(s.info)
        # radius to kpc
        gcats['r'] *= s.info.boxtokpc
        for gdat in gcats[::-1][:2000]:
            print("gdat[r']", gdat['r'])
            gdat['r'] = min((gdat['r'],50/s.info.h*100)) # kpc
            print("gdat[r']", gdat['r'])
            try:
                make_ramski.make_ram_input(s, tc, gdat, fsave,
                                           dir_out=out_dir,
                                           fn_template=fn_template,
                                           ski_params=None,
                                           smooth = 50, 
                                           nvec_frac=0.3,
                                           plot_stellar=True,
                                           more_props=False,
                                           r_fixed=radius)
            except:
                print("Failed, gid=",gdat['id'])
    print(f"nout {nout} done.")

