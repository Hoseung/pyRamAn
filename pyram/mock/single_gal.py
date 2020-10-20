import os
import numpy as np
import pickle
import pyram
from pyram.tree import halomodule as hmo 
from pyram.mock import make_ramski
from pyram.utils import match
from pyram.utils.cosmology import Timeconvert as TC

def get_contam_fname(contam_dir, nout, purity="100"):
    cdir = os.path.join(contam_dir, purity+"PerCentPurity/")
    fname = cdir + f"list_gal_{nout:05d}.dat_nocontam"
    return fname

def get_contam(contam_dir, nout, gids):
    contam = np.genfromtxt(get_contam_fname(contam_dir, nout, purity="100"), 
                dtype=[("id", int), ("level", int), ("mass",float), 
                       ("x", float), ("y", float), ("z",float)])
    try:
        len(gids)
        inds = match.match_list_ind(contam['id'], gids, allow_swap=False) 
    except:
        inds = np.where(contam['id'] == gids)[0]
    return contam[inds]

def replace_name(arr, ind, name):
    temp = list(arr.dtype.names)
    temp[ind] = name
    arr.dtype.names = tuple(temp)


fn_template = './template_zubko.ski'
nout_fi=925
out_dir = f'gal_{nout_fi:05d}/'
if not os.path.isdir(out_dir): os.mkdir(out_dir)

sim_dir = './'
#gcat = hmo.Halo(base=sim_dir, nout=nout_fi, is_gal=True, double=True)

#pure = np.genfromtxt(f"contam/100PerCentPurity/list_gal_{nout:05d}.dat_nocontam")
#good_ids = pure[:,0].astype(int)

tree = np.load("MJ_tree/morp_1.npy" )
#gal = gcat[50]
#radius = 20/(100*1e3/0.704) 20kpc fixed radius
#ind_good = match.match_list_ind(gcat.data['id'], good_ids, allow_swap=False)
#gcats = gcat.data[ind_good[1:50]]#[:200]
#gcats['r'] *=10 # double the radius only for z=4 galaxies

r_multiply = 2 # FoV = 120% of R90 

contam_dir = './contam/'
galid_fi = tree['galid'][0]
gdats = np.zeros(len(tree), dtype=[('nout', int), ("id", int), ("level", int), ("mass",float),
                       ("r", float), ("x", float), ("y", float), ("z",float)])

for i, tt in enumerate(tree):
    nout = tt['nout']
    galid = tt['galid']
    contam = get_contam(contam_dir, nout, galid)
    for nn in contam.dtype.names:
        gdats[i][nn] = contam[nn]
    gdats[i]['r'] = tt['R90'] * r_multiply
    gdats[i]['nout'] = nout

pickle.dump(gdats, open(f"gtree_{galid_fi}_{nout_fi}.pickle", "wb"))

with open(out_dir+f"centers_{galid_fi}_{nout_fi}.txt", 'a+') as fsave:
    for gdat in gdats:
        print(gdat)
        nout = gdat['nout']
        # x,y,z in code unit, mass in Msun, r in kpc
        s=pyram.load.sim.Sim(nout, base=sim_dir, sim_type='nh')
        tc = TC(s.info)
        make_ramski.make_ram_input(s, tc, gdat, fsave,
                                       dir_out=out_dir,
                                       fn_template=fn_template,
                                       ski_params=None,
                                       smooth = 50, 
                                       nvec_frac=0.3,
                                       plot_stellar=True,
                                       more_props=True,
                                       r_fixed=None)
