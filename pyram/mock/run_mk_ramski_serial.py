import os
import numpy as np
import pickle
import pyram
from pyram.tree import halomodule as hmo
from pyram.mock import make_ramski
from pyram.utils import match

sim_dir = './'#'/storage6/NewHorizon/'
fn_template = './template_zubko.ski'
gid_fi = 1 # ID at nout = 876
gtree = np.load(f"./MJ_tree/morp_{gid_fi}.npy")

#gtree = gtree[3:]

# change fields
names = list(gtree.dtype.names)
# galid -> id
names[1] = 'id'
# : vrot -> r
names[6] = 'r' 
gtree.dtype.names = names

# we don't have #867 yet...
gid_fi = gtree['id'][0]
nout_fi = gtree['nout'][0]

out_dir = f'gal_{nout_fi:05d}_{gid_fi}/'
if not os.path.isdir(out_dir): os.mkdir(out_dir)

ski_params=None
inc=0
with open(out_dir+f"centers_{gid_fi}.txt", 'w') as fsave:
    for i, gt in enumerate(gtree):
        nout=gt['nout']
        s = pyram.load.sim.Sim(base=sim_dir, nout=nout)
        gcat = hmo.Halo(base=sim_dir, nout=nout, is_gal=True, double=True)
    
        i_this = np.where(gcat.data['id'] == gt['id'])
        gcats = gcat.data[i_this]#.reshape(-1,1)
        gcats['r'] = 2 * gt['R90']/s.info.boxtokpc# 


        #gcats=gt.reshape(-1,1) # as a iterable
        if i ==1:
            if inc > 180:
                inc-=180
                inc_ed = inc +90
            elif inc > 90:
                inc_ed = inc -90
            else:
                inc_ed = inc +90
            ski_params=[dict(name='face_on',  pixel_scale=50, nphoton=5e7, incli_off = inc,    azim_off = 0, roll_off = 0),
             dict(name='d30',      pixel_scale=50, nphoton=5e7, incli_off = inc+30*(inc<150)-30*(inc>150), azim_off = 0, roll_off = 0),
             dict(name='d60',      pixel_scale=50, nphoton=5e7, incli_off = inc+60*(inc<150)-60*(inc>150), azim_off = 0, roll_off = 0),
             dict(name='edge_on',  pixel_scale=50, nphoton=5e7, incli_off = inc_ed, azim_off = 0, roll_off = 0)]
    
        parameters = make_ramski.make_ram_input(sim_dir, nout, gcats, fsave,
                            dir_out=out_dir,
                            fn_template=fn_template,
                            ski_params=ski_params,
                            smooth = 50,
                            nvec_frac=0.3,
                            plot_stellar=True,
                            more_props=False,
                            r_fixed=None)
        if i==0:
            inc = parameters['inc']
            azim = parameters['azim']
            roll = parameters['roll']
