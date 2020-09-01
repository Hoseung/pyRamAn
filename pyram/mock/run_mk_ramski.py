import os
import numpy as np
import pickle
import pyram
from pyram.tree import halomodule as hmo


fn_template = './template_zubko.ski'
nout=925
out_dir = f'gal_{nout:05d}/'
if not os.path.isdir(out_dir): os.mkdir(out_dir)
#gal_cat = np.array([13, 6.854e+10, 0.48733889, 0.47836283, 0.49842872, 1.5930e-04])
#gcats = np.genfromtxt(f"./centers_{nout}.txt",
#                       dtype=[('id','int'), ('x','float'),('y','float'),('z','float')])

sim_dir = './'

gcat = hmo.Halo(base=sim_dir, nout=nout, is_gal=True, double=True)

#gal = gcat[50]
#radius = 20/(100*1e3/0.704) 20kpc fixed radius

with open(out_dir+f"centers_{nout}.txt", 'w') as fsave:
    for gdat in gcat.data[:50]:
        make_ram_input(sim_dir, nout, gcats, fsave,
                          dir_out=out_dir,
                          fn_template=fn_template,
                          ski_params=None,
                          smooth = 50,
                          nvec_frac=-1,
                          plot_stellar=False,
                          more_props=False,
                          r_fixed=None)


        #f.write("{} {} {} {} {}\n".format(gal['id'], gal['x'],gal['y'],gal['z'], radius))

