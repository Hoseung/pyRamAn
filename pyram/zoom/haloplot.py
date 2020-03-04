# -*- coding: utf-8 -*-
"""
1) Load halo_sub, which is the complete halo list avobe a mass cut. 
2) Draw density map of halos and print halo ID.
3) Plots cover 20Mpc cube so that I can find merging halos.
--------
To do :
For selected halos, check merger history (next step.)


Created on Fri May 15 01:32:03 2015

@author: hoseung
"""

import load
import utils.sampling as smp
from draw import pp
import zoom
import pickle

work_dir = '/home/hoseung/Work/data/DMO/'
#work_dir = './'

rhaloscale = 3.0 
npix = 800
show = False
denmap = False

nout_fi = 80
nout_ini = 1
nout_list = [i for i in range(0, nout_fi + 1, 5)]
nout_list[0] = 1

# nout_list is the list of snapshots actually taken into account.
nnout = nout_fi- nout_ini + 1
# final output are in 80 nouts. 
# Missing nouts will have 0 value. 
# This is clearer. 
snout = str(nout_fi).zfill(3)

f_halo = "halo_sub.pickle"
with open(work_dir + f_halo, 'rb') as f:
    halos = pickle.load(f)

info = load.info.Info(nout=nout_fi, base=work_dir)
info.read_info()

s = load.sim.Sim(nout=nout_fi, base=work_dir, dmo=True)
ptype = ["dm mass id pos"]
s.add_part(ptype)

region = smp.set_region(centers=[0.5]*3, radius=0.4)
s.set_ranges(ranges=region["ranges"])
#s.show_ranges()

s.part.load(verbose=False)
part = s.part.dm
part["m"] = part["m"] *info.msun

# Draw particle map to check.
scale = 2.0
for i in len(halos.data.id):
    h_region = smp.set_region(xc=halos.data.x[i], yc=halos.data.y[i], zc=halos.data.z[i],
                            radius = halos.data.rvir[i] * scale)
    img = zoom.utils.part2den(part, info, region=h_region, npix=800)
    if img is False:
        continue

    sid = str(halos.data.id[i]).zfill(5)
    f_save = work_dir + snout + ptype[0][0] + "halo_" + sid + img.proj + ".pickle"
    img.save_data(fout=f_save)
    fout = work_dir + snout + ptype[0][0] + "halo_" + sid + img.proj + ".png"
    img.plot_2d_den(save=fout, show=show, vmin=1e7, vmax=3e10, dpi=100, zposition=True)

