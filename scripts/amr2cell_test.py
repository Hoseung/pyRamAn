# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 02:44:05 2015

@author: hoseung
"""
import ..load
#from load import info
import ..utils.sampling as smp
import ..utils.match as mtc
import time
import ..utils
import ..draw
base = '/home/hoseung/Work/data/Lmax/'
nout=187
snout = str(nout).zfill(3)
fnml = 'cosmo_200.nml'
frefine = 'refine_parmas.txt'

#utils.util.reimport (draw)
npix = 800
ptype=["dm id pos mass"]
"""
nn=load.info.Nml()
nn.loadNml(base + fnml)

rr = load.info.RefineParam()
rr.loadRegion(base + frefine)

aexp = nn.aout[nout]
i_aexp = mtc.closest(aexp, rr.aexp)

x_refine = rr.x_refine[i_aexp]
y_refine = rr.y_refine[i_aexp]
z_refine = rr.z_refine[i_aexp]
r_refine = rr.r_refine[i_aexp] * 0.5

region = smp.set_region(xc = x_refine, yc = y_refine, zc = z_refine,
                        radius = r_refine * scale)
"""
region = smp.set_region(xr=[0.25,0.3], yr=[0.25,0.3], zr=[0.25,0.3])
s = load.sim.Sim(nout=nout,base=base, dmo=True)
s.add_info()
s.set_ranges(ranges=region["ranges"])
#%%
s.add_part(ptype)
s.part.load()
part = getattr(s.part, s.part.pt[0])

x = part['x']
y = part['y']  # These are views. right?
z = part['y']
m = part['m'] * s.info.msun # part must be normalized already!

imgp = draw.img_obj.MapImg(info=s.info, proj='z', npix=npix, ptype=ptype)
imgp.set_data(draw.pp.den2d(x, y, z, m, npix, s.info, cic=True, norm_integer=True))
imgp.set_region(region)
#    imgp.show_data()
imgs = draw.img_obj.MapSet(info=s.info, region=region)    
imgs.ptden2d = imgp

#%%
s.add_hydro()
t0 = time.time()
s.hydro.amr2cell(lmax=12)
print(time.time() - t0)
#%%

field = draw.pp.pp_cell(s.hydro.cell, npix, s.info, verbose=True)


utils.util.reimport (draw.img_obj)

ptype = 'gas_den'
img = draw.img_obj.MapImg(info=s.info, proj='z', npix=100, ptype=ptype)
img.set_data(field)
img.set_region(region)
img.show_data()

#%%
import pickle
fdump = base + snout + ptype +'map.pickle'
with open(fdump, 'wb') as f:
    pickle.dump(img, f)

fin = base + snout + ptype +'map.pickle'
with open(fin, 'rb')as f:
    img = pickle.load(f)
    fout = base + snout + ptype + "map_" + img.proj + ".png"
    img.plot_2d_den(save=fout, show=True, vmin=1e7, vmax=5, dpi=500, log=False)
