# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 15:12:04 2015

@@author: hoseung
"""
import os
base = os.path.abspath('./')
cluster_name = base.split('/')[-1]
base = base + '/'
frefine= 'refine_params.txt'
fnml = input("type namelist file name (enter = cosmo_200.nml):")
if fnml =="":
    fnml = 'cosmo_200.nml'
nout_ini=int(input("Starting nout?"))
nout_fi=int(input("ending nout?"))
scale = input("Scale?: ")
if scale=="":
    scale = 0.3
scale = float(scale)

npix = input("npix (enter = 400)")
if npix == "":
    npix = 400
npix = int(npix)

ptype=["star pos mass"]
refine_params = True
dmo=False
draw=True
draw_halos=False
draw_part = False
draw_hydro = True
if draw_hydro:
    lmax=input("maximum level")
    if lmax=="":
        lmax=19
    lmax = int(lmax)



import load
import utils.sampling as smp
import utils.match as mtc
import draw
import pickle

nouts = range(nout_ini,nout_fi)

for nout in nouts:
    snout = str(nout).zfill(3)

    if refine_params:
        rr = load.info.RefineParam()
        rr.loadRegion(base + frefine)
      
        nn = load.info.Nml()
        nn.load(base + fnml)
        
        aexp = nn.aout[nout]
        i_aexp = mtc.closest(aexp, rr.aexp)
        
        x_refine = rr.x_refine[i_aexp]
        y_refine = rr.y_refine[i_aexp]
        z_refine = rr.z_refine[i_aexp]
        r_refine = rr.r_refine[i_aexp] * 0.5
        
        region = smp.set_region(xc = x_refine, yc = y_refine, zc = z_refine,
                                radius = r_refine * scale)
    else:   
        region = smp.set_region(xc=0.5, yc=0.5, zc=0.5, radius=0.1)                            
    #region = s.part.search_zoomin(scale=0.5, load=True)
    s = load.sim.Sim(nout, base, dmo=dmo)
    s.set_ranges(region["ranges"])
    imgs = draw.img_obj.MapSet(info=s.info, region=region)


    imgp = draw.img_obj.MapImg(info=s.info, proj='z', npix=npix, ptype=ptype)
    imgp.set_region(region)

    #%%
    if draw_part:
        s.add_part(ptype)
        s.part.load()
        part = getattr(s.part, s.part.pt[0])
    
        x = part['x']
        y = part['y']
        z = part['z']
        m = part['m'] * s.info.msun # part must be normalized already!
    
        imgp.set_data(draw.pp.den2d(x, y, z, m, npix,
                                    region=region, cic=True, norm_integer=True))
        imgs.ptden2d = imgp
#    imgp.show_data()

    #%%
    if draw_hydro:
        s.add_hydro()
        s.hydro.amr2cell(lmax=lmax)
        field = draw.pp.pp_cell(s.hydro.cell, npix, s.info, verbose=True)
        ptype = 'gas_den'
        imgh = draw.img_obj.MapImg(info=s.info, proj='z', npix=npix, ptype=ptype)
        imgh.set_data(field)
        imgh.set_region(region)
    #    imgh.show_data()
        imgs.hydro = imgh
    
    #%%
    fdump = base + snout + 'map.pickle'
    with open(fdump, 'wb') as f:
        pickle.dump(imgs, f)
    if draw:
        if draw_part:
            imgs.ptden2d.plot_2d_den(save= base + cluster_name + snout +'star.png', dpi=400, show=False)
        if draw_hydro:    
            imgs.hydro.plot_2d_den(save= base + cluster_name +snout + 'hydro.png',vmax=15,vmin=10, show=False,
                                   dpi=400)
        
