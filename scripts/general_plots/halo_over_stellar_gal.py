# -*- coding: utf-8 -*-
"""
Created on Thu May 28 22:10:36 2015

halos over stellar density map


@author: hoseung
"""
def extract_halos_within(halos, ind_center, scale=1.0):
    import numpy as np
    import utils.sampling as smp
    '''
    Returns indices of halos within SCALE * Rvir of the central halo.

    halos : halo finder output (single snapshot)
    ind_center : index of central halo
    scale : multiplying factor to the Rvir of the central halo
    '''
    xc = halos.data['x'][i_center]
    yc = halos.data['y'][i_center]
    zc = halos.data['z'][i_center]
    rvir= halos.data['rvir'][i_center]

    xx = halos.data['x']
    yy = halos.data['y']
    zz = halos.data['z']
#    m = np.array(halos.data['mvir'])

    dd = smp.distance_to([xc,yc,zc], [xx,yy,zz])
    print(dd)

#    Mcut = 1e10
#    i_m = m > Mcut
#    i_ok = np.logical_and(dd < (rvir * scale), i_m)
    return np.where(dd[0] < (rvir * scale))[0]


import numpy as np
import load
import utils.sampling as smp
import utils.match as mtc
import draw
import pickle
import tree


base = './29172/'
frefine= 'refine_params.txt'
fnml = 'cosmo_200.nml'
scale = 0.3
ptype=["star pos mass"]
hydro = True
npix = 800
refine_params = False
dmo=False
draw=True



nouts = range(187,188)

for nout in nouts:
    snout = str(nout).zfill(3)
        
    hall = tree.halomodule.Halo(nout=nout, base=base, halofinder="HM")  

    if refine_params:
        rr = load.info.RefineParam()
        rr.loadRegion(base + frefine)
      
        nn = load.info.Nml()
        nn.loadNml(base + fnml)
        
        aexp = nn.aout[nout]
        i_aexp = mtc.closest(aexp, rr.aexp)
        
        xc = rr.x_refine[i_aexp]
        yc = rr.y_refine[i_aexp]
        zc = rr.z_refine[i_aexp]
        radius = rr.r_refine[i_aexp] * 0.5
    else:
        cluster = hall.data[np.argmax(hall.data["np"])]
        xc=cluster["x"]#/s.info.cboxsize
        yc=cluster["y"]#/s.info.cboxsize
        zc=cluster["z"]#/s.info.cboxsize
        radius=cluster["rvir"]#/s.info.cboxsize
        #print(xc,yc,zc,radius)
    region = smp.set_region(xc = xc, yc = yc, zc = zc,
                            radius = radius * scale)
                            
    # convert to code unit.
    #hall.normalize()
    
    # subset of halos ONLY inside zoom-in region
    i_center = np.argmax(hall.data['np'])
    h_ind = extract_halos_within(hall, i_center)
    
    #%%
    s = load.sim.Sim(nout, base)
    region = smp.set_region_multi(xc=hall.data[h_ind]["x"],
                            yc=hall.data[h_ind]["y"],
                            zc=hall.data[h_ind]["z"],
                            radius =hall.data[h_ind]["rvir"]* scale)
    s.set_ranges(region["ranges"])
    imgs = draw.img_obj.MapSet(info=s.info, region=region)

    #%%
    s.add_part(ptype)
    s.part.load()
    part = getattr(s.part, s.part.pt[0])
    
    x = part['x']
    y = part['y']
    z = part['y']
    m = part['m'] * s.info.msun # part must be normalized already!
    
    #%%
    imgp = draw.img_obj.MapImg(info=s.info, proj='z', npix=npix, ptype=ptype)
    imgp.set_data(draw.pp.den2d(x, y, z, m, npix, s.info, cic=True, norm_integer=True))
    imgp.set_region(region)
    imgs.ptden2d = imgp
#    imgp.show_data()

    #%%
    if hydro:
        s.add_hydro()
        s.hydro.amr2cell(lmax=14)
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
        imgs.ptden2d.plot_2d_den(save='067star.png', dpi=400, show=False)
        imgs.hydro.plot_2d_den(save='067hydro.png',vmax=12,vmin=10, show=False,
                               dpi=400)
