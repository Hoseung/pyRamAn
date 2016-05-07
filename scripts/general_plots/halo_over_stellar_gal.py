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

    def extract_halos_within(halos, ind_center, scale=1.0)
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

#    Mcut = 1e10
#    i_m = m > Mcut
#    i_ok = np.logical_and(dd < (rvir * scale), i_m)
    print(dd[0])
    print(rvir)
    i_ok = np.where(dd[0] < (rvir * rscale))[0]

    return i_ok


base = './'
frefine= 'refine_params.txt'
fnml = 'cosmo_200.nml'
scale = 0.3
ptype=["star pos mass"]
hydro = True
npix = 400
refine_params = True
dmo=False
draw=True

import load
import utils.sampling as smp
import utils.match as mtc
import draw
import pickle
import tree


nouts = range(36,38)

for nout in nouts:
    snout = str(nout).zfill(3)

    if refine_params:
        rr = load.info.RefineParam()
        rr.loadRegion(base + frefine)
      
        nn = load.info.Nml()
        nn.loadNml(base + fnml)
        
        aexp = nn.aout[nout]
        i_aexp = mtc.closest(aexp, rr.aexp)
        
        x_refine = rr.x_refine[i_aexp]
        y_refine = rr.y_refine[i_aexp]
        z_refine = rr.z_refine[i_aexp]
        r_refine = rr.r_refine[i_aexp] * 0.5
        
        region = smp.set_region(xc = x_refine, yc = y_refine, zc = z_refine,
                                radius = r_refine * scale)
                                
        hall = tree.halomodule.Halo(nout=nout, base=wdir, halofinder="HM", info=info)    
        hall.load()
        # convert to code unit.
        hall.normalize()
        
        #%%
        
        
        #%%
        # subset of halos ONLY inside zoom-in region
        i_center = np.where(hall.data['np'] == max(hall.data['np']))
        h_ind = extract_halos_within(hall, i_center)
        #print(h_ind)
        h = tree.halomodule.Halo()
        h2 = tree.halomodule.Halo()
        h2.derive_from(hall, h_ind)
        h.derive_from(hall, [h_ind])
        
        
        #%%
        s = load.sim.Sim(nout, wdir)
        region = set_region(h.data.x, yc=h.data.y, zc=h.data.z, radius = h.data.rvir * rscale)

#%%%%%%%%%%%%%%%%%%%%%%

    else:   
        region = smp.set_region(xc=0.5, yc=0.5, zc=0.5, radius=0.1)                            
    #region = s.part.search_zoomin(scale=0.5, load=True)
    s = load.sim.Sim(nout, base, dmo=dmo)
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
