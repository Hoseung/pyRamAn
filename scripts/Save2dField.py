# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 17:11:10 2015

Calculates 2D density fields and save in .pickle files.

Keywords:
nout = a list of nouts.
base = './'
ptype = 'dm'
npix = 800,
scale = 0.5
frefine = 'refine_params.txt'
fnml = 'cosmo_200.nml'

:example:

No example.

@author: hoseung
"""
def Save2dField(nout, base='./', ptype='dm', npix=800,
               scale = 0.5, frefine = 'refine_params.txt',
               fnml = 'cosmo_200.nml'):
    import load
    import utils.sampling as smp
    import utils.match as mtc
    from draw import pp
    from draw import img_obj
#    base = '/home/hoseung/Work/data/C01605_mid/'
    snout = str(nout).zfill(3)

    s = load.sim.Sim(nout, base)
    s.add_part(ptype)

    rr = load.info.RefineParam()
    rr.loadRegion(base + frefine)

    nn= load.info.Nml()
    nn.loadNml(base + fnml)

    aexp = nn.aout[nout]
    i_aexp = mtc.closest(aexp, rr.aexp)

    x_refine = rr.x_refine[i_aexp]
    y_refine = rr.y_refine[i_aexp]
    z_refine = rr.z_refine[i_aexp]
    r_refine = rr.r_refine[i_aexp] * 0.5

    region = smp.set_region(xc = x_refine, yc = y_refine, zc = z_refine,
                            radius = r_refine * scale)
    #region = s.part.search_zoomin(scale=0.5, load=True)

    s.set_ranges(region["ranges"])

    # part
    s.part.load()

    part = getattr(s.part, s.part.pt[0])
    x = part['x']
    y = part['y']  # These are views. right?
    z = part['y']
    m = part['m'] * s.info.msun # part must be normalized already!

    img = img_obj.Prj2D(info=s.info, proj='z', npix=npix, ptype=ptype)
    img.set_data(pp.den2d(x, y, z, m, npix, s.info, cic=True, norm_integer=True))
    img.set_region(region)
    #%%
    foutput = base + snout + ptype + 'map.pickle'
    import pickle
    with open(foutput, 'wb') as f:
        pickle.dump(img, f)


#%%
base = "/home/hoseung/Work/data/036370/"
Save2dField(37, base=base, ptype=["star pos mass id"], npix = 800)