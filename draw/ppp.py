# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 23:20:08 2015

@author: hoseung
"""

def den2d_parallel(x, y, z, m, npix, region=None,
                    ngp=False, cic=False, tsc=False,
                    vmin=None, vmax=None, verbose=False,
                    norm_integer=False):
    """
        accept view of x,y,z, charge arrays and return 2D projected density map.

        assignment funciton (ngp, cic, tsc) takes these arrays.
        But no modification is required.

    """
    from utils import assign
    import multiprocessing
#    from importlib import reload
#    reload(assign)

    # parameters
    mass_test = True
    nprocs = 2

   # physical size
    lbox = max([x.ptp(), y.ptp()])

    dx = lbox * s.part.info.pboxsize
    # As info is per snapshot, it is better to let part infer info; part.info
    ppdx = npix / dx  # pixel per dx. physical density is density * npix^2 / dx^2
    if verbose:
        print("lbox:", lbox)
        print("dx:", dx)

    if region is not None:
        ind_ok = np.where((x > region.xr[0]) & (x < region.xr[1]) &
                          (y > region.yr[0]) & (y < region.yr[1]) &
                          (z > region.zr[0]) & (z < region.zr[1]) )
        x, y, z = x[ind_ok], y[ind_ok], z[ind_ok]

    # scaile to 0 ~ npix.
    y = (y - y.min()) / lbox * npix  # 0 to < npix
    x = (x - x.min()) / lbox * npix  # 0 to < npix
    #z = (z - z.min()) / z.ptp() * npix  # peak to peak

    # mass in solar mass unit.
    # Test if that's true.
    if not mass_test:
        if m.min() < 1:
            print("mass seems not in solar mass unit.")
            print("Minmum mass is {} and maximum is {}".format(m.min(), m.max()))
    #        return False

    # mass assigning in parallel!


    for i in range(nprocs):
        p = multiprocessing.Process(
                target=assign.cic,
                args=(nums[chunksize * i:chunksize * (i + 1)],
                      out_q))
        procs.append(p)
        p.start()

    field = assign.cic(m, x, npix, y, npix, wraparound=True,
                           average=False, norm_integer=False)

    print("minmax field", field.min(), field.max())

    # From this line on, operations on values above the cut only.
    # vlaues into physical units. solar mass / kpc^2
    val_ok = np.logical_and(field > vmin, field < vmax)
    field[val_ok] = field[val_ok] * ppdx * ppdx
    # normalize into solar mass / kpc^2
    print("minmax field after crop and converting into physical unit", field[val_ok].min(), field[val_ok].max())

    # scale to log
    minval = field[val_ok].min()
    maxval = field[val_ok].max()
    #print(minval, maxval)
    return(field)

