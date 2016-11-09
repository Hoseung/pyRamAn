# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 14:01:14 2015

@author: hoseung
"""

def circle_scatter(ax, x_array, y_array, radii_array,
                   colors=None,
                   cmap='RdYlBu_r', **kwargs):
    """
    draws circles of gievn x,y, and radii
    
    To Do:
        zip only works with iterables (not with a single value, int or float.)
        Can I make it more general? 
    """
    import matplotlib.pylab as plt
    from matplotlib.collections import PatchCollection

    mypatches = []
    try:
        for (x, y, r) in zip(x_array, y_array, radii_array):
            circle = plt.Circle((x,y), radius=r)#, **kwargs)
            mypatches.append(circle)
            
    except:
        circle = plt.Circle((x_array,y_array), radius=radii_array)#, **kwargs)
        mypatches.append(circle)
        
    color_map = plt.get_cmap(cmap)
 
    p = PatchCollection(mypatches, alpha=1.0, **kwargs)
                        #facecolors = color_map(colors))
    if colors is not None:
        p.set_edgecolor(color_map(colors))
        #p.set_facecolor('none')
        #p.set_clim([np.min(colors),np.max(colors)])
    #plt.colorbar(p,shrink=0.5)
    ax.add_collection(p)
    return True

def square_scatter(ax, x_array, y_array, size_array, **kwargs):
    """
    draws square of given x,y, and lengths of side.
    """
    import matplotlib.pylab as plt
    for (x, y, size) in zip(x_array, y_array, size_array):
        square = plt.Rectangle((x-size/2,y-size/2), size, size, **kwargs)
        ax.add_patch(square)
    return True

def polygon_scatter(ax, x_array, y_array, resolution=5, radius=0.5, **kwargs):
    ''' resolution is number of sides of polygon '''
    from matplotlib.patches import CirclePolygon
    for x, y in zip(x_array, y_array):
        polygon = CirclePolygon((x,y), radius=radius, resolution=resolution, **kwargs)
        ax.add_patch(polygon)
    return True

def part2den(part, info, region=None, proj='z', npix=800, ptype=None,
               offset=[0,0,0], hist=False, **kwargs):
    """
    creates img object, calculate 2d-projection map and return the img object.
    given part, info, region, it passes x,y,z,m,npix,info arrays to pp.den2d.
    mass [Msun]
    """

    import numpy as np
    from draw import img_obj
    
    if region is None:
        import utils.sampling as smp
    
        region = smp.set_region(xr=[part['x'].min() + offset[0], part['x'].max() + offset[0]],
                                yr=[part['y'].min() + offset[1], part['y'].max() + offset[1]],
                                zr=[part['z'].min() + offset[2], part['z'].max() + offset[2]])       
        ind_ok = np.arange(len(part['x']))
    else:
        ind_ok = np.where((part.x > region["xr"][0] - offset[0])
                        & (part.x < region["xr"][1] - offset[0])
                        & (part.y > region["yr"][0] - offset[1])
                        & (part.y < region["yr"][1] - offset[1])
                        & (part.z > region["zr"][0] - offset[2])
                        & (part.z < region["zr"][1] - offset[2]))
    
    if len(ind_ok) > 0:
        img = img_obj.MapImg(info=info, proj=proj, npix=npix, ptype=ptype)
        img.set_region(region)
        img.set_data(den2d(part["x"][ind_ok] + offset[0],
                           part["y"][ind_ok] + offset[1],
                           part["z"][ind_ok] + offset[2],
                           part["m"][ind_ok], npix, region=None, proj=proj,
                           cic=True, norm_integer=True, hist=hist, **kwargs))
        return img
    else:
        return False


def den2d(x, y, z, m, npix, region=None, proj='z',
                    ngp=False, cic=False, tsc=False,
                    vmin=None, vmax=None, verbose=False,
                    norm_integer=False, mass_test=True,
                    hist=True):
    """
    Return 2D density map of given particle distribution in (npix,npix) array.
    
    Parameters
    ----------
    x : float array
        x position of particles
    y : float array
        y position of particles            
    z : float array
        z position of particles        
    m : float array
        mass (or any types of charge) of particles
    npix : int
        number of pixel on a side.
    region : Optional[dict; see utils.sampling.set_region]
        If given, only particles inside the region is taken into account.
    vmin : float
        minimum value of map, below which are cut-off.
    vmax : float
        maximum value of map, above which are cut-off.
    
    NOTE
    ----
    Region must be in double precision.
    If single is given max(x) - min(x) > 2 * region['radius'] is possible,
    causing error in assignment calculation.
    """
    import pyximport; pyximport.install()
    from utils import assign
    import numpy as np
#   vmax = 1e12# in solar mass / kpc^2

#   range = in data unit.   # prange = in physical unit
    if region is not None:
        lbox = 2.0 * region['radius']        
        
        buffer = 0.0001 * lbox        
        
        ind_ok = np.where((x > region["xr"][0] + buffer) & (x < region["xr"][1] - buffer) &
                          (y > region["yr"][0] + buffer) & (y < region["yr"][1] - buffer) &
                          (z > region["zr"][0] + buffer) & (z < region["zr"][1] - buffer))[0]
    
        print(len(ind_ok))        
        if len(ind_ok) < 10:
            print("Too few particles to construct a density map")
            return False
    
        m = m[ind_ok]        
        #   image center = region center
        #   position, region in the same unit. 
        if proj == 'z':
            x = ((x[ind_ok] - region['xc']) / lbox + 0.5) * npix
            y = ((y[ind_ok] - region['yc']) / lbox + 0.5) * npix
        elif proj == 'x':
            x = ((x[ind_ok] - region['xc']) / lbox + 0.5) * npix
            y = ((z[ind_ok] - region['zc']) / lbox + 0.5) * npix
        elif proj == 'y':
            x = ((z[ind_ok] - region['zc']) / lbox + 0.5) * npix
            y = ((y[ind_ok] - region['yc']) / lbox + 0.5) * npix

    else:
        if len(x) < 10:
            print("Too few particles to construct a density map")
            return False        
        
        lbox = max([x.ptp(), y.ptp(), z.ptp()])
        # image center = 0.5 * ptp
        if proj == 'z':
            x = (x - min(x)) / lbox * npix
            y = (y - min(y)) / lbox * npix
        elif proj == 'x':
            x = (x - min(x)) / lbox * npix
            y = (z - min(z)) / lbox * npix
        elif proj == 'y':
            x = (z - min(z)) / lbox * npix
            y = (y - min(y)) / lbox * npix


        # mass in solar mass unit.
    # Test if that's true.
    if not mass_test:
        if m.min() < 1:
            print("mass seems not in solar mass unit.")
            print("Minmum mass is {} and maximum is {}".format(m.min(), m.max()))

    # mass assigning
    if hist:
        field, xedges, yedges = np.histogram2d(x, y, bins=[npix,npix])
        field = field.T * m[0]
#        img.set_data(H * part['m'][0]) # Assuming all stellar particle have the same mass
    else:
        if ngp is False and cic is False and tsc is False: cic = True
        if ngp:
            field = assign.ngp(m, x, npix, y, npix, wraparound=True)
        
        print("x.ptp()", x.ptp(), "npix", npix)
 
        if cic:
            field = assign.cic(m, x, npix, y, npix, wraparound=True,
                               average=False, norm_integer=False)
 
        if tsc:
            field = assign.tsc(m, x, npix, y, npix, wraparound=True)

    if verbose:
        print("minmax field", field.min(), field.max())

    
    ppdx = npix / lbox  # pixel per dx. physical density is density * npix^2 / dx^2
    ppddxx = ppdx * ppdx
    
    if vmin is None:
        vmin = field.min() * ppddxx
    if vmax is None:
        vmax = field.max() * ppddxx
    
    val_ok = (field > vmin/ ppddxx) * (field < vmax / ppddxx)
    field[val_ok] = field[val_ok] * ppddxx
    # normalize into solar mass / kpc^2
    if verbose:
        print("minmax field after crop and converting into physical unit", field[val_ok].min(), field[val_ok].max())

    return(field)

def pp_halo(h, npix, rscale=1.0, region=None, ind=None, ax=None,
            name=False, radius="rvir",
            verbose=False, new_ax=False,
            fontsize=10, linewidth=1.0,
            color_field=None,
            color_log=False,
            vmin=None, vmax=None,
            keep_clim=True,
            cmap="RdYlBu_r",**kwargs):
    """
    plot halos as circles on the current/given/new ax.
    
    Parameters
    ----------
    h : tree.halomodule.Halo instance.
        NOT Halo.data. 
    npix : int
        number of pixels. (assuming square image)
    rscale : float, default = 1.0
        Radii of halos are magnified by rscale. 
    region : Optional[dict; see utils.sampling.set_region]
        If region is given, only halos inside the region are plotted.
        If region and ind are both given, 
    ind : int array
        If ind is given, 
        only selected halos are plotted out of the bigger halo sample.
    ax : pyplot ax instance
        If ax is given, halos are plotted on the ax.
    new_ax : Boolean, default = False
        If True, make a new ax and plot halos onto it. 
        If ax is None and new_ax is False, then ax = plt.gca()
    
    Notes
    -----
    1. Region does NOT modify x,y labels. 
       It would be better to modify labels outside.
    2. Unless only halos are being plotted, it is better to pass a region. 


    """
    import matplotlib.pyplot as plt
    from draw import pp
    import numpy as np

    class ax_meta():
        def __init__(self):
            self.region=None
            pass

	# h can be either a halo.data or just a halo.
    try:
        h.data['x']
        hd = h.data
    except:
        hd = h

    # use pp_halo rather then the script below.
    if ax is None:
#        if new_ax:
#            fig = plt.figure()
#            # is it OK to created a new figure object and not return it?
#            ax = fig.add_subplot(111)
#        else:
        ax = plt.gca()
    
    if not hasattr(ax, "pp_hal_meta"):
        ax.pp_hal_meta = ax_meta()
    elif ax.pp_hal_meta.region is not None:
        if region is None:
            # If the given ax already has a region defined, and
            # no new region is explicitly given, use the region 
            # from the ax.
            region =ax.pp_hal_meta.region 

    if ind is None:
        if region is None:
            # If no region, plot all 
            ind = np.arange(len(hd))
            xmin = min(hd['x'][ind])
            ymin = min(hd['y'][ind])
            xspan = np.ptp(hd['x'][ind])
            yspan = np.ptp(hd['y'][ind])
        else:
            # If reion is given, only plot halos inside the region.
            # The size of region is retained. 
            # image area does not shrink to fit only valid halos.
            ind = np.where((hd['x'] > region["xr"][0]) &
                    (hd['x'] < region["xr"][1]) &
                    (hd['y']> region["yr"][0]) & 
                    (hd['y'] < region["yr"][1]) &
                    (hd['z'] > region["zr"][0]) &
                    (hd['z'] < region["zr"][1]))[0]
            xmin = region["xr"][0]
            ymin = region["yr"][0]
            xspan = np.ptp(region["xr"])
            yspan = np.ptp(region["yr"])
    else:
        # if ind is a boolean array, convert it to an index array.
        if ind.dtype == 'bool':
            ind = np.arange(len(ind))[ind]
        
        if region is None:
            xmin = min(hd['x'][ind])
            ymin = min(hd['y'][ind])
            xspan = np.ptp(hd['x'][ind])
            yspan = np.ptp(hd['y'][ind])
        else:
            xmin = region["xr"][0]
            ymin = region["yr"][0]
            xspan = np.ptp(region["xr"])
            yspan = np.ptp(region["yr"])

    x = (hd["x"][ind] - xmin) / xspan * npix 
    y = (hd["y"][ind] - ymin) / yspan * npix 
    r =  hd[radius][ind]/xspan * npix * rscale # Assuing xspan == yspan
    
    import utils.sampling as smp
    
    if ax.pp_hal_meta.region is None:
        # Keep a physical region so that I can plot multiple set of halos 
        # in a ax consistently. 
        ax.pp_hal_meta.region = smp.set_region(xr = (xmin, xmin+xspan),
                                                yr = (ymin, ymin+yspan))

    if verbose:
        print("# of halos to plot:", len(ind))
        print(x)
        print(y)
        print(r)
        print(xmin, ymin, xspan, yspan, npix, rscale)

    if color_field is not None:
        #colors = hd[color_field][ind]
        # normalize!
        # If ax is given, decide whether to keep colormap range or update.
        # Of course, this does not modify patches already added to the ax.
        if color_log:
            ccc = np.log10(hd[color_field][ind])
        else:
            ccc = hd[color_field][ind]
        # normalize colors
        ccc = 256 * (ccc-ccc.min()) /ccc.ptp()
        if not (hasattr(ax, "clim") and keep_clim):
            #ax.clim=None
            ax.clim = ccc.min(), ccc.max()
                    
        print("MinMax ccc", ccc.min(), ccc.max())
        kwargs.update({"colors": ccc})
    #else:
        #kwargs.update({"colors": None})
    pp.circle_scatter(ax, x, y, r, facecolors='none',
                      linewidth=linewidth, cmap=cmap, **kwargs)

    ax.set_xlim([min(x), max(x)])
    ax.set_ylim([min(y), max(y)])
    if name:
        for i, ii in enumerate(ind):
            ax.annotate(str(hd["id"][ii]), (x[i],y[i]),
                              fontsize=fontsize)
    return ax


def resize(X,shape=None):
    """
    resizes N-D array X into shape.
    For example,
    2D 100 * 100 array -> 200 * 200
    new = resize(img, [200,200])
    """
    import numpy as np
    if shape==None:
        return X
    m,n = shape
#    print(shape)
#    print(m,n)
    Y = np.zeros((m,n),dtype=type(X[0,0]))
    k = len(X)
    p,q = k/m,k/n
    for i in range(m):
        Y[i,:] = X[i*p,np.int_(np.arange(n)*q)]
    return Y


def pp_colden(cell, npix, info, proj="z", verbose=False, autosize=False):
    """
    Warning
    -------
    incomplete. 
    column density is simpler (no need to divide by projected mass.)
    """
    import numpy as np
#    xh = np.asarray([0.5, 0.5, 0.5])
    hvar = 1
    sig = 1.0
    sigrange = sig * 2# what is sigrange?

    x = cell.x
    y = cell.y
#    z = cell.z
#    di = [0, 1]  # What is this?
#    zh = xh[2]  # what is this?

    # range
    xmi0 = min(x)
    xma0 = max(x)
    ymi0 = min(y)
    yma0 = max(y)

    xl = x - cell.dx/2*sigrange # array as long as x
    xr = x + cell.dx/2*sigrange
    yl = y - cell.dx/2*sigrange
    yr = y + cell.dx/2*sigrange

    maxdx = max(cell.dx)

# Assuming no slice.
    tol = maxdx #
    val= np.where((xr >= xmi0-tol) & (xl <= xma0+tol) &
                  (yr >= ymi0-tol) & (yl <= yma0+tol))[0]

# Check for valid cell existing.
    if len(val) == 0:
        print("No cell is selected")
        return False

    scale_d = info.unit_d
    scale_l = info.unit_l
#    scale_nH = info.unit_nH
    scale_T2 = info.unit_T2

    dx = cell.dx[val]

# No max, no column
    # mass/L**2 = rho*L
    # sden differs with hydro variable type.
    if hvar == 1:
        sden = cell.var0[val]*dx*(scale_d*scale_l)*0.76/1.66e-24
        #sden = cell.var0[val]**2*dx*scale_nH
    if hvar == 5:
        sden = cell.var4[val]*scale_T2
    if hvar == 6:
        sden = cell.var0[val]*dx*cell.var5[val]/0.02
#    mass = cell.var0[val]*dx

    mindx = min(dx)

    xmi = np.floor(xmi0/mindx)*mindx
    xma = np.ceil(xma0/mindx)*mindx
    nx = np.round((xma-xmi-mindx)/mindx)

    ymi = np.floor(ymi0/mindx)*mindx
    yma = ymi + mindx*nx + mindx
    ny = np.round((yma-ymi-mindx)/mindx)

    if verbose:
        print(" ... working resolution % npix ", nx, ny)
        print(" ... given npix:", npix)

    # adjust npix (IDL Line 173)

    dx_dot = (xma0-xmi0)/npix
#    lv_skip = np.floor(-1. * np.log10(dx_dot)/np.log10(2)) - 1
    # Assum no smoothing.

    ixmi = round(xmi/mindx)
    iymi = round(ymi/mindx)

    xl   = xl[val]
    xr   = xr[val]
    yl   = yl[val]
    yr   = yr[val]

    ixl = np.round(xl / mindx) - ixmi
    ixr = np.round(xr / mindx) - ixmi -1
    iyl = np.round(yl / mindx) - iymi
    iyr = np.round(yr / mindx) - iymi -1

    """
    fd = np.where(ixl < 0)[0]
    if len(fd) > 0:
        ixl[fd] = 0
    fd = np.where(ixr > nx - 1)[0]
    if len(fd) > 0:
        ixr[fd] = nx -1
        fd = np.where(ixl < 0)[0]
    if len(fd) > 0:
        ixl[fd] = 0
    fd = np.where(ixl < 0)[0]
    if len(fd) > 0:
        ixl[fd] = 0
    """
    ixl[ixl < 0] = 0
    ixr[ixr > nx -1] = nx -1
    iyl[iyl < 0] = 0
    iyr[iyr > ny - 1] = ny -1

    fd = np.where((ixl >= 0) & (ixr >=0) & (ixl > ixr))[0]
    if len(fd):
        ixr[fd] = ixl[fd]
    fd = np.where((iyl >= 0) & (iyr >=0) & (iyl > iyr))[0]
    if len(fd):
        iyr[fd] = iyl[fd]

#    mass.transpose()
    sden.transpose()
#    print("shape of sden", sden.shape)
#    finalmap = ppc.col_over_denom(iin, ixl, ixr, iyl, iyr, mass, sden, nx, ny)
    
    # original size != npix * npix. 
    # rescale it to npix * npix image. 
    return(resize(sden.reshape(npix,npix), [npix,npix]))


def pp_cell(cell, npix, info, proj="z", verbose=False, autosize=False,
            column=False,
            region=None,
            xmin=None, xmax=None, ymin=None, ymax=None,
            hvar="rho", field_var=None):
    """
    Accepts cell data and returns 2D projected gas map.
     *position and dx must be in the same unit.

    example
    -------
    >>> gas_map = pp_cell(gal.cell, 200, info)
    >>> plt.imshow(gas_map, origin="lower")
    >>> plt.show()
    -> Without range or region, all cells are taken. 

    >>> region = smp.set_region(centers=[0.,0.,0.], radius=gal.region['radius'])
    >>> gas_map = pp_cell(gal.cell, 200, info, region=region)
    -> cells only inside the region are taken into account.
       *It's a cubic region, not a sphere.


    Note
    ----
    To do: Currently only column sum is supported. maximum value along the column, or column average option are needed.

    """
    import numpy as np
    from draw import ppc
    

    sig = 1.0
    sigrange = sig * 2# what is sigrange?

    x = cell["x"]
    y = cell["y"]
#    z = cell.z
#    di = [0, 1]  # What is this?
#    zh = xh[2]  # what is this?
    if (xmin is None) & (xmax is None) & (ymin is None) & (ymax is None): 
        if region is not None:
            xmi0, xma0 = region['xr'][0], region['xr'][1]
            ymi0, yma0 = region['yr'][0], region['yr'][1]
        else:
            xmi0, xma0, ymi0, yma0 = min(x), max(x), min(y), max(y)
    if xmi0 is None:
        xmi0 = min(x)
    if xma0 is None:
        xma0 = max(x)        
    if ymi0 is None:
        ymi0 = min(y)
    if yma0 is None:
        yma0 = max(y)
       
    xl = x - cell['dx']/2*sigrange # array as long as x
    xr = x + cell['dx']/2*sigrange
    yl = y - cell['dx']/2*sigrange
    yr = y + cell['dx']/2*sigrange

    maxdx = max(cell['dx'])

# Assuming no slice.
    tol = maxdx #
    val = np.where((xr >= xmi0-tol) & (xl <= xma0+tol) &
                  (yr >= ymi0-tol) & (yl <= yma0+tol))[0]

# Check for valid cell existing.
    if len(val) == 0:
        print("No cell is selected")
        return False

    scale_d = info.unit_d
    scale_l = info.unit_l
    scale_nH = info.unit_nH
    scale_T2 = info.unit_T2

    dx = cell["dx"][val]

# No max, no column
    # mass/L**2 = rho*L
    # sden differs with hydro variable type.

    if column:
        sden = cell[hvar][val]*dx*(scale_d*scale_l)*0.76/1.66e-24
    else:
        if field_var is None:
            if hvar == "rho":
                hvar = "var0"
                sden = cell[hvar][val]**2*dx*scale_nH
            if hvar == "temp":
                sden = cell[hvar][val]*scale_T2
            if hvar == "metal":
                sden = cell[havr][val]*dx*cell.var5[val]/0.02            
        else:
            sden = cell[field_var][val]

    mass = cell["var0"][val]*dx
    mass.transpose()

    mindx = min(dx)

    xmi = np.floor(xmi0/mindx)*mindx
    xma = np.ceil(xma0/mindx)*mindx
    nx = np.round((xma-xmi-mindx)/mindx).astype(np.int32)

    ymi = np.floor(ymi0/mindx)*mindx
    yma = ymi + mindx*nx + mindx
    ny = np.round((yma-ymi-mindx)/mindx).astype(np.int32)

    if verbose:
        print(" ... working resolution % npix ", nx, ny)
        print(" ... given npix:", npix)
    # adjust npix (IDL Line 173)

    #dx_dot = (xma0-xmi0)/npix
#    lv_skip = np.floor(-1. * np.log10(dx_dot)/np.log10(2)) - 1
    # Assum no smoothing.
    ixmi = np.round(xmi/mindx).astype(np.int32) # int
    iymi = np.round(ymi/mindx).astype(np.int32)

    xl   = xl[val]
    xr   = xr[val]
    yl   = yl[val]
    yr   = yr[val]

    ixl = np.round(xl / mindx).astype(np.int32) - ixmi
    ixr = np.round(xr / mindx).astype(np.int32) - ixmi -1
    iyl = np.round(yl / mindx).astype(np.int32) - iymi
    iyr = np.round(yr / mindx).astype(np.int32) - iymi -1
    iin = np.where((ixr >= 0) & (ixl <= nx-1) & (iyr >= 0) & (iyl <= ny-1))[0].astype(np.int32)
    # What does it mean?
    
    fd = ixl < 0
    if len(fd) > 0:
        ixl[fd] = 0
    fd = ixr > nx - 1
    if len(fd) > 0:
        ixr[fd] = nx -1
        fd = ixl < 0
    if len(fd) > 0:
        ixl[fd] = 0
    fd = ixl < 0
    if len(fd) > 0:
        ixl[fd] = 0
    
    ixl[ixl < 0] = 0
    ixr[ixr > nx -1] = nx -1
    iyl[iyl < 0] = 0
    iyr[iyr > ny - 1] = ny -1

    fd = np.where((ixl >= 0) & (ixr >=0) & (ixl > ixr))[0]
    if len(fd):
        ixr[fd] = ixl[fd]
    fd = np.where((iyl >= 0) & (iyr >=0) & (iyl > iyr))[0]
    if len(fd):
        iyr[fd] = iyl[fd]
    
#    sden.transpose()
#    colden = np.zeros((nx,ny), dtype=np.float32)
#    denom =  np.zeros((nx,ny), dtype=np.float32)

    if verbose:
        print("dimensions of \n mass = {}, \n sden = {}, and nx, ny are {}, {}".format(\
        mass.shape, sden.shape, nx, ny))
        print(all(np.isfinite(mass)), all(np.isfinite(sden)))
        print(mass[100:110])
        print(sden[100:110])
        
    #print(nx, ny, npix)
    # if ppc.col_over_denom throw an type missmatch error,
    # compile the ppc module locally once more.
    return resize(ppc.col_over_denom(iin,
            ixl, ixr, iyl, iyr,
            mass, sden,
            nx, ny, column), [npix,npix])

