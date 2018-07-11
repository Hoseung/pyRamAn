# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 14:01:14 2015

@author: hoseung
"""
import numpy as np
import matplotlib.pyplot as plt

def circle_scatter(ax, x_array, y_array, radii_array,
                   colors=None,
                   cmap='RdYlBu_r', **kwargs):
    """
    draws circles of gievn x,y, and radii

    To Do:
        zip only works with iterables (not with a single value, int or float.)
        Can I make it more general?
    """
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
    from draw import img_obj

    if region is None:
        import utils.sampling as smp
        region = smp.Region( **{"xr":(part['x'].min() + offset[0], part['x'].max() + offset[0]),
                                "yr":(part['y'].min() + offset[1], part['y'].max() + offset[1]),
                                "zr":(part['z'].min() + offset[2], part['z'].max() + offset[2])})
        ind_ok = np.arange(len(part['x']))
    else:
        ind_ok = np.where((part["x"] > region.xr[0] - offset[0])
                        & (part["x"] < region.xr[1] - offset[0])
                        & (part["y"] > region.yr[0] - offset[1])
                        & (part["y"] < region.yr[1] - offset[1])
                        & (part["z"] > region.zr[0] - offset[2])
                        & (part["z"] < region.zr[1] - offset[2]))

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
    #import pyximport; pyximport.install()
    from utils import assign
#   vmax = 1e12# in solar mass / kpc^2

#   range = in data unit.   # prange = in physical unit
    if region is not None:
        lbox = 2.0 * region.radius

        buffer = 0.0001 * lbox

        ind_ok = np.where((x > region.xr[0] + buffer) & (x < region.xr[1] - buffer) &
                          (y > region.yr[0] + buffer) & (y < region.yr[1] - buffer) &
                          (z > region.zr[0] + buffer) & (z < region.zr[1] - buffer))[0]

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


def update_tick_labels(ax):
    proj = ax.pp_hal_meta.proj
    npix = ax.pp_hal_meta.npix

    if proj == "x":
        axis1 = "z"
        axis2 = "y"
        xrn, yrn, zrn ="zr", "yr", "xr"
    elif proj == "y":
        axis1 = "x"
        axis2 = "z"
        xrn, yrn, zrn="xr", "zr", "yr"
    elif proj == "z":
        axis1 = "x"
        axis2 = "y"
        xrn, yrn, zrn="xr", "yr", "zr"

    axis1_range = ax.pp_hal_meta.region[xrn]
    axis2_range = ax.pp_hal_meta.region[yrn]

    ax.set_xlabel(axis1 + " [Mpc]")
    ax.set_xticks(np.linspace(0,npix,5))
    xticks = ["{:.1f}".format(x) for x in np.linspace(axis1_range[0],
                                                      axis1_range[1],
                                                      num=5)]
    ax.set_xticklabels(xticks)

    ax.set_ylabel(axis2 + " [Mpc]")
    ax.set_yticks(np.linspace(0,npix,5))
    yticks = ["{:.1f}".format(y) for y in np.linspace(axis2_range[0],
                                                      axis2_range[1],
                                                      num=5)]
    ax.set_yticklabels(yticks)


def pp_halo(h, npix, rscale=1.0, region=None, ind=None, ax=None,
            name=False, radius="rvir",
            verbose=False, new_ax=False,
            fontsize=10, linewidth=1.0,
            color_field=None,
            color_log=False,
            vmin=None, vmax=None,
            keep_clim=True,
            cmap="RdYlBu_r",
            proj="z",
            **kwargs):
    """
    plot halos as circles on the current/given/new ax.

    Parameters
    ----------
    h : tree.halomodule.Halo instance OR Halo.data is also acceptable.
    npix : int
        number of pixels. (assuming square image)
    rscale : float, default = 1.0
        Radii of halos are magnified by rscale.
    region : Optional[dict; see utils.sampling.set_region]
        If region is given, only halos inside the region are plotted.
        If region and ind are both given,
    radius : str, (rvir by default)
        Name of 'radius' field.
    ind : int array
        If ind is given,
        only selected halos are plotted out of the bigger halo sample.
    ax : pyplot ax instance
        If ax is given, halos are plotted on the ax.
    new_ax : Boolean, default = False
        If True, make a new ax and plot halos onto it.
        If ax is None and new_ax is False, then ax = plt.gca()
    color_field: str
        Name of field by whch value halos are colored.
    color_log : Boolean
        If True, coloring by color_field is in log scale.
    keep_clim: Boolean
        If True, keep the color range of given axis
        so that halos overplotted on top of a prexisting axis
        has a consistent color scheme.)

    Both raw halo catalog (x=[0,1]) and tree-like catalog (x=[0,pboxsize]) are OK
    as long as position and radius are in the same unit.
    ->  xrange, yrange are determined relatively.

    Example
    -------
    >>> import tree.halomodule as hmo
    >>> import matplotlib.pyplot as plt
    >>> import draw
    >>> hh = hmo.Halo(nout=187, is_gal=True)
    >>> ax = draw.pp.pp_halo(hh, 400, edgecolor='blue')
    >>> plt.show()
    >>> import tree.treemodule as tmo
    >>> gt = tmo.load_tree('./', is_gal=True)
    >>> ax = draw.pp.pp_halo(gt.data[gt.data["nout"]==187], 400, edgecolor='green', rscale=1e-3)# rvir in kpc unit.
    >>> plt.show()

    Notes
    -----
    1. Region does NOT modify x,y labels.
       better to modify labels outside.
    2. Unless only halos are being plotted, it is better to pass a region.

    TODO
    ----
    Convert length unit into a specific unit so that halos and trees can be overplotted.
    But, what about tree being dependant on the expansion factor?
    """

    class ax_meta():
        def __init__(self, proj, npix):
            self.region=None
            self.proj=proj
            self.npix=npix
            self.xrange=None
            self.yrange=None
            pass

	# h can be either a halo.data, a halo instance, or a tree array.
    try:
        try:
            h.data['x']
        except:
            h.data['xc']
        hd = h.data
    except:
        hd = h

    # both x and xc are fine.
    posname=0
    try:
        hd["x"]
        xn_org, yn_org, zn_org = "x", "y", "z"
    except:
        hd["xc"]
        xn_org, yn_org, zn_org = "xc", "yc", "zc"

    if proj == "x":
        xn,yn,zn = zn_org, yn_org, xn_org
        xrn, yrn, zrn ="zr", "yr", "xr"
    elif proj == "y":
        xn,yn,zn = xn_org, zn_org, yn_org
        xrn, yrn, zrn="xr", "zr", "yr"
    elif proj == "z":
        xn,yn,zn = xn_org, yn_org, zn_org
        xrn, yrn, zrn="xr", "yr", "zr"

    # use pp_halo rather then the script below.
    if ax is None:
        ax = plt.gca()

    if not hasattr(ax, "pp_hal_meta"):
        ax.pp_hal_meta = ax_meta(proj, npix)
    elif ax.pp_hal_meta.region is not None:
        if region is None:
            # If the given ax already has a region defined, and
            # no new region is explicitly given, use the region
            # from the ax IF the projection is the same.
            region = ax.pp_hal_meta.region
        else:
            print("ax has region, but a new region is given")

    if ind is None:
        if region is None:
            # If no region, plot all
            ind = np.arange(len(hd))
            xmin = min(hd[ind][xn] - hd[ind][radius])
            ymin = min(hd[ind][yn] - hd[ind][radius])
            zmin = min(hd[ind][zn] - hd[ind][radius])
            xmax = max(hd[ind][xn] + hd[ind][radius])
            ymax = max(hd[ind][yn] + hd[ind][radius])
            zmax = max(hd[ind][zn] + hd[ind][radius])

            xspan= xmax - xmin
            yspan= ymax - ymin
            zspan= zmax - zmin

        else:
            # If a reion is given, plot only  the halos inside the region.
            # The size of region is retained.
            # image area does not shrink to fit only valid halos.
            ind = np.where( (hd[xn] > getattr(region,xrn)[0]) &
                            (hd[xn] < getattr(region,xrn)[1]) &
                            (hd[yn] > getattr(region,yrn)[0]) &
                            (hd[yn] < getattr(region,yrn)[1]) &
                            (hd[zn] > getattr(region,zrn)[0]) &
                            (hd[zn] < getattr(region,zrn)[1]))[0]
            #print("number of new halos to plot", ind)
            # Make iterable
            if ind is 0:
                hd = np.array([hd])

            xmin = getattr(region,xrn)[0]
            ymin = getattr(region,yrn)[0]
            zmin = getattr(region,zrn)[0]
            xspan = np.ptp(getattr(region,xrn))
            yspan = np.ptp(getattr(region,yrn))
            zspan = np.ptp(getattr(region,zrn))

    else:
        # if ind is a boolean array, convert it to an index array.
        if ind.dtype == 'bool':
            ind = np.arange(len(ind))[ind]

        if region is None:
            # in the original direction.
            xmin = min(hd[ind][xn] - hd[ind][radius])
            ymin = min(hd[ind][yn] - hd[ind][radius])
            zmin = min(hd[ind][zn] - hd[ind][radius])
            xmax = max(hd[ind][xn] + hd[ind][radius])
            ymax = max(hd[ind][yn] + hd[ind][radius])
            zmax = max(hd[ind][zn] + hd[ind][radius])

            xspan= xmax - xmin
            yspan= ymax - ymin
            zspan= zmax - zmin

        else:
            xmin = getattr(region,xrn)[0]
            ymin = getattr(region,yrn)[0]
            zmin = getattr(region,zrn)[0]
            xspan = np.ptp(getattr(region,xrn))
            yspan = np.ptp(getattr(region,yrn))
            zspan = np.ptp(getattr(region,zrn))


    if ax.pp_hal_meta.region is None:
        import utils.sampling as smp
        # Keep a physical region so that I can plot multiple set of halos
        # in a ax consistently.
        ax.pp_hal_meta.region = smp.Region(**{xrn:(xmin, xmin+xspan),
                                              yrn:(ymin, ymin+yspan),
                                              zrn:(zmin, zmin+zspan)})

    x = (hd[ind][xn] - xmin) / xspan * npix
    y = (hd[ind][yn] - ymin) / yspan * npix
    r =  hd[radius][ind]/xspan * npix * rscale # Assuing xspan == yspan

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
            ccc = np.log10(hd[ind][color_field])
        else:
            ccc = hd[ind][color_field]
        if vmin is None:
            vmin = ccc.min()
        if vmax is None:
            vmax = ccc.max()
        # normalize colors
        ccc = 256 * (ccc-vmin) /(vmax-vmin)
        if not (hasattr(ax, "clim") and keep_clim):
            ax.clim = vmin, vmax

        print("MinMax ccc", ccc.min(), ccc.max())
        kwargs.update({"colors": ccc})

    circle_scatter(ax, x, y, r, facecolors='none',
                      linewidth=linewidth, cmap=cmap, **kwargs)

    #if hasattr(ax,"pp_hal_meta"):
    #    if ax.pp_hal_meta.xrange is None:
    #        ax.pp_hal_meta.xrange = [min(x), max(x)]
    #        ax.pp_hal_meta.yrange = [min(y), max(y)]
    ax.set_xlim([0,npix])
    ax.set_ylim([0,npix])

    if name:
        for i, ii in enumerate(ind):
            ax.annotate(str(hd[ii]["id"]), (x[i],y[i]),
                              fontsize=fontsize)
    return ax


def resize(X,shape=None):
    from scipy.ndimage import zoom
    shape = [shape[0]/X.shape[0], shape[1]/X.shape[1]]

    return zoom(np.nan_to_num(X), shape)

def resize_deprecated(X,shape=None):
    """
    Deprecated since 2017.06.26

    resizes N-D array X into shape.
    For example,
    2D 100 * 100 array -> 200 * 200
    new = resize(img, [200,200])

    NOTE
    ----
        Doesn't work if the image is not rectangular.
    """
    if shape==None:
        return X
    m,n = shape
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
            hvar="rho", field_var=None,
            do_resize=True):
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
    from draw import ppc

    sig = 1.0
    sigrange = sig * 2# what is sigrange?

    if proj=="z":
        dim1 = "x"
        dim1r = "xr"
        dim2 = "y"
        dim2r = "yr"
        dim2_rev = 1

    elif proj == "y":
        dim1 = "x"
        dim1r = "xr"
        dim2 = "z"
        dim2r = "zr"
        dim1_rev = 1
        dim2_rev = -1

    elif proj == "x":
        dim1 = "z"
        dim1r = "zr"
        dim2 = "y"
        dim2r = "yr"
        dim1_rev = -1
        dim2_rev = 1
        # apply reverse in imshow stage.

    x = cell[dim1]
    y = cell[dim2]
#    di = [0, 1]  # What is this?
#    zh = xh[2]  # what is this?
    xmi0 = xmin
    xma0 = xmax
    ymi0 = ymin
    yma0 = ymax

    if (xmin is None) & (xmax is None) & (ymin is None) & (ymax is None):
        if region is not None:
            xmi0, xma0 = getattr(egion,dim1r)#[0], region[dim1r][1]
            ymi0, yma0 = getattr(egion,dim2r)#region[dim2r][0], region[dim2r][1]
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

    xl = x - cell['dx']*0.5*sigrange # array as long as x
    xr = x + cell['dx']*0.5*sigrange
    yl = y - cell['dx']*0.5*sigrange
    yr = y + cell['dx']*0.5*sigrange

    maxdx = max(cell['dx'])

# Assuming no slice.
    tol = maxdx #
    #print("maxdx", tol)
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
                sden = cell["var4"][val]*dx*scale_T2#/cell["var0"][val]
            if hvar == "metal":
                sden = cell["var5"][val]*dx*cell["var0"][val]/0.02
        else:
            sden = cell[field_var][val]

    mass = cell["var0"][val]*dx
    mass.transpose()

    mindx = min(dx)
    if verbose:
        print("mindx", mindx)

    xmi = np.floor(xmi0/mindx)*mindx
    xma = np.ceil(xma0/mindx)*mindx
    nx = np.round((xma-xmi-mindx)/mindx).astype(np.int32)

    ymi = np.floor(ymi0/mindx)*mindx
    yma = np.ceil(yma0/mindx)*mindx#ymi + mindx*nx + mindx
    ny = np.round((yma-ymi-mindx)/mindx).astype(np.int32)

    #nx=ny=max([nx,ny])
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
    iin = np.where((ixr >= 0) * (ixl <= nx-1) * (iyr >= 0) * (iyl <= ny-1))[0].astype(np.int32)
    #print(len(iin))
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

    # if ppc.col_over_denom throw an type missmatch error,
    # compile the ppc module locally once more.
    if do_resize:
        return resize(ppc.col_over_denom(iin,
                ixl, ixr, iyl, iyr,
                mass, sden,
                nx, ny, column), [npix,npix])
    else:
        return ppc.col_over_denom(iin,
            ixl, ixr, iyl, iyr,
            mass, sden,
            nx, ny, column)
