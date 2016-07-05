import utils.cosmology
import tree.ctutils as ctu
import load
import numpy as np
import matplotlib.pyplot as plt
import pickle

def gal_file_list(prg_this, galidx, GM_dir = './GalaxyMaker/'):
    """
        Parameters
        ----------
        prgt : 
            Progenitor only tree
        galidx :
            idx of the interested galaxy at nout_fi.
    """
    sgalidx=str(galidx)
    gal_dir = GM_dir + "gal_" + sgalidx + '/' 

    file_list = []
    for nout in prg_this['nout']:
        nout_dir = "GAL_" + str(nout).zfill(5) +"/"
        this_gal = prg_this[prg_this['nout'] == nout]
        if this_gal['phantom'] == 1:
            continue
        galid = this_gal['Orig_halo_id'][0]
        file_list.append(gal_dir + "gal_stars_" + str(nout).zfill(3) + "_" + str(galid).zfill(7))

    return file_list


def convert_gm_to_halo(cat, info):
    """
        tree: 
        pos in Mpc/h [0-200]
    """
    cat['x'] /= info.cboxsize
    cat['y'] /= info.cboxsize
    cat['z'] /= info.cboxsize
    
    cat['rvir'] /= info.boxtokpc
    

def nout2lbt(nout, nout_fi=187):
    import astropy.cosmology as ac
    aexp = 1 - (nout_fi - nout)*0.005
    
    return ac.WMAP7.lookback_time(1/aexp -1).value
        
    
#for ac.WMAP7.lookback_time(1).value

def angular_momentum(star):
    x = star['x']
    y = star['y']
    z = star['z']
    
    vx = star['vx']
    vy = star['vy']
    vz = star['vz']
    
    m = star['m']
    
    lx = np.sum((y*vz - z*vy)*m)
    ly = np.sum((z*vx - x*vz)*m)
    lz = np.sum((x*vy - y*vx)*m)
    
    return lx,ly,lz


def gal_from_GM(gcat, info, nout=None, idgal=None, fname=None, **kwargs):
    from galaxy import galaxy
    if fname is None:
        gm = load.rd_GM.rd_gal(nout, idgal)
    else:
        gm = load.rd_GM.rd_gal(info.nout,0,fname=fname)

    if gm.star is None:
        print("no stellar particles, skipping")
        return False

    gm.dm = None
    gm.cell = None
    gm.star['time'] = utils.cosmology.time2gyr(gm.star['time'],
                                 z_now = info.zred,
                                 info=info)

    gal = galaxy.Galaxy(halo = gcat,
                        info=info)
    
    good_gal = gal.mk_gal(gm.star, gm.dm, gm.cell,
                          unit_conversion="GM",
                          **kwargs)
    return gal


def draw_plot(gal, cat, ax, ax_hist, smoothed_lambda=True):
    #ax = ax.ravel()
    t_min = min(gal.star['time'])
    t_max = max(gal.star['time'])
    bins = np.linspace(t_min, t_max, num=4)
    i_bins = np.digitize(gal.star['time'],bins)
    ax[0].plot(nout2lbt(cat.nouts), cat.data['mstar'])
    ax_hist.hist(gal.star['time'], histtype='step')
    if smoothed_lambda:
        ax[2].plot(nout2lbt(cat.nouts), cat.smoothed_lambda)# cat.data['lambda_r'])
    else:
        ax[2].plot(nout2lbt(cat.nouts), cat.data['lambda_r'])

    ax[3].plot(nout2lbt(cat.nouts), cat.data['reff'])
    if cat.merger is not None:
        for mr, xx, x2 in zip(cat.merger.mr,
                              nout2lbt(cat.merger.nout),
                              nout2lbt(cat.merger.nout_ini)):
            ax[0].axvline(xx, linestyle=':')
            ax[0].annotate("{:.1f}".format(mr), xy=(xx,0.4))
            ax[2].axvline(xx, linestyle=':')
            ax[2].annotate("{:.1f}".format(mr), xy=(xx,0.4))
            ax[2].axvline(xx, linestyle=':')
            ax[2].axvline(x2, linestyle=':', c='g')


    ax[0].set_ylabel(r"$M_{*} [M_{\odot}]$", fontsize=26)
    ax[1].set_ylabel(r"$ j_{x,y,z} [kpc km s^{-1}]$", fontsize=26)
    ax[2].set_ylabel(r"$\lambda$", fontsize=26)
    ax[3].set_ylabel(r"$R_{eff} [kpc]$", fontsize=26)
    ax[4].set_ylabel("Mstar")

#    ax[0].locator_params(nbins=3, axis='y')
#    ax[1].locator_params(nbins=3, axis='y')
#    ax[2].locator_params(nbins=3, axis='y')
#    ax[3].locator_params(nbins=3, axis='y')
#    ax[4].locator_params(nbins=3, axis='y')



def new_old_gals(galid_at_187, prgt):
    galidx = prgt["id"][(prgt["nout"] == 187) * (prgt["Orig_halo_id"] == galid_at_187)][0]
    prg_this = ctu.extract_main_tree(prgt, galidx)
    #prg_this['Orig_halo_id'][[0,100]]
    #prg_this['nout'][[0,100]]

    info1 = load.info.Info(187)
    info2 = load.info.Info(87)

    #gg = prg_now[prg_now["Orig_halo_id"] == 707]
    gg = prg_this[0]
    convert_gm_to_halo(gg, info1)

    kwargs = {"debug":False, "verbose":False}
    gal1 = gal_from_GM(187, gg['Orig_halo_id'], gg, info1, **kwargs)

    gg = prg_this[100]
    convert_gm_to_halo(gg, info2)
    gal2 = gal_from_GM(87, gg['Orig_halo_id'], gg, info2)
    
    return gal1, gal2


def draw_lxyz(ax, t_mean, lvec, **kwargs):
    ax.scatter(t_mean,lvec[0], marker='>', edgecolor='red', **kwargs)
    ax.scatter(t_mean,lvec[1], marker='^', edgecolor='green', **kwargs)
    ax.scatter(t_mean,lvec[2], marker='v', edgecolor='blue', **kwargs)

def add_lxyz_point(ax, gal, info):
    lvec = angular_momentum(gal.star)
    draw_lxyz(ax, nout2lbt(info.nout), lvec, s=100, facecolor='none')
    
def load_prg_gal_data222(gg, info):
    convert_gm_to_halo(gg, info)
    inout = max(prg_this['nout']) - gg['nout']
    idgal = prg_this['Orig_halo_id'][inout]
    if prg_this['phantom'][inout] == 1:
        return False
    nout = prg_this['nout'][inout]
    return gal_from_GM(gg, info, nout=nout, idgal=idgal)


def load_prg_gal_data(gg, info, fname=None):
    convert_gm_to_halo(gg, info)
    if prg_this['phantom'][inout] == 1:
        return False
#    nout = prg_this['nout'][inout]
    return gal_from_GM(gg, info, nout=nout, idgal=idgal, fname=fname)

def gal_from_GM(gcat, info, nout=None, idgal=None, fname=None,
                dm=False, cell=False, **kwargs):
    from galaxy import galaxy
    if fname is None:
        gm = load.rd_GM.rd_gal(nout, idgal)
    else:
        gm = load.rd_GM.rd_gal(info.nout,0,fname=fname)

    if gm.star == None:
        print("no stellar particles, skipping")
        return False

    gm.star['time'] = utils.cosmology.time2gyr(gm.star['time'],
                                 z_now = info.zred,
                                 info=info)
    gcat['x'] /= info.cboxsize
    gcat['y'] /= info.cboxsize
    gcat['z'] /= info.cboxsize
    #print("Cat x", gcat['x'])
    gal = galaxy.Galaxy(halo = gcat,
                        info=info)

    
    if dm:
        fname_dm = fname.replace("stars", "dms")
        gm.dm = load.rd_GM.rd_dm(info.nout, 0, fname = fname_dm)
    else:
        gm.dm = None
    
    if cell:
        fname_cell = fname.replace("stars", "cells")
        gm.cell = load.rd_GM.rd_cell(info.nout, 0, fname = fname_cell)
        # unit of "rho" is not touched.
    else:
        gm.cell = None
        
    good_gal = gal.mk_gal(gm.star, gm.dm, gm.cell,
                          unit_conversion="GM",
                          verbose=False, debug=False,
                          mstar_min = 1e8,
                          **kwargs)
    return gal


def gas_cell_mass(cell, info):
    "Returns mass of each cell in Msun unit"
    msun = 1.98892e33
    rho = cell['rho'] * info.unit_d
    dx = cell['dx'] / (info.pboxsize * 1e3) * info.unit_l
    return rho * dx**3 / msun
    

# Measure gas Mass 
def gas_mass(gal, info):
    from scipy.stats import binned_statistic

    cell = gal.cell
    gasmass = gal.gasmass
    gmass = gas_cell_mass(cell, info)
    
    d_cut = np.percentile(cell['rho'], 90)
    ind = cell['rho'] > d_cut
    gasmass["dense"] = np.sum(gmass[ind])

    dx_unique = np.unique(cell['dx'])
    dx_cut = dx_unique[min(len(dx_unique), 2) -1]
    ind = cell['dx'] <= dx_cut
    gasmass["highlevel"] = np.sum(gmass[ind])

    rr = np.sqrt(np.square(cell['x']) +
                 np.square(cell['y']) +
                 np.square(cell['z']))
    
    ind = np.where(rr < gal.meta.rgal)[0]
    gasmass["rgal"] = np.sum(gmass[ind])
    result = binned_statistic(rr, gmass, 'sum', bins=8)
    gasmass["radial"] = result.statistic
    gasmass['loc'] = 0.5 * (result.bin_edges[:-1] + result.bin_edges[1:])
    
    #return gasmass


def pp_cell(cell, npix, info, proj="z", verbose=False, autosize=False,
            column=False, region=None, sigrange=1,
            xmin=None, xmax=None, ymin=None, ymax=None,
            hvar="rho", field_var=None):
    """
    Accepts cell data and returns 2D projected gas map.
     *position and dx must be in the same unit.

    Parameters
    ----------
    havr: {"rho", "temp", "metal"}
    sigrange : int
        if > 1, smoothing. 

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
    from draw.pp import resize

    field_x, field_y, field_z, \
    field_dx, field_rho,\
    field_vx, field_vy, field_vz, \
    field_temp, field_metal = cell.dtype.names

    sigrange = np.float64(sigrange) #smoothing scale. 1.0 = None

    x = cell[field_x]
    y = cell[field_y]
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

    xl = x - cell.dx/2*sigrange # array as long as x
    xr = x + cell.dx/2*sigrange
    yl = y - cell.dx/2*sigrange
    yr = y + cell.dx/2*sigrange

    maxdx = max(cell.dx)

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

    dx = cell[field_dx][val]

# No max, no column
    # mass/L**2 = rho*L
    # sden differs with hydro variable type.

    if column:
        sden = cell[hvar][val]*dx*(scale_d*scale_l)*0.76/1.66e-24
    else:
        if field_var is None:
            if hvar == "rho":
                sden = cell[field_rho][val]**2*dx*scale_nH
            if hvar == "temp":
                sden = cell[field_temp][val]*scale_T2
            if hvar == "metal":
                sden = cell[field_temp][val]*dx*cell.var5[val]/0.02
        else:
            sden = cell[field_var][val]

    mass = cell[field_rho][val]*dx
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
    return resize(ppc.col_over_denom(iin,
            ixl, ixr, iyl, iyr,
            mass, sden,
            nx, ny, column), [npix,npix])


	
