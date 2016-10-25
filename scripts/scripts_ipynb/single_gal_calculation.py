
# coding: utf-8

# In[27]:

#Now, Add auto size.

def pp_cell(cell, npix, info, proj="z", npix_max = 1200, scale=1.0,
            column=False, region=None, sigrange=1,
            xmin=None, xmax=None, ymin=None, ymax=None,
            hvar="rho", smoothing=False, autosize=True,
            verbose=False, debug=False):
    """
    Accepts cell data and returns 2D projected gas map.
     *position and dx must be in the same unit.

    Parameters
    ----------
    cell: 
    npix: int
    info: info instance
    proj: {"x", "y", "z"}
    column : logical
    region : region dict
        If not None, only cells inside the region are drawn.
    xmin : float
    xmax : float
    ymin : float
    ymax : float
        min, max range of x,y
    hvar: {"rho", "temp", "metal"}
        type of hydro variable to draw
    sigrange : int
        if > 1, smoothing. 
    autosize : logical
        if True, npix is modified so that 1 pixel corresponds to one highest level cell.
    npix_max : int
        maximum possible npix to prevent huge image array, 1200 by default. 
    scale : float
        scale the image on returning. Does not involve additional calculation.  
    verbose: logical
        

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

    field_x, field_y, field_z,     field_dx, field_rho,    field_vx, field_vy, field_vz,     field_temp, field_metal = cell.dtype.names

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

    if hvar == "rho" and column:
        sden = cell[hvar][val]*dx*(scale_d*scale_l)*0.76/1.66e-24
    else:
        if hvar == "rho":
            sden = cell[field_rho][val]**2 * dx * scale_nH
        elif hvar == "temp":
            field = field_temp
            sden = cell[field][val]*scale_T2
        elif hvar == "metal":
            field = field_metal
            sden = cell[field_rho][val] * dx**3 * cell[field][val]/0.02
        elif hvar == 'vx':
            field = field_vx
            sden = cell[field_rho][val] * dx**3 * cell[field][val]
        elif hvar == 'vy':
            field = field_vy
            sden = cell[field_rho][val] * dx**3 * cell[field][val]
        elif hvar == 'vz':
            field = field_vz
            sden = cell[field_rho][val] * dx**3 * cell[field][val]
        else:
            sden = cell[hvar][val]
            # without column density, sden is divided by column mass.
            # pixel-wise average velocity = the sum of density-weighted velocity / density sum
    if hvar in ["rho", "temp"]:
        mass = cell[field_rho][val] * dx
    else:
        mass = cell[field_rho][val] * dx**3
    #mass.transpose()

    mindx = min(dx)
    while True:
        xmi = np.floor(xmi0/mindx)*mindx
        xma = np.ceil(xma0/mindx)*mindx
        nx = np.round((xma-xmi-mindx)/mindx).astype(np.int32)
        ymi = np.floor(ymi0/mindx)*mindx
        yma = ymi + mindx*nx + mindx
        ny = np.round((yma-ymi-mindx)/mindx).astype(np.int32)
        if nx < npix:
            break
        mindx *= 2

    if verbose:
        print(" ... working resolution % npix ", nx, ny)
        print(" ... given npix:", npix)

    if autosize:
        npix = min([nx, npix_max])
        
    # adjust npix (IDL Line 173)

    dx_dot = (xma0-xmi0)/npix
    lv_skip = np.floor(-1. * np.log10(dx_dot)/np.log10(2)) - 1


    if smoothing:
        lv = -np.round(np.log10(dx)/np.log10(2))
        levelmin = min(lv)
        levelmax = max(lv)

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


    if smoothing:
        ixc = (ixl + ixr)/2
        iyc = (iyl + iyr)/2
        idxc = (ixr - ixl+1)

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

    if debug:
        print("dimensions of \n mass = {}, \n sden = {}, and nx, ny are {}, {}".format(        mass.shape, sden.shape, nx, ny))
        #print(all(np.isfinite(mass)), all(np.isfinite(sden)))
        print("mass", mass[10:30])
        print("field name=", field)
        print("field value", cell[field][val][10:30])

    return resize(ppc.col_over_denom(iin,
            ixl, ixr, iyl, iyr,
            mass, sden,
            nx, ny, column), [np.int64(scale * npix), np.int64(scale * npix)])

def get_nout_merger_org(prg_this,  nout_merger):
    nout_merger_org = nout_merger
    while True:
        prg_now = prg_this[prg_this["nout"] == nout_merger_org]
        if prg_now["phantom"] == 0:
            break
        nout_merger_org -=1
    return nout_merger_org, prg_now


# In[24]:

def angular_momentum_gas(cell, info):
    """
    
    """
    x = cell['x']
    y = cell['y']
    z = cell['z']

    vx = cell['vx']
    vy = cell['vy']
    vz = cell['vz']

    m = gas_cell_mass(cell, info)

    lx = np.sum((y*vz - z*vy)*m)
    ly = np.sum((z*vx - x*vz)*m)
    lz = np.sum((x*vy - y*vx)*m)

    return lx,ly,lz


def gas_scatter_cbar(xx,yy,colors,fig,ax, cdx):
    cax = ax.scatter(xx, yy, s=9 * cdx**2, c=colors, edgecolor='none', alpha=0.2)
    fig.colorbar(cax, ax=ax)

def draw_gas_map(gal, region, info, axs, proj="y", method="scatter",
                 verbose=False, npix = 300):
    #xdx = width
    #ydx = height
    #zdx = depth

    #i = (abs(xx) < xdx) * (abs(yy) < ydx) * (abs(zz) < zdx)
    i = (gal.cell['x'] > region['ranges'][0][0]) *         (gal.cell['x'] < region['ranges'][0][1]) *         (gal.cell['y'] > region['ranges'][1][0]) *         (gal.cell['y'] < region['ranges'][1][1]) *         (gal.cell['z'] > region['ranges'][2][0]) *         (gal.cell['z'] < region['ranges'][2][1])

    if sum(i) < 100:
        return
    
    cdx = gal.cell['dx'][i]

    # scatter version
    if method == "scatter":
        if proj == "z":
            xx = gal.cell['x']
            yy = gal.cell['y']
            zz = gal.cell['z']
        if proj == "x":
            xx = gal.cell['y']
            yy = gal.cell['z']
            zz = gal.cell['x']
        if proj == "y":
            xx = gal.cell['z']
            yy = gal.cell['x']
            zz = gal.cell['y']
        
        #cz = zz[i]
        cell_ok = gal.cell[i]
        gas_scatter_cbar(xx[i], yy[i], np.log10(cell_ok["rho"]), fig, axs[0], cdx)
        gas_scatter_cbar(xx[i], yy[i], cell_ok["vx"], fig, axs[1], cdx)
        gas_scatter_cbar(xx[i], yy[i], cell_ok["vy"], fig, axs[2], cdx)
        gas_scatter_cbar(xx[i], yy[i], cell_ok["vz"], fig, axs[3], cdx)
        gas_scatter_cbar(xx[i], yy[i], np.log10(cell_ok["temp"]), fig, axs[4], cdx)
    elif method == "pp":
        cell_ok = gal.cell[i]
        region=None
        gas_map = pp_cell(cell_ok, npix, info, region=region, column=0, hvar="rho")
        cax = axs[0].imshow(np.log10(gas_map))
        cbar = fig.colorbar(cax, ax=axs[0])

        gas_map = pp_cell(cell_ok, npix, info, region=region,
                          column=0, hvar = "vx", verbose=verbose)
        cax = axs[1].imshow(gas_map)
        cbar = fig.colorbar(cax, ax=axs[1])

        gas_map = pp_cell(cell_ok, npix, info, region=region,
                          column=0, hvar = "vy", verbose=verbose)
        cax = axs[2].imshow(gas_map)
        cbar = fig.colorbar(cax, ax=axs[2])

        gas_map = pp_cell(cell_ok, npix, info, region=region,
                          column=0, hvar = "vz", verbose=verbose)
        cax = axs[3].imshow(gas_map)
        cbar = fig.colorbar(cax, ax=axs[3])

        gas_map = pp_cell(cell_ok, npix, info, region=region, column=0,
                         hvar = "temp", verbose=verbose)
        cax = axs[4].imshow(np.log10(gas_map))
        cbar = fig.colorbar(cax, ax=axs[4])

    
    axs[0].set_title("log density")
    axs[1].set_title("vx")
    axs[2].set_title("vy")
    axs[3].set_title("vz")
    axs[4].set_title("log temp")


# In[4]:

def ang_evol_gal(prg_this, nouts, t_merger, draw_method="scatter"):
    old_in, old_ex, new_sb, new_sf, all_stars,     m_oi, m_oe, m_sb, m_sf, m_all, n_all,     lambda_r, reff, gmass, nout_ok,    ang_gas = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

    oldid = gal_before.star['id']

    for i, nout in enumerate(nouts):
        print(nout)
        prg_now = prg_this[prg_this["nout"] == nout]
        haloid = prg_now["Orig_halo_id"].squeeze()
        afile = gal_dir_in + "gal_stars_" + str(nout).zfill(3) + "_" + str(haloid).zfill(7)
        info = load.info.Info(nout=nout)

        try:
            gal_after = gal_from_GM(prg_now, info = info, fname=afile, cell=True)
        except:
            continue

        gal_after.gasmass=dict(dense=0, highlevel=0, total=0,
                            radial=[], loc=[], rgal=0)# rgal = 4*reff?

        fig, axs = plt.subplots(2,3)
        axs = axs.ravel()
        fig.set_size_inches(14,6, forward=True)
        axs[0].hist2d(gal_after.star['x'], gal_after.star['y'], bins=[100,100])
        radi = gal_after.region['radius']
        region = smp.set_region(xr=[-3*radi, 3*radi], yr=[-3*radi,3*radi], zr=[-1 * radi, radi])
        if len(gal_after.cell['x']) > 500:                
            #gas_map = pp_cell(gal_after.cell, 100, info, region=region, column=1)#, field_var='rho')
            #gas_map = pp_cell(gal_after.cell, 100, info, region=region, column=1, field_var='vz')
            draw_gas_map(gal_after, region, info, axs[1:], method=draw_method)
        #gas_map = draw.pp.pp_cell(gal_before.cell[ind], 100, info, verbose=False)
            #ax[1].imshow(np.log10(gas_map),origin="lower")

        plt.savefig(gal_dir_out + "star_gasmap_" + str(nout).zfill(3) +".png")
        plt.close()

        # before and after merger
        #i_old = gal_after.star['time'] > t_merger
        #i_new = ~i_old

        newid = gal_after.star['id']

        # indices of old galaxy stars that are also found in the new galaxy.
        i_old_remaining = mtc.match_list_ind(oldid, newid, allow_swap=False)
        # indices of old galaxy stars that are not found from the new galaxy.
        i_old_gone = np.setdiff1d(np.arange(len(oldid)), i_old_remaining, assume_unique=True)

        # indices of new galaxy stars that are young
        i_new_young = np.where(gal_after.star['time'] < t_merger)[0]
        # indices of new galaxy stars that are also found in the old galaxy
        i_new_remaining = mtc.match_list_ind(newid, oldid, allow_swap=False)

        # old stars, but not in the old galaxy.
        i_new_old_ext = np.intersect1d(np.where(gal_after.star['time'] > t_merger)[0],                               
                                   np.setdiff1d(np.arange(len(newid)), i_new_remaining))    

        old_in.append(angular_momentum(gal_after.star[i_new_remaining]))
        old_ex.append(angular_momentum(gal_after.star[i_new_old_ext]))

        i_new_sb = i_new_young[np.where(gal_after.star['time'][i_new_young] > t_merger - 1)[0]]
        i_new_sf = i_new_young[np.where(gal_after.star['time'][i_new_young] < t_merger - 1)[0]]

        new_sb.append(angular_momentum(gal_after.star[i_new_sb]))
        new_sf.append(angular_momentum(gal_after.star[i_new_sf]))

        all_stars.append(angular_momentum(gal_after.star))

        m_oi.append(np.sum(gal_after.star['m'][i_new_remaining]))
        m_oe.append(np.sum(gal_after.star['m'][i_new_old_ext]))
        m_sb.append(np.sum(gal_after.star['m'][i_new_sb]))
        m_sf.append(np.sum(gal_after.star['m'][i_new_sf]))
        m_all.append(np.sum(gal_after.star['m']))
        n_all.append(len(gal_after.star['m']))

        ang_gas.append(angular_momentum_gas(gal_after.cell, info))

        gal_after.cal_lambda_r_eps(mge_interpol=False)
        gal_after.plot_gal(fn_save=gal_dir_out + str(nout) + '_' + str(haloid)                                    + "_" + str(prg_now['id'].squeeze()) +'.png')
        lambda_r.append(gal_after.meta.lambda_r)
        reff.append(gal_after.meta.reff)
        gas_mass(gal_after, info)
        gmass.append(gal_after.gasmass)
        nout_ok.append(nout)

    old_in = np.array(old_in)
    old_ex = np.array(old_ex)
    new_sf = np.array(new_sf)
    new_sb = np.array(new_sb)
    all_stars = np.array(all_stars)

    m_oi = np.array(m_oi)
    m_oe = np.array(m_oe)
    m_sb = np.array(m_sb)
    m_sf = np.array(m_sf)
    m_all = np.array(m_all)
    n_all = np.array(n_all)
    lambda_r = np.array(lambda_r)
    reff = np.array(reff)
    nouts = np.array(nout_ok)
    gmass = pd.DataFrame(gmass).to_records()
    ang_gas = np.array(ang_gas)
    
    return old_in, old_ex, new_sf, new_sb, all_stars,            m_oi, m_oe, m_sb, m_sf, m_all, n_all,            lambda_r, reff, nouts, gmass, ang_gas

        
def final_plot(old_in, old_ex, new_sf, new_sb, all_stars,     m_oi, m_oe, m_sb, m_sf, m_all, n_all,     lambda_r, reff, nouts, gmass, ang_gas, ax):

    ax[0][0].plot(nouts,old_in[:,0], label=r"$j_{x}$")
    ax[0][0].plot(nouts,old_in[:,1], label=r"$j_{y}$")
    ax[0][0].plot(nouts,old_in[:,2], label=r"$j_{z}$")
    leg = ax[0][0].legend()
    leg.get_frame().set_alpha(0.5)
    ax[0][0].set_ylabel("old original")

    ax[0][1].plot(nouts,old_ex[:,:])
    ax[0][1].set_ylabel("old accreted")

    ax[0][2].plot(nouts,new_sb[:,:])
    ax[0][2].set_ylabel("new starburst")

    ax[0][3].plot(nouts,new_sf[:,:])
    ax[0][3].set_ylabel("new smooth SF")

    ax[1][0].plot(nouts, all_stars[:,:])
    ax[1][0].set_ylabel("Total")

    ax[1][1].plot(nouts, m_all, label="all")
    ax[1][1].plot(nouts, m_oi, label="old in")
    ax[1][1].plot(nouts, m_oe, label="old ex")
    ax[1][1].plot(nouts, m_sf, label="new SF")
    ax[1][1].plot(nouts, m_sb, label="new SB")
    ax[1][1].set_ylabel(r"$M_{*}$")
    leg = ax[1][1].legend()
    leg.get_frame().set_alpha(0.5)
    #ax[1][1].set_zorder(0)

    # I am not sure if it is smooth accretion or 'numerically' delayed star burst..
    lns1 = ax[1][2].plot(nouts, lambda_r, 'b-', label="lambda_r")
    ax_reff = ax[1][2].twinx()
    lns2 =ax_reff.plot(nouts, reff, 'r-', label="reff")
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    leg = ax[1][2].legend(lns, labs, loc=0)
    leg.get_frame().set_alpha(0.5)
    ax[1][2].set_ylabel("Lambda r")


    ax[1][3].plot(nouts, gmass["rgal"]) # gas mass inside rgal
    #ax[1][3].plot(nouts, m_all, 'r-')
    #ax[1][3].plot(nouts, n_all, 'b-')

    ax[1][3].set_ylabel("Gas mass")

    ax[0][0].set_xlim([nout_merger, nout_fi + int((nout_fi - nout_merger) * 0.3) ])


# In[2]:

import numpy as np
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import pickle
import os
import load
import utils.cosmology
import utils.match as mtc
from galaxy import galaxy
from analysis.single_gal_module import *
#from analysis.single_gal_module import angular_momentum
import pandas as pd
import utils.sampling as smp
import draw
import tree.ctutils as ctu

cluster = "29176"
nout_fi = 187
prgt = pickle.load(open("./prg_only_tree.pickle", "rb"))
#galidx = 347991
# FR : 265877(189)
# nout_merger = None
#set1 = dict(galidx = 265840, galid = 152, nout_merger=173)
#set2 = dict(galidx = 265710, galid = 22, nout_merger=75)
#set3 = dict(galidx = 265877, galid = 189, nout_merger=None)

#gal_set = set2
# SR : 265710(22)
#galidx = gal_set['galidx']
#nout_merger = gal_set['nout_merger']


mpgs = pickle.load(open("./mpgs" + cluster + ".pickle", "rb"))

gal = mpgs[5]
galidx = gal.idxs[0]
galid = gal.ids[0]

nout_merger = 120
#nout_merger=gal.merger.nout[-1] -1

gal_dir_in = "./GalaxyMaker/gal_" + str(galidx) + "/"
gal_dir_out = "gal_" + str(galidx) + "/"

if not os.path.isdir(gal_dir_out):
    os.mkdir(gal_dir_out)

prg_this = ctu.extract_main_tree(prgt, idx=galidx)

# Galaxy before merger
nout_merger_org = nout_merger
#while True:
prg_now = prg_this[prg_this["nout"] == nout_merger_org]
#    if prg_now["phantom"] == 0:
#        break
#    nout_merger_org -=1

info = load.info.Info(nout=nout_merger_org)

haloid = prg_now["Orig_halo_id"].squeeze()
afile = gal_dir_in + "gal_stars_" + str(nout_merger_org).zfill(3) + "_" + str(haloid).zfill(7)

gal_before = gal_from_GM(prg_now, info = info, fname=afile, cell=True)
gal_before.gasmass=dict(dense=0, highlevel=0, total=0,
                        radial=[], loc=[], rgal=0) # rgal = 4*reff?
gas_mass(gal_before, info)


# In[66]:

results = []
for igal, gal in enumerate(mpgs[5:6]):
    if gal.merger is not None:
        # Galaxy ID
        print(igal, "-th, ", gal_dir_out)
        galidx, galid = gal.idxs[0], gal.ids[0]
        
        # dir
        gal_dir_in = "./GalaxyMaker/gal_" + str(galidx) + "/"
        gal_dir_out = "gal_" + str(galidx) + "/"

        if not os.path.isdir(gal_dir_out):
            os.mkdir(gal_dir_out)

        prg_this = ctu.extract_main_tree(prgt, idx=galidx)

        # Galaxy before merger
        nout_merger=gal.merger.nout[-1] -1
        nout_merger_org, prg_now = get_nout_merger_org(prg_this, nout_merger)
        
        # global data at this nout
        info = load.info.Info(nout=nout_merger_org)

        # file name
        haloid = prg_now["Orig_halo_id"].squeeze()
        afile = gal_dir_in + "gal_stars_" + str(nout_merger_org).zfill(3) +                 "_" + str(haloid).zfill(7)

        gal_before = gal_from_GM(prg_now, info = info, fname=afile, cell=True)
        gal_before.gasmass=dict(dense=0, highlevel=0, total=0,
                                radial=[], loc=[], rgal=0) # rgal = 4*reff?
        gas_mass(gal_before, info)
        if len(gal_before.cell) > 200:
            fig, ax = plt.subplots(1,2)
            # forward = True: change the size of the interactive window
            fig.set_size_inches(14,6, forward=True)
            radi = gal_before.region['radius']
            region = smp.set_region(xr=[-3*radi, 3*radi], yr=[-3*radi,3*radi], zr=[-1 * radi, radi])
            
            ax[0].hist2d(gal_before.star['x'], gal_before.star['y'], bins=[100,100])
            
            dxs = np.unique(gal_before.cell['dx'])
            cell = gal_before.cell[gal_before.cell['dx'] <= dxs[-3]]           
            gas_map = pp_cell(cell, 300, info, region=region, column=1)#, field_var='rho')
            
            ax[1].imshow(np.log10(gas_map), origin="lower")
            plt.savefig(gal_dir_out + "star_gasmap_" +                         str(nout_merger_org).zfill(3) +".png")
            plt.close()

            
        # Loop over subsequent nouts.
        #outs = range(nout_merger, nout_fi + 1, 1)
        nouts = [89,91]
        nnouts = len(nouts)

        t_merger = nout2lbt(nout_merger)        
        
#       old_in, old_ex, new_sf, new_sb, all_stars, \
#       m_oi, m_oe, m_sb, m_sf, m_all, n_all, \
#       lambda_r, reff, nouts, gmass, ang_gas = \
        all_outs = ang_evol_gal(prg_this, nouts, t_merger, draw_method="pp")
        
        
        # Final plot
        fig, ax = plt.subplots(4,2, sharex=True)
        fig.set_size_inches(16,12, forward=True)
        ax = ax.T
#       final_plot(old_in, old_ex, new_sf, new_sb, all_stars, \
#                   m_oi, m_oe, m_sb, m_sf, m_all, n_all, \
#                   lambda_r, reff, nouts, gmass, ang_gas, ax)
        final_plot(*all_outs, ax)

        fig.suptitle(galidx)

        plt.savefig(gal_dir_out + "plots.png")
        plt.close()
        results.append([all_outs])
pickle.dump(results, open("single_gal_results.picle", "wb"))


# In[56]:

plt.scatter(nouts, signed_log(ang_gas[:,0]), c='r', edgecolor="none")
plt.scatter(nouts, signed_log(ang_gas[:,1]), c='b', edgecolor="none")
plt.scatter(nouts, signed_log(ang_gas[:,2]), c='g', edgecolor="none")

plt.scatter(nouts, signed_log(all_stars[:,0]), c='r', marker="*", edgecolor="none")
plt.scatter(nouts, signed_log(all_stars[:,1]), c='b', marker="*", edgecolor="none")
plt.scatter(nouts, signed_log(all_stars[:,2]), c='g', marker="*", edgecolor="none")

plt.show()


# In[55]:

def signed_log(val):
    return np.sign(val) * np.log10(abs(val))


# In[12]:

gal_dir_out


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# Old version
