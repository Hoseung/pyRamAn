#import utils.cosmology
from .. import load
import numpy as np

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
        file_list.append(gal_dir + "gal_stars_" + str(nout).zfill(3) \
                         + "_" + str(galid).zfill(7))

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


def draw_plot(gal, cat, ax, ax_hist, smoothed_lambda=True):
    #ax = ax.ravel()
    t_min = min(gal.star['time'])
    t_max = max(gal.star['time'])
    bins = np.linspace(t_min, t_max, num=4)
    i_bins = np.digitize(gal.star['time'],bins)
    ax[0].plot(nout2lbt(cat.nouts), cat.data['mstar'])
    ax_hist.hist(gal.star['time'], histtype='step')
    if smoothed_lambda:
        ax[2].plot(nout2lbt(cat.nouts), cat.smoothed_lambda)
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
    import tree.ctutils as ctu
    galidx = prgt["id"][(prgt["nout"] == 187) * (prgt["Orig_halo_id"] \
             == galid_at_187)][0]
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
    from ..utils.cosmology import Timeconvert
    tc = Timeconvert(s.info)
    from galaxymodule import galaxy
    if fname is None:
        gm = load.rd_GM.rd_gal(nout, idgal)
    else:
        gm = load.rd_GM.rd_gal(info.nout,0,fname=fname)
        galid = gm.header['my_number']

    if gm.star == None:
        print("no stellar particles, skipping")
        return False

    gg.star['time'] = tc.time2gyr(gg.star['time'],
                                    z_now = gg.info.zred)
    #replace_field_name(gg.star, "time", "age")
    gcat['x'] /= info.cboxsize
    gcat['y'] /= info.cboxsize
    gcat['z'] /= info.cboxsize
    #print("Cat x", gcat['x'])
    gal = galaxy.Galaxy(halo = gcat, info=info)
    
    if dm:
        fname_dm = fname.replace("stars", "dms")
        gm.dm = load.rd_GM.rd_dm(info.nout, 0, fname = fname_dm)
    else:
        gm.dm = None
    
    if cell == "GM":
        fname_cell = fname.replace("stars", "cells")
        gm.cell = load.rd_GM.rd_cell(info.nout, 0, fname = fname_cell)
        # unit of "rho" is not touched.
    elif cell == "raw":
        from ..load.hydro import Hydro
        import utils.sampling as smp
        centers = gm.header['xg'] /info.pboxsize + 0.5
         
        radius = 0.5 * max([gm.star['x'].ptp(), gm.star['y'].ptp(), \
                           gm.star['z'].ptp()])/info.pboxsize

        region = smp.set_region(centers=centers, radius=1.5*radius)
        hh = Hydro(info, region=region, load=True)
        gm.cell = hh.cell
        gm.cell["x"] = (gm.cell["x"] + 0.5) * info.pboxsize
        gm.cell["y"] = (gm.cell["y"] + 0.5) * info.pboxsize
        gm.cell["z"] = (gm.cell["z"] + 0.5) * info.pboxsize
        gm.cell['dx'] *= info.pboxsize

        
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
    if len(cell) < 100:
        return 0
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


