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

    if gm.star == None:
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

	
