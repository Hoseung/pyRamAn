# coding: utf-8

import matplotlib.pyplot as plt

def save_result(savedir, galid, data):
    """    
    parameters
    ----------
    savedir : str
    
    galid : int
    
    """
    # savedir = 'simulation/galaxyXXXX/'
    header = "NOUT      ID        X        Y       Z[Mpc] \
    Vx      Vy     Vz[km/s]    Reff[kpc]   Mstar    Mgas[Msun]   Lambda_r"
    fmt = "  %3i  %5i  %.5f  %.5f  %.5f  %.4f  %.4f  %.4f  %.3f  %.f  %.f  %.6f"
    
    fname = savedir + "galaxy_" + str(galid).zfill(5) + '.txt'
    with open(fname, 'wb') as fout:
        np.savetxt(fout, [],  header=header) # Header
        for data in dictout:
            np.savetxt(fout, np.c_[data['nout'], data['id'],
                            data['xc'], data['yc'], data['zc'],
                            data['vx'], data['vy'], data['vz'],
                            data['rgal'], data['mstar'],
                            data['mgas'], data['lambda_r']],
                            fmt=fmt) # Append
    print("Save_result Done")

def load_gal_txt(savedir, galid):
    fname = savedir+ str(galid).zfill(5) + '.txt'    
    return np.loadtxt(fname, dtype = 'float, int', skiprows=0, unpack=True)    

def plot_result_pickle(fname, savedir):
    import pickle
    with open(fname, 'rb') as f:
        dd = pickle.load(f)
        # plot data
    fig, host = plt.subplots()
    fig.subplots_adjust(right=0.75)
    par1 = host.twinx()
    host.plot(nouts[::-1], dd['lambda_r'])    
    par1.plot(nouts[::-1], dd['mstar'], 'r-')
    par1.set_ylabel("Stellar mass [M_{\odot}]")
    #plt.show()
    plt.savefig(savedir +  + str(dd[0]).zfill(5) + '.png')
    plt.close()
    
    

def plot_result_txt(savedir, galid):   
    # Load data
    index, gal_id, xc, yc, zc, vx, vy,vz, rgal, mstar, mgas, lambda_r = \
                                                load_gal_txt(savedir, galid)

    # plot data
    fig, host = plt.subplots()
    fig.subplots_adjust(right=0.75)
    par1 = host.twinx()
    host.plot(nouts[::-1], lambda_r[i])    
    par1.plot(nouts[::-1], mstar[i], 'r-')
    par1.set_ylabel("Stellar mass [M_{\odot}]")
    #plt.show()
    plt.savefig(savedir +  + str(galid).zfill(5) + '.png')
    plt.close()


def test_galaxy(galaxy, npart_min=1000):
    # Is it a galaxy? 
    # minimum particles
    npart = len(galaxy.star.x)
    print("number of particles:", npart)
    if npart < npart_min:
        print("!! It has only {} particles. Not a galaxy.")
        return
    
    # Test coordinate origin


    # Test projection


    # Test population
        # Stellar particle

        # DarkMatter particle

            # Sink particle in DM?

        # Gas cell




    # Test unit
        # position

        # velocity

        # mass

    # Test for existances
        # normal vection

        # rotation matrix


#%%
def mk_gal(halodata, out_q, info, i, final_gal,
            save=False, rscale=0.3, verbose=False, galaxy_plot_dir='./',
            rscale_lambda=2.0, npix_lambda=50, npix=400):
    """
    Direct plot,
    Create galaxy, 
    Calculate lambda_r (using Cappellari 2003)
    Draw ed map of galaxy.
    
    """
    from galaxymodule import galaxy
    import utils.sampling as smp
    import draw
    import matplotlib.pyplot as plt    
#    print("IDs:", id(star), id(dm), id(cell))
    # 
    id_gal = halodata['id']
    sid_gal = str(id_gal).zfill(5)
    sid_fgal = str(final_gal).zfill(5)
    snout = str(info.nout).zfill(5)

    gal_out = {"id":0, "xc":0.0, "yc":0.0, "zc":0.0,
               "vx":0.0, "vy":0.0, "vz":0.0,
               "mstar":0.0, "nstar":0.0, "mgas":0.0,
               "lambda_arr":[], "lambda_r":0, "rgal":0,
               "final":final_gal, "mhal":0.0, "rhal":0.0}

    # Direct plot ---------------------------------------------------------
    extent = (0, npix, 0, npix)
    region = smp.set_region(xc=halodata['x'],
                            yc=halodata['y'],
                            zc=halodata['z'],
                            radius = halodata['rvir'])
    star_map = draw.pp.den2d(star['x'],
                             star['y'],
                             star['z'],
                             star['m'], 
                             npix,
                             region=region,
                             cic=True,
                             norm_integer=False)
    if star_map is not False:
        ls = np.zeros((npix,npix))
        ii = star_map > 0
        ls[ii] = np.log10(star_map[ii]) # Stellar map HAS empty pixels.
        ls[star_map <= 0] = np.floor(ls.min())
        im1 = plt.imshow(ls,
                         cmap="gray",
                         interpolation='nearest',
                         extent=extent)       
    
    # One of two should be transposed.
    # But which one?
    gas_map = draw.pp.pp_cell(cell, npix, info, region=region, verbose=False)
    im2 = plt.imshow(np.transpose(np.log10(gas_map)),
                     cmap="CMRmap",
                     alpha=.5,
                     interpolation='bilinear',
                     extent=extent)

    rgal = region['radius'] * s.info.pboxsize * 1000

    ax = plt.gca()
    ax.set_xlabel("position [kpc]")
    ax.set_xticks(np.linspace(0,npix,5))
    xticks = ["{:.2f}".format(x) for x in np.linspace(-rgal, rgal, num=5)]
    ax.set_xticklabels(xticks)
    ax.set_ylabel("position [kpc]")
    ax.set_yticks(np.linspace(0,npix,5))
    yticks = ["{:.2f}".format(y) for y in np.linspace(-rgal, rgal, num=5)]
    ax.set_yticklabels(yticks)

    fn_suffix = snout + "_" + sid_gal + "_" + sid_fgal +'.png'    
    plt.savefig(galaxy_plot_dir + "2dmap_" + fn_suffix, dpi=144)
    plt.close()
    #Create galaxy ---------------------------------------------------------
    gal = galaxy.Galaxy(halodata, radius_method='simple', info=info)
    is_gal = gal.mk_gal(star=star, dm=dm, cell=cell,
               rscale=rscale, verbose=verbose)
    #-----------------------------------------------------------------------    
    print(gal.id, "IS_GAL",is_gal)
    if not is_gal:
        print(gal.id, " Not a good galaxy")
        out_q.put(gal_out)
    else:
        # Save to catalog -----------------------------------------------------
        gal.cal_lambda_r(npix=npix_lambda, method=1, rscale=rscale_lambda)
        # Calculate lambda_r -------------------------------------------------
        gal.plot_gal(fn_save = galaxy_plot_dir + "galaxyplot" + fn_suffix,
                     ioff=True)
    #        gal.save_gal(base=wdir)

        # Instead of galaxy class, save them in a dict.    
        gal_out['mstar'] = gal.mstar
        gal_out['mgas'] = gal.mgas
        gal_out['nstar'] = gal.nstar
        gal_out['id'] = gal.id
        gal_out['xc'] = gal.xc * info.pboxsize
        gal_out['yc'] = gal.yc * info.pboxsize
        gal_out['zc'] = gal.zc * info.pboxsize
        gal_out['vx'] = gal.vxc * info.kms
        gal_out['vy'] = gal.vyc * info.kms
        gal_out['vz'] = gal.vzc * info.kms        
        gal_out['lambda_arr'] = gal.lambda_arr
        gal_out['lambda_r'] = gal.lambda_r
        gal_out['rgal'] = gal.reff# * info.pboxsize * 1000.0 # in kpc
        gal_out['mhal'] = halodata['mvir']
        gal_out['rhal'] = halodata['rvir']
        out_q.put(gal_out)
    print("mk_gal Done")
    return
    
def plot_lambda(catalog, i_early, i_late, i_bad, out_dir='./'):
    import matplotlib.pyplot as plt
    plt.ioff()
    f = plt.figure()
    ax = f.add_subplot(111)
    #for i, val in enumerate(lambdar_arr):
    for i in i_early:
        a = np.asarray(catalog['lambda_arr'][i])
        ax.plot(a, 'r-', alpha=0.5) # Red = Early
    for i in i_late:
        ax.plot(catalog['lambda_arr'][i], 'b-', alpha=0.3) # Red = Early
    
    #plt.xlabel() # in the unit of Reff
    ax.set_title(r"$\lambda _{R}$") 
    ax.set_ylabel(r"$\lambda _{R}$") 
    ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
    ax.set_xlim(right=9)
    ax.set_xticks([0, 4.5, 9])
    ax.set_xticklabels(["0", "0.5", "1"])
    plt.savefig(out_dir + "lambdar_disk.png")
    plt.close()    
    
#%%    
#    import pickle
#    with open(fcat, 'rb') as f:
#        dd = pickle.load(f)
        # plot data
def _get_data_dict(dictout):
    nouts = []
    lambda_r = []
    mstar =[]
    for dd in dictout:
        nouts.append(dd['nout'])
        lambda_r.append(dd['lambda_r'])
        mstar.append(dd['mstar'])
    return nouts, lambda_r, mstar

def _get_data_catalog(catalog):
    return catalog['nout'],catalog['lambda_r'], catalog['mstar']

def plot_growth_history(data, nouts, idgal=None):
    zreds=[]
    aexps=[]
    import load
    nnouts = len(nouts)
    for nout in nouts:
        info = load.info.Info(nout=nout, base=wdir, load=True)
        aexps.append(info.aexp)
        zreds.append(info.zred)
    aexps = np.array(aexps)
    zreds = np.array(zreds)
    #%%
    def aexp2zred(aexp):
        return [1.0/a - 1.0 for a in aexp]
    
    def zred2aexp(zred):
        return [1.0/(1.0 + z) for z in zred]
    
    def lbt2aexp(lts):
        import astropy.units as u
        from astropy.cosmology import WMAP7, z_at_value
        zreds = [z_at_value(WMAP7.lookback_time, ll * u.Gyr) for ll in lts]
        return [1.0/(1+z) for z in zreds]
    
    # For a given list of nouts, 
    # calculate a nice-looking set of zreds.
    # AND lookback times
    z_targets=[0, 0.2, 0.5, 1, 2, 3]
    z_target_str=["{:.2f}".format(z) for z in z_targets]
    a_targets_z = zred2aexp(z_targets)
    z_pos =  [nout_ini + (1 - (max(aexps) - a)/aexps.ptp()) * nnouts for a in a_targets_z]
    
    lbt_targets=[0.00001,1,3,5,8,12]
    lbt_target_str=["{:.0f}".format(l) for l in lbt_targets]
    a_targets_lbt = lbt2aexp(lbt_targets)
    lbt_pos = [nout_ini + (1 - (max(aexps) - a)/aexps.ptp()) * nnouts for a in a_targets_lbt]


    nouts, lambda_r, mstar = _get_data_dict(data)
    fig, ax1 = plt.subplots()
    fig.subplots_adjust(right=0.75)
    fig.suptitle("ID: " + str(idgal).zfill(5), fontsize=18)#, y=1.01)
    fig.subplots_adjust(right=0.75)
    # label needed for legend()
    lns1 = ax1.plot(nouts, lambda_r, label=r"$\lambda_{R}$")
    ax1.set_xticks(z_pos)
    ax1.set_xticklabels(z_target_str)
    ax1.set_xlim([0,187])
    ax1.set_ylim([0,1.0])
    ax1.set_ylabel(r"$\lambda_{R}$")
    ax1.set_xlabel("redshift")

    ax2 = ax1.twinx()   
    lns2 = ax2.plot(nouts, mstar, 'r-', label="stellar mass")
    ax2.set_ylim([0, 1.3*max(mstar)])
    ax2.set_ylabel(r"Stellar mass $[M_{\odot}]$")

    ax3 = ax1.twiny()
    ax3.set_xlabel("Lookback time", labelpad=10)
    ax3.set_xticks(lbt_pos)
    ax3.set_xticklabels(lbt_target_str)            
    # Because there are two axes superimposed, 
    # manually gather labels of objects first.
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=0)
   
    plt.savefig(dir_gal + str(sidgal_final).zfill(5) + '.png')
    #plt.show()
    plt.close()       
#%%

"""
The processing pool needs to be instantiated in the main 
thread of execution. 
"""
import multiprocessing as mp    
import load
from tree import tmtree
import numpy as np
import utils.sampling as smp
#import ctypes
import tree.halomodule as hmo 
import os

if __name__ == '__main__':
#    ncore = int(input("How many cores? \n"))
#    wdir = input("Working directory \n")
    #nout = int(input("nout? \n"))
    #wdir = './'
    wdir = '/home/hoseung/Work/data/05427/'
    ncore = 1
    #ncore = 16
    
    # 27 : z=4;  37 : z=3;  20 : ~ z=5
    nout_ini = 186
    nout_fi = 187
    
    # 05427 galaxies:
    # [1600, 1612, 1645, 1648, 1664, 1665, 1669, 1681, 1686]
    gal_final = int(input("Galaxy ID"))
    
    mstar_min_plot = 1e9
    rscale = 0.8
    r_cluster_scale = 2.0 # maximum radius inside which galaxies are searched for
    npix=800
    rscale_lambda = 3.0
    npix_lambda = int(15 * rscale_lambda)
    lmax = 19
    ptypes=["star id pos mass vel time metal", "dm id pos mass vel"]
    frefine= 'refine_params.txt'
    fnml = 'cosmo_200.nml'
    
    ## halo part ###
    dir_out = wdir + 'catalog/'
       
    # nout_halo = 122 == nout 10, nout_halo = 0   == nout 132
    nouts = range(nout_fi, nout_ini -1, -1) 
    Nnouts = len(nouts)
    
    tt = tmtree.load(work_dir=wdir, filename="halo/TMtree.fits")
    tfin = tt[np.where(tt['NOUT'] == 0)]
    tini = tt[np.where(tt['NOUT'] == nout_fi - nout_ini)]
    #%%
    info = load.info.Info(nout=nout_fi, base=wdir, load=True)

    hh = hmo.Halo(base=wdir, nout=nout_fi,
                  halofinder='HM', info=info, load=True)

    i_center = np.where(hh.data['np'] == max(hh.data['np']))[0]
    i_satellites = smp.extract_halos_within(hh.data,
                                            i_center,
                                            scale=r_cluster_scale)
    print("Total {0} halos \n{1} halos are selected".format(
          len(i_satellites),sum(i_satellites)))
    #%%
    # ALL halos inside the cluster and have tree back to nout_ini
    halo_list = hh.data['id'][i_satellites]
    h_ind_ok, halo_ok = tmtree.check_tree_complete(tt, 0,
                                                   nout_fi - nout_ini,
                                                   halo_list)
    print(len(halo_ok), "halos left")
    
    igal = np.where(halo_ok[:,0] == gal_final)[0]
#    gal_final = halo_ok[igal,0]
    sidgal_final = str(gal_final).zfill(5)
    dir_gal = wdir + 'galaxy_' + str(gal_final).zfill(5) + '/'

    fcat = dir_gal +"catalog_" + sidgal_final + ".pickle"    
    
    if os.path.isdir(dir_gal) is False:
        os.mkdir(dir_gal)

    m = mp.Manager()
    out_q = m.Queue()
    dictout=[]

    for inout, nout in enumerate(nouts):
        print(inout, nout)
        snout = str(nout)
        
        info = load.info.Info(nout=nout, base=wdir, load=True)
               
        hh = hmo.Halo(base=wdir, nout=nout,
                      halofinder='HM', info=info, load=True)
  
        hind = np.where(hh.data['id'] == halo_ok[igal,inout])

        region = smp.set_region_multi(xc=hh.data['x'][hind],
                                      yc=hh.data['y'][hind],
                                      zc=hh.data['z'][hind],
                                      radius = hh.data['rvir'][hind] * rscale)

        gal_id_now = hh.data['id'][hind]

        s = load.sim.Sim()
        s.setup(nout, wdir)    
        s.set_ranges(region["ranges"])
        s.show_cpus()
        s.add_part(ptypes)
        s.part.load(fortran=True)
        s.add_hydro()
        s.hydro.amr2cell(lmax=19)

        star = s.part.star
        dm = s.part.dm
        cell = s.hydro.cell
    
        nh = len(hind)
    
        mk_gal(hh.data[hind][0], out_q, s.info, 1, gal_final,
               galaxy_plot_dir=dir_gal,
               verbose=False,
               rscale_lambda=rscale_lambda,
               npix_lambda=npix_lambda)

        print("----------Done---------")

        dd =  out_q.get(timeout=1)
        dd['nout'] = nout
        dictout.append(dd)

    save_result(dir_gal, gal_final, dictout)
    import pandas as pd
    catalog = pd.DataFrame(dictout).to_records()    

    import pickle
    with open(fcat, 'wb') as f:
        pickle.dump(catalog, f)  
    
    plot_growth_history(catalog, nouts, idgal=gal_final)
    
#    plot_merger_tree