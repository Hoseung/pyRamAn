# coding: utf-8

#%%
def mk_gal(halodata, out_q, info, i, final_gal,
            save=False, rscale=0.3, verbose=False, galaxy_plot_dir='./',
            rscale_lambda=2.0, npix_lambda=50, npix=400, galaxy_plot=False):
    """
    Direct plot,
    Create galaxy, 
    Calculate lambda_r (using Cappellari 2003)
    Draw ed map of galaxy.
    
    """
    from galaxy import galaxy
    import utils.sampling as smp
    import draw
    import matplotlib.pyplot as plt    
#    print("IDs:", id(star), id(dm), id(cell))

    gal_out = {"id":0, "xc":0.0, "yc":0.0, "zc":0.0,
               "vx":0.0, "vy":0.0, "vz":0.0,
               "mstar":0.0, "nstar":0.0, "mgas":0.0,
               "lambda_arr":[], "lambda_r":0, "rgal":0, "final":final_gal,
                "rhalo":halodata['rvir'], "boxtokpc":info.pboxsize*1000}
    # Direct plot ---------------------------------------------------------
    region = smp.set_region(xc=halodata['x'],
                            yc=halodata['y'],
                            zc=halodata['z'],
                            radius = halodata['rvir'])
    if galaxy_plot:
        extent = (0, npix, 0, npix)        
        star_map = draw.pp.den2d(star['x'],star['y'],star['z'],star['m'], npix,
                                 region=region, cic=True, norm_integer=False)
        if star_map is not False:
            ls = np.zeros((npix,npix))
            ii = star_map > 0
            ls[ii] = np.log10(star_map[ii]) # Stellar map HAS empty pixels.
            ls[star_map <= 0] = np.floor(ls.min())
            im1 = plt.imshow(ls, cmap="gray", interpolation='nearest', extent=extent)       
        
        # One of two should be transposed.
        # But which one?
#        gas_map = draw.pp.pp_cell(cell, npix, info, region=region, verbose=False)
#        im2 = plt.imshow(np.transpose(np.log10(gas_map)), cmap="CMRmap", alpha=.5, interpolation='bilinear', extent=extent)
    
        rgal = region['radius'] * s.info.pboxsize * 1000

        ax = plt.gca()
        ax.set_xlabel("position [kpc]")
        ax.set_xticks(np.linspace(0,npix,5))
        xticks = ["{:.2f}".format(x) \
                    for x in np.linspace(-rgal, rgal, num=5)]
        ax.set_xticklabels(xticks)
        ax.set_ylabel("position [kpc]")
        ax.set_yticks(np.linspace(0,npix,5))
        yticks = ["{:.2f}".format(y) \
                    for y in np.linspace(-rgal, rgal, num=5)]
        ax.set_yticklabels(yticks)
        
        plt.savefig(galaxy_plot_dir+"2dmap_"+str(halodata['id']).zfill(5)+'.png', dpi=144)
        plt.close()

    #-----------------------------------------------------------------------    
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
        # Save to catalog -------------------------------------------------------
#        print("Good galaxy, R eff:", gal.reff)

        gal.cal_lambda_r(npix=npix_lambda, method=1, rscale=rscale_lambda) # calculate within 1.0 * reff    
        # Calculate lambda_r ---------------------------------------------------

        gal.plot_gal(save_dir= galaxy_plot_dir, ioff=True)
    #            gal.save_gal(base=wdir)

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
        out_q.put(gal_out)
    print("mk_gal done")
    return
    

def plot_lambda(catalog, i_early, i_late, i_bad, fn_out='./'):
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
    plt.savefig(fn_out)
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
import tree.halomodule as hmo 
from utils import match
import time
import os

ncore = int(input("How man cores? \n"))
wdir = input("Working directory \n")
#nout = int(input("nout? \n"))
#wdir = './'
#wdir = '/home/hoseung/Work/data/05427/'
#ncore = 1
#ncore = 16

#----------------------------------------------------------------------
# 27 : z=4;  37 : z=3;  20 : ~ z=5
nout_ini = 37
nout_fi = 187
nouts = range(nout_fi, nout_ini -1, -1) 
#----------------------------------------------------------------------

mstar_min_plot = 1e9
mk_gal_rscale = 1.0 # unit of Rvir,halo
rscale = 1.5
r_cluster_scale = 2.0 # maximum radius inside which galaxies are searched for
npix=800
rscale_lambda = 2.0 # Reff unit.
npix_lambda = int(10 * rscale_lambda)
lmax = 19
ptypes=["star id pos mass vel time metal", "dm id pos mass vel"]

## halo part -----------------------------------------------------------
m_halo_min = 1e10 # minimum halo mass above which galaxies are searched for. 
dir_out = wdir + 'catalog/'

# optional parameters ----------------------------------------------------
lambda_plot = False 


# nout_halo = 122 == nout 10, nout_halo = 0   == nout 132
tt = tmtree.load(work_dir=wdir, filename="halo/TMtree.fits")
tfin = tt[np.where(tt['NOUT'] == 0)]
tini = tt[np.where(tt['NOUT'] == nout_fi - nout_ini)]


#%%
info = load.info.Info(nout=nout_fi, base=wdir, load=True)
hh = hmo.Halo(base=wdir, nout=nout_fi, halofinder='HM', info=info, load=True)

#halo = hmu.load_data(nout_fi, work_dir=work_dir, normalize=True)
i_center = np.where(hh.data['np'] == max(hh.data['np']))
i_satellites = smp.extract_halos_within(hh.data, i_center, scale=r_cluster_scale)
print("Total {0} halos \n{1} halos are selected".format(
      len(i_satellites),sum(i_satellites)))
#%%
# halos found inside the cluster and has tree back to nout_ini
halo_list = hh.data['id'][i_satellites]
h_ind_ok, halo_ok = tmtree.check_tree_complete(tt, 0, nout_fi - nout_ini, halo_list)
print(len(halo_ok), "halos left")


#%%------------------------------------------------------------
# r_mod array
nnouts = halo_ok.shape[1]
rarr = np.zeros(halo_ok.shape)
for j in range(nnouts):
    i_now = np.where(tt['NOUT'] == j)
    ttmp = tt[i_now]
    for i in range(halo_ok.shape[0]):
        i_hal = ttmp['HALNUM'] == halo_ok[i,j]
        rarr[i,j] = ttmp["RVIR"][i_hal]
    
r_mod = np.zeros(halo_ok.shape)
for i in range(nnouts):
    for j in range(1,halo_ok.shape[1]):
        rhal_min_possible = max([rarr[i,j], 0.9 * rarr[i,j-1]])
        r_mod[i,j] = rhal_min_possible/rarr[i,j]

#%%
import matplotlib.pyplot as plt
for i in r_mod:
    plt.plot(i)
    
#%%
#igal = [25] # must be a sequence 
igal = np.where(halo_ok[:,0] == 1442)[0]

#%%
#import galaxy
#import utils.util
#utils.util.reimport(galaxy.galaxy)

final_gal = halo_ok[:,0]
#for i in final_gal:
#    galaxy_plot_dir = wdir + 'galaxy_' + str(final_gal).zfill(5) + '/'    
    #    galaxy_plot_dir = wdir + 'galaxy_plot' + snout + '/'
#    if os.path.isdir(galaxy_plot_dir) is False:
#        os.mkdir(galaxy_plot_dir)

#%%
for inout, nout in enumerate(nouts):
    print(inout, nout)
    snout = str(nout)

    fcat = dir_out +"catalog" + snout + ".pickle"    
    info = load.info.Info(nout=nout, base=wdir, load=True)

    galaxy_plot_dir = wdir + 'galaxy_plot' + snout + '/'
    if os.path.isdir(galaxy_plot_dir) is False:
        os.mkdir(galaxy_plot_dir)

    # Load all halo
    hh = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info, load=True)
    hind = match.match_list_ind(hh.data['id'], halo_ok[:,inout])   
    h = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info)
    h.derive_from(hh, hind)
#    region = smp.set_region_multi(h.data.x, yc=h.data.y, zc=h.data.z,
#                                  radius = h.data.rvir * rscale)

#    print(h.data['id'])
    s = load.sim.Sim(nout=nout, base=wdir, setup=True)#,ranges=region["ranges"])

    s.show_cpus()
    t0 = time.time()    
    s.add_part(ptypes, load=True, fortran=True)
    assert s.part.nstar > 0, "Not enough stellar particles in given cpus"
        
    t1 = time.time()
    s.add_hydro(load=True, lmax=lmax)
    t2 = time.time()

    print("Loading particle took {}, \n and loading hydro took {}".format(t1-t0, t2-t1))    
    star = s.part.star
    dm = s.part.dm
    cell = s.hydro.cell

#   Multiprocessing -----------------------------------------------------------
    m = mp.Manager()
    out_q = m.Queue()
    nh = len(h.data)
    t3 = time.time()

    keywords = dict(galaxy_plot_dir=galaxy_plot_dir,
                    rscale = mk_gal_rscale,
                    verbose=True, rscale_lambda=rscale_lambda,
                    npix_lambda=npix_lambda, galaxy_plot = True)

    print("Looking for galaxies inside {} halos".format(nh))
    pool = mp.Pool(processes=ncore)
    for i in range(nh):
        pool.apply_async(mk_gal, args=(h.data[i], out_q,
                                       s.info, i, final_gal[i]), kwds=keywords)
    # Pass optional parameters in dict!s
    pool.close()
    pool.join()
    
    print("----------Done---------")
    print(" Took", time.time() - t3, "Seconds")
    
    dictout=[]
    try:
        if not os.path.isdir(dir_out):
            os.mkdir(dir_out)
        f = open(dir_out + 'galaxies' + snout + '.txt', 'w')
    except:
        print("No filename is given.\n ")
    
    f.write(" Final ID:" + str(final_gal) + "\n")
    f.write(" nout     ID        x          y       z[Mpc]       vx      vy     vz[km/s]")
    f.write("    Reff[kpc]     Mstar    Mgas[Msun]  Rhalo[kpc]  boxtokpc  \n")    
    
    for i in range(nh):
        print("dict out", i)
        try:
            tmp =  out_q.get(timeout=1)
            if tmp['id'] == 0:
                continue
            dd = tmp
            f.write("{:<4}   {:<4}   {:.5f}  {:.5f}  {:.5f}".format(i,
                    dd['id'],dd['xc'],dd['yc'],dd['zc']))
            f.write("  {:.3f}  {:.3f}  {:.3f}".format(dd['vx'],dd['vy'],dd['vz']))
            f.write("  {:.6f}  {:.0f}  {:.0f}".format(dd['rgal'],dd['mstar'], dd['mgas']))
            f.write("  {:.5f}  {:.5f}    \n".format(dd['rhalo'],dd['boxtokpc']))
            dictout.append(dd)
        except:
            continue
    
    f.close()    
    print("Text file written")
    import gc
    gc.collect()
    import pandas as pd
    catalog = pd.DataFrame(dictout).to_records()    
    # pickling also works!
    # The problem was memory being locked!
    import pickle
    with open(fcat, 'wb') as f:
        pickle.dump(catalog, f)  
    
    i_early = np.where(catalog['mstar'] > mstar_min_plot)[0]
    i_late = []
    i_bad = np.where(catalog.id == 0)[0]
    if lambda_plot:
    #    if not os.path.isdir(wdir + snout + '/'):
    #        os.mkdir(wdir + snout + '/')
        plot_lambda(catalog, i_early, i_late, i_bad,
                    fn_out = snout + "_lambdar_disk.png")

