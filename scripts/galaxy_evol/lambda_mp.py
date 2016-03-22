# coding: utf-8

#%%
import time

def extract_data(halo, rscale=0.25):
    xc_tmp0 = halo['x']
    yc_tmp0 = halo['y']
    zc_tmp0 = halo['z']
    
    rr_tmp0 = min([halo['rvir'] * rscale, 0.0002]) 
    # arbitrary! < 20kpc
    rr_tmp0 = max([rr_tmp0, 0.000025])
    # When merger occurs, larger radius is likely to include 
    # companion galaxy resulting center to be in the middle of nowhere.
    # If you want a larger galaxy, # increase rgal_tmp instead. 
    #        
    # xx is easier to search for than x.

    if star_all is not None:
        ind_s = np.where((star_all['x'] - xc_tmp0)**2 + (star_all['y'] - yc_tmp0)**2 
                        + (star_all['z'] - zc_tmp0)**2 < rr_tmp0**2)[0]
    if dm_all is not None:
        ind_d = np.where((dm_all['x'] - xc_tmp0)**2 + (dm_all['y'] - yc_tmp0)**2 
                        + (dm_all['z'] - zc_tmp0)**2 < rr_tmp0**2)[0]
    if cell_all is not None:
        ind_c = np.where((cell_all['x'] - xc_tmp0)**2 + (cell_all['y'] - yc_tmp0)**2 
                        + (cell_all['z'] - zc_tmp0)**2 < rr_tmp0**2)[0]
    else:
        return star_all[ind_s], dm_all[ind_d], None
        
#    print(len(ind_s), len(ind_d), len(ind_c))    

    return star_all[ind_s], dm_all[ind_d], cell_all[ind_c]

def mk_gal(halodata, out_q, info, i, final_gal,
           save=False, rscale=0.3, verbose=False, galaxy_plot_dir='./',
           rscale_lambda=2.0, npix_lambda=50, npix=400, galaxy_plot=False,
           method_com=2, mstar_min=5e9):
    """
    Direct plot,
    Create galaxy, 
    Calculate lambda_r (using Cappellari 2003)
    Draw ed map of galaxy.
    
    """
    #t = time.time()
    #print(i, time.time() - t, "seconds --- 1")
    
#    print(" !!!!!!!!!!" )
    print("This is {}-th halo".format(i))
#    print(" !!!!!!!!!! \n" )
    from galaxy import galaxy

    #print("IDs:", id(star), id(dm), id(cell))

    gal_out = {"id":0, "xc":0.0, "yc":0.0, "zc":0.0,
               "vx":0.0, "vy":0.0, "vz":0.0,
               "mstar":0.0, "nstar":0.0, "mgas":0.0,
               "lambda_arr":[], "lambda_r":0, "rgal":0, "final_gal":final_gal,
               "rhalo":halodata['rvir'], "boxtokpc":info.pboxsize*1000}
               
    star, dm, cell = extract_data(h.data[i])

    if sum(star['m']) * info.msun < mstar_min:
        print("(1)Not enough stars: {:.2f} Msun".format(sum(star['m']) * info.msun))
        print("Aborting... \n")
        print(" Not a good galaxy")
        out_q.put(gal_out)
        return
               
               
    # Direct plot ---------------------------------------------------------                                
    if galaxy_plot:
        import utils.sampling as smp
        import draw
        import matplotlib.pyplot as plt            
        region = smp.set_region(xc=halodata['x'],
                            yc=halodata['y'],
                            zc=halodata['z'],
                            radius = halodata['rvir'])

        extent = (0, npix, 0, npix)        
        star_map = draw.pp.den2d(star['x'],star['y'],star['z'],star['m'], npix,
                                 region=region, cic=True, norm_integer=False)
        if star_map is not False:
            ls = np.zeros((npix,npix))
            ii = star_map > 0
            ls[ii] = np.log10(star_map[ii]) # Stellar map HAS empty pixels.
            ls[star_map <= 0] = np.floor(ls[ii].min())
            plt.imshow(ls, cmap="CMRmap", interpolation='nearest', extent=extent)
        
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

    #Create galaxy ----------------------------------------------
    gal = galaxy.Galaxy(halodata, radius_method='eff', info=info)
    #print(i, time.time() - t, "seconds ---2")
    is_gal = gal.mk_gal(star, dm, cell,
                        mstar_min=mstar_min,
               rscale=rscale, verbose=verbose, method_com=method_com)
    #print(i, time.time() - t, "seconds ---3")               
    #-----------------------------------------------------------------------    
#    print(gal.id, "IS_GAL",is_gal)
    if not is_gal:
        print(gal.id, " Not a good galaxy")
        out_q.put(gal_out)
    else:
        # Save to catalog -------------------------------------------------------
#        print("Good galaxy, R eff:", gal.reff)
        gal.cal_lambda_r(npix=npix_lambda, method=1, rscale=rscale_lambda) # calculate within 1.0 * reff    
        #print(i, time.time() - t, "seconds ---4")
        # Calculate lambda_r ---------------------------------------------------

#        gal.plot_gal(fn_save = galaxy_plot_dir + str(nout).zfill(3) \
#                             + "_" + str(final_gal).zfill(5) + "_"  \
#                             + str(gal.id) + ".png", ioff=True)
#       gal.save_gal(base=wdir)

        # save in a dict.
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

#    print("mk_gal done \n")
    

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
    
  

def set_affinity_on_worker():
    import os
    """When a new worker process is created, the affinity is set to all CPUs"""
    print("I'm the process %d, setting affinity to all CPUs." % os.getpid())
    os.system("taskset -p 0xff %d" % os.getpid())    
    
    

def worker(halodata, out_q, info, inds, final_gal, **kwargs):        
    worker_q = Queue()
    if type(inds) == int:
        inds = [inds]
    for i in inds:
        mk_gal(h.data[i], worker_q, s.info, i, final_gal[i], **kwargs)
    for i in inds:
        out_q.put(worker_q.get())


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
import os
multi = True # 
hydro = False

#wdir = input("Working directory \n")
wdir = './'
#wdir = '/home/hoseung/Work/data/05427/'
#ncore = 16
if multi: ncore = int(input("How many cores? \n"))
nout_ini = input("First snapshot: (default = 37 , z=3) \n")
nout_end = input("Last snapshot: (default = 187, z=0) \n")
#ncore=1
#nout_ini=186
#nout_end=186
#----------------------------------------------------------------------
# 27 : z=4;  37 : z=3;  20 : ~ z=5

if nout_ini == "":
    nout_ini = 37
else:
    nout_ini = int(nout_ini)

if nout_end == "":
    nout_end = 187
else:
    nout_end = int(nout_end)

nout_ini0 = 37
nout_fi = 187
nouts = range(nout_fi, nout_ini -1, -1) 
#----------------------------------------------------------------------
mstar_min = 5e9
# Only galaxies above this stellar mass at the final snapshot are considered.
mstar_min_plot = 5e9
mk_gal_rscale = 1.0 # unit of Rvir,halo
rscale = 1.5
r_cluster_scale = 2.0 # maximum radius inside which galaxies are searched for
npix=800
rscale_lambda = 2.0 # Reff unit.
npix_lambda = int(10 * rscale_lambda)
lmax = 19
ptypes=["star id pos mass vel time metal", "dm id pos mass vel"]

## halo part -----------------------------------------------------------
m_halo_min = 2e10 # minimum halo mass above which galaxies are searched for. 
dir_out = wdir + 'catalog/'

# optional parameters ----------------------------------------------------
lambda_plot = False 

# nout_halo = 122 == nout 10, nout_halo = 0   == nout 132
tt = tmtree.load(work_dir=wdir, filename="halo/TMtree.fits")
tfin = tt[np.where(tt['NOUT'] == 0)]
tini = tt[np.where(tt['NOUT'] == nout_fi - nout_ini0)]


info = load.info.Info(nout=nout_fi, base=wdir, load=True)

#import galaxy
#import utils.util
#utils.util.reimport(galaxy.galaxy)
if nout_end == nout_fi:
    hh = hmo.Halo(base=wdir, nout=nout_fi, halofinder='HM', info=info, load=True)
    i_center = np.where(hh.data['np'] == max(hh.data['np']))
    i_satellites = smp.extract_halos_within(hh.data, i_center, scale=r_cluster_scale)
    print("Total {0} halos \n{1} halos are selected".format(
          len(i_satellites),sum(i_satellites)))
    
    # halos found inside the cluster and has tree back to nout_ini
    large_enugh = hh.data['mvir'] > m_halo_min
    halo_list = hh.data['id'][i_satellites * large_enugh]
    h_ind_ok, halo_ok = tmtree.check_tree_complete(tt, 0, nout_fi - nout_ini0, halo_list)
    print(len(halo_ok), "halos left")
    final_gal = halo_ok[:,0]
    ngals = len(final_gal)


from queue import Queue
import pandas as pd
import pickle
import time
#os.system("taskset -p 0xff %d" % os.getpid())
# http://stackoverflow.com/questions/15414027/multiprocessing-pool-makes-numpy-matrix-multiplication-slower
#%%
t0 =  time.time()


for inout, nout in enumerate(nouts):
    print(inout, nout, nout_end)
    print("Mstar min:", mstar_min, )
    if nout > nout_end:
        mstar_min = 0
        continue

    # If nout != 187,
    # load the final galaxy list from galaies187.txt
    # and the list remains for the rest of loop
    if (nout == nout_end) and (nout < nout_fi):
        from astropy.io import ascii
        data = ascii.read(dir_out + 'galaxies187.txt')
        halo_list = data['final_ID']
        mstar_min = 0
    
        # hlo_okk = complete-tree halos at all nouts.
        h_ind_ok, halo_ok = tmtree.check_tree_complete(tt, 0, nout_fi - nout_ini0, halo_list)
        
        print(len(halo_ok), "halos left")
        final_gal = halo_ok[:,0]
        ngals = len(final_gal)


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
    print("Final target halos:", h.data['id'], len(h.data['id']))

    s = load.sim.Sim(nout=nout, base=wdir, setup=True)#,ranges=[[0.4,0.5],[0.4,0.5],[0.4,0.5]])
    s.add_part(ptypes, load=True, fortran=True)
    #assert s.part.nstar > 0, "Not enough stellar particles in given cpus"
    if hydro:
        s.add_hydro(load=True, lmax=lmax)
        cell_all = s.hydro.cell
    else:
        cell_all = None

    star_all = s.part.star
    dm_all = s.part.dm
    

    nh = len(h.data)
    keywords = dict(galaxy_plot_dir=galaxy_plot_dir,
                rscale = mk_gal_rscale,
                verbose=False, rscale_lambda=rscale_lambda,
                npix_lambda=npix_lambda, galaxy_plot = False,
                mstar_min=mstar_min)

    if multi == 1:
#   Multiprocessing -----------------------------------------------------------
        m = mp.Manager()
        out_q = m.Queue()
        print("Looking for galaxies inside {} halos".format(nh))
        inds=[]
        [inds.append([]) for i in range(ncore)]
        
        for i in range(ngals):
            j = i % ncore
            inds[j].append(i)

        print(inds)
        processes = [mp.Process(target=worker, args=(h.data, out_q,
                    s.info, inds[i], final_gal), kwargs=keywords) for i in range(ncore)]

        # Run processes
        for p in processes:
            p.start()
        
        # Exit the completed processes
        for p in processes:
            p.join()
            
    elif multi == 2:
        m = mp.Manager()
        out_q = m.Queue()
        print("Looking for galaxies inside {} halos".format(nh))
        
        pool = mp.Pool(processes=ncore)
        for i in range(nh):
            pool.apply_async(mk_gal, args=(h.data[i], out_q,
                              s.info, i, final_gal[i]), kwds=keywords)
        pool.close()
        pool.join()
    else:
        for i in range(nh):
            out_q = Queue()
            mk_gal(h.data[i], out_q, s.info, i, final_gal[i], **keywords)
    
    print("----------Done---------")
    
    dictout=[]
    try:
        if not os.path.isdir(dir_out):
            os.mkdir(dir_out)
        f = open(dir_out + 'galaxies' + snout + '.txt', 'w')
    except:
        print("No filename is given.\n ")
    
    f.write(" nout    ID        x          y       z[Mpc]       vx      vy     vz[km/s]")
    f.write("    Reff[kpc]     Mstar    Mgas[Msun]  Rhalo[kpc]  boxtokpc  final_ID \n")    
    
    for i in range(nh):
        try:
            tmp =  out_q.get(timeout=0.1)
#            print(tmp)
            if tmp['id'] == 0:
                continue
            dd = tmp
            f.write("{:<4}   {:<4}   {:.5f}  {:.5f}  {:.5f}".format(nout,
                    dd['id'],dd['xc'],dd['yc'],dd['zc']))
            f.write("  {:.3f}  {:.3f}  {:.3f}".format(dd['vx'],dd['vy'],dd['vz']))
            f.write("  {:.6f}  {:.0f}  {:.0f}".format(dd['rgal'],dd['mstar'], dd['mgas']))
            f.write("  {:.5f}  {:.5f}  {:<4}   \n".format(dd['rhalo'],dd['boxtokpc'],dd['final_gal']))
            dictout.append(dd)
        except:
            continue
    
    f.close()    
    print("Text file written")

    catalog = pd.DataFrame(dictout).to_records()    
    # pickling also works!
    # The problem was memory being locked!

    with open(fcat, 'wb') as f:
        pickle.dump(catalog, f)
        
    star_all = 0
    dm_all = 0
    cell_all = 0
    s = 0
    # minimum stellar mass check only for the final snapshot galaxies,
    #N o more mstar_min test.
    print("------------------")
    print("main loop took ", time.time() - t0, "seconds")
