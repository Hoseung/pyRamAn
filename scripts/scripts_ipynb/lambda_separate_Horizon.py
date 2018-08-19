"""
    MODIFICATIONS
    
    2015.12.03  
        final_gal is removed as non-main-progenitors are also analized,
        there is no 1:1 correlation between final galaxies and earlier galaxies.
        
    2015.12.21
        Galaxies are saved in HDF format at selected nouts.
        (for example, 
         nouts_dump_gal = [187, 154, 130, 112,  98,  87,  67,  54,  37,  20]
         nouts close to zred = [0,0.2,0.4,0.6,0.8,1.0, 1.5, 2.0, 3.0, 5.0].)
         
    2016.01.01
        idx is stored in catalog along with id.
        
    2016.01.02
        lambdar measured after re-oriented.

    2016.03.10
        galaxy dump and catalog in the same directory.
        optional suffix to the catalog file.

    2016.03.26
        Only close galaxies were accounted before (i_gal_near), - why...??
        and that restriction is now gone. 

    2016.06.01
        No more need to read raw file every time.
        GalaxyMaker dump file (Galaxy, DM), and my own dump file (cell) are used.
        Currently script runs fine without cells file. 
        
    2016.06.04
        Not all of star, DM, and cell dump files are needed. 
        catalog now stored in dict, before transformed into pd.Dataframe().to_records()
        All modules are now moved to another file.
 
    2016.07.18
        worker(rscale_extract_cell=1.0) added. 
        At least all the gas inside the DM halo of a galaxy should be considered.
        gas inside the size of "stellar" galaxy has little correlation to the
        total amount of gas physically related to a galaxy. 
        Gas is colisional, and thus pressure can hold the gas even more extended than
        the virial radius.        

"""
import matplotlib
matplotlib.use("Agg")

import numpy as np
import collections
from galaxymodule import galaxy
import pickle
import load
from load.hydro import Hydro

from analysis.cal_lambda import *

def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx


def nouts_from_zreds(zreds, base='./'):
#if True:
    """
    Look for info in wdir/snapshots/outpu*/
    """
    import glob
    from load.info import Info
    finfos = glob.glob(wdir + "snapshots/output*/info*")
    nouts_info = []
    zreds_info = []
    for fn in finfos:
        ii = Info(fn=fn)
        nouts_info.append(ii.nout)
        zreds_info.append(ii.zred) 
    
    nouts_info = np.array(nouts_info)
    zreds_info = np.array(zreds_info)
    
    isort = np.argsort(zreds_info)
    nouts_info = nouts_info[isort]
    zreds_info = zreds_info[isort]

    nouts=[]
    for zz in zreds:
        nouts.append(nouts_info[find_closest(zreds_info, zz)])

    return nouts


def worker(gcatdata, additional_info,  hals, out_q, info, inds,
           reorient=False,
           galaxy_plot=False,
           galaxy_plot_dir='galaxy_plot/',
           region_plot=False,
           wdir = './',
           with_DM=False,
           with_cell=False,
           mk_gal_params={},
           cal_lambda_params={},
           rscale_extract_cell=1.0,
           **kwargs):
    
    import galaxy.make_gal as mkg
    from galaxymodule import galaxy
    nout = info.nout
    for i in inds:
        gcat_single = gcatdata[i]
        galid = gcat_single["id"]
        galidx = additional_info["idx"][i]
        tree_root_id = additional_info["tree_root_id"][i]

        gal = galaxy.Galaxy(halo=gcat_single, info=info)
        gm = load.rd_GM.rd_gal(nout, galid, base=wdir)

        # with_cell, with_DM are not implemented yet!
        if with_cell:
            gm.cell = extract_cell(hh, info, wdir, rscale=rscale_extract_cell)
        else:
            gm.cell = None
        if with_DM:
            if halid > 0 :
                gm.dm = load.rd_GM.rd_dm(info.nout, halid, base=wdir)
            else:
                gm.dm = extract_dm(hh, info, wdir)
                # The	se are not in the final units that I want,
                # but in the same units as GM outputs to make later conversion steps simpler.
                gm.dm['x'] = (gm.dm['x'] - 0.5) * info.pboxsize
                gm.dm['y'] = (gm.dm['y'] - 0.5) * info.pboxsize
                gm.dm['z'] = (gm.dm['z'] - 0.5) * info.pboxsize
                gm.dm['vx'] *=info.kms
                gm.dm['vy'] *=info.kms
                gm.dm['vz'] *=info.kms
                gm.dm['m'] * info.msun/1e11
        else:
            gm.dm = None
        good_gal = gal.mk_gal(gm.star, gm.dm, gm.cell,
                          **mk_gal_params)
        #do other things
        if not good_gal:
            print(galid, " Not a good galaxy")
            #out_q.put(gal.meta)
        else:
            print("galaxy {} is made \n".format(galid))

            try:
                lambdas = gal.cal_lambda_r_eps(**cal_lambda_params)
 
                gal.meta.lambda_arr, gal.meta.lambda_arrh, gal.meta.lambda_12kpc= lambdas[0]
                gal.meta.lambda_r,   gal.meta.lambda_rh,   gal.meta.lambda_r12kpc = lambdas[1]
                if galaxy_plot:
                     gal.plot_gal(fn_save = galaxy_plot_dir + str(nout).zfill(3) \
                                      + "_" + str(galid) + ".png", ioff=True)
 
                 #gal.meta.__dict__['idx'] = galidx
                gal.meta.idx = galidx
                gal.meta.tree_root_id = tree_root_id
                #gal.meta.__dict__['rhalo'] = hals.data['rvir'][i]
 
                print("   lambda calculation done. ID, IDx", galid, gal.meta.idx, gal.meta.tree_root_id)
 
                if reorient:
                    reori_pops = ["star"]
                    if with_cell:
                        reori_pops.append("cell")
                    if with_DM:
                        reori_pops.append("dm")
                        
                    gal.cal_norm_vec(["star"], dest=[0.,1.,0.])
                    gal.cal_rotation_matrix(dest = [0., 1., 0.])
                    gal.reorient(dest=[0., 1., 0.], pop_nvec = ['star'],
                                 pops=reori_pops, verbose=True)
                    lambdas = gal.cal_lambda_r_eps(**cal_lambda_params)
                                    #npix_per_reff=npix_lambda,
                                    #rscale=rscale_lambda, method='ellip')
## Be Because lambda is measured for ALL stars,
## B/ B/T must also be measured for ALL stars, not only for bound stars.
                    gal.meta.lambda_arr2, gal.meta.lambda_arr2h, gal.meta.lambda_arr2q = lambdas[0]
                    gal.meta.lambda_r2,   gal.meta.lambda_r2h  , gal.meta.lambda_r2q   = lambdas[1]
                    if galaxy_plot:
                        gal.plot_gal(fn_save = galaxy_plot_dir + str(nout).zfill(3) \
                                         + "_" + str(galid) + "_reori.png", ioff=True)
                out_q.put(gal.meta.__dict__)
            except:
                return

def main(base='./',
         ncore=1,
         nout_ini=20,
         nout_end=788,
         read_halo_list=True,
         out_dir="",
         cat_suffix="",
         is_gal = True,
         n_pseudo=1,
         w_reorient=False,
         w_galaxy_plot=True,
         w_galaxy_plot_dir= 'galaxy_plot/',
         w_region_plot=False,
         w_base='./',
         w_with_DM=False,
         w_with_cell=False,
         verbose=False,
         log_run_detail = False,
         sample_set=0,
         out_base='/scratch10/sukyi/'
         ):

    #from tree import treemodule
    import os
    import pickle
    #import tree.ctutils as ctu
    #import pandas as pd
    

    # worker options
    worker_params = {"reorient":       w_reorient,
                     "galaxy_plot":    w_galaxy_plot,
                     "galaxy_plot_dir":out_dir + w_galaxy_plot_dir,
                     "region_plot":    w_region_plot,
                     "wdir":           w_wdir,
                     "with_DM":        w_with_DM,
                     "with_cell":      w_with_cell
                    }
    multi = False                                             
    if ncore > 1: multi=True
    if multi: 
        import multiprocessing as mp
    else:
        import queue
 
    dir_cat = out_dir

    masscut_a = 1256366362.16
    masscut_b = -20583566.5218

    print("Check1")
    if out_dir == "":
        out_dir = str(sample_set) + '/'
        cat_suffix = str(sample_set)
        dir_cat = out_dir
        #out_dir = "out_default/"

    if out_dir[-1] != '/':
        out_dir = out_dir + '/'
        # need a tailing slash to be joined with sub-directories

    if not os.path.isdir(out_base + out_dir):
        os.mkdir(out_base + out_dir)

    out_base = out_base + out_dir

    # 27 : z=4;  37 : z=3;  20 : ~ z=5

    nout_fi = 788
    mstar_min = 1 # No mstar_min.
    npix=800
    lmax = 17
    ptypes=["star id pos mass vel time metal", "dm id pos mass vel"]

    if worker_params["galaxy_plot"]:
        if not os.path.isdir(worker_params["galaxy_plot_dir"]):
            os.mkdir(worker_params["galaxy_plot_dir"])


    voronoi_dict = dict(targetSN=15, plot=False, quiet=True)
    # voronoi tesselation need many points,
    # always use voronoi with pseudo particle option
    cal_lambda_params = dict(npix_per_reff=5,
                             rscale=3.0, 
                             method='ellip',
                             n_pseudo=n_pseudo,
                             verbose=False,
                             voronoi=None,#voronoi_dict,
                             mge_interpol = True)

    ## halo part -----------------------------------------------------------
    out_dir_cat = out_base + dir_cat + '/'

    # Load complete tree -----------------------------------------------------
    #if is_gal:
        # Galaxy tree
        #tree_path = 'GalaxyMaker/Trees/'
        #m_halo_min = 5e9 # minimum galaxy mass above which galaxies are searched for. 


##----------------------------------------------------------------------------------
    prgs = pickle.load(open(wdir + 'prgs_HM' + str(sample_set) + '.pickle', 'rb'))
    #tt = pickle.load(open(wdir + "Tree.pickle", "rb"))
    print("Check2")
    dd = np.genfromtxt(fname = "nstep_zred.txt", dtype=[("nstep", np.int),
                                                        ("zred", np.float32)],
                                                        usecols=[0,1])
    nsteps = dd["nstep"]
    zreds = dd["zred"]
    #nsteps = np.unique(tt["nstep"])[::-1]
    #zreds = np.unique(tt["zred"])
    print("Check2.5, loading halo")
    nouts = nouts_from_zreds(zreds, base=wdir)
    print(nouts)
    print("Check2.7, loading halo")
    for nout, nstep in zip(nouts, nsteps):
        if nout > nout_end:
            continue

        if nout < nout_ini:
            break
        print("Check3, loading halo")
        gcat = hmo.Halo(nout=nout, base=wdir, is_gal=is_gal)
        print("Check4, loading halo done")
        info = load.info.Info(nout=nout, base=wdir)
        ids_now = prgs["id"][prgs["nstep"]==nstep]
        idxs_now = prgs["idx"][prgs["nstep"]==nstep]
        trid_now = prgs["tree_root_id"][prgs["nstep"]==nstep]

        inds = ids_now[ids_now > 0] -1
        gcat_ok = gcat.data[inds]
        additional_info=np.zeros(len(gcat_ok), dtype=[("idx", "<i8"), ("tree_root_id", "<i8")])
        additional_info["idx"]= idxs_now[ids_now > 0]
        additional_info["tree_root_id"]= trid_now[ids_now > 0]

        print(nout, nout_fi)

        snout = str(nout)
        galaxy_plot_dir = out_base + 'galaxy_plot' + snout + '/'
        if not os.path.isdir(galaxy_plot_dir):
            os.mkdir(galaxy_plot_dir)

        mk_gal_params = dict(verbose=verbose,
                             mstar_min=mstar_min,
                             unit_conversion="GM")

        print("Worker")
    #   Multiprocessing --------------------------------------------------
        if multi:
            m = mp.Manager()
            out_q = m.Queue()
    
            # Distribute galaxies among cores
            new_inds=[]
            [new_inds.append([]) for i in range(ncore)]            
            for i in range(len(gcat_ok)):
                j = i % ncore
                new_inds[j].append(i)
            print(new_inds)
            processes = [mp.Process(target=worker, args=(gcat_ok, additional_info, None, out_q,
                        info, new_inds[i],
                        w_reorient,
                        w_galaxy_plot,
                        galaxy_plot_dir, 
                        w_region_plot,
                        w_wdir,
                        w_with_DM,
                        w_with_cell,
                        mk_gal_params,
                        cal_lambda_params)) for i in range(ncore)] 
            # Run processes
            for p in processes:
                p.start()
    
            # Exit the completed processes
            for p in processes:
                p.join()
        else:
            out_q = queue.Queue()
            new_inds = range(len(gcat_ok))
            worker(gcat_ok, additional_info, None, out_q, info, new_inds,
                   mk_gal_params=mk_gal_params,
                   cal_lambda_params=cal_lambda_params,
                   **worker_params)

        gcat = 0
        gcat_ok = 0
        print("----------Done---------")
    
        dictout=[]
        try:
            if not os.path.isdir(out_dir_cat):
                os.mkdir(out_dir_cat)
            f = open(out_dir_cat + 'galaxies' + snout + '.txt', 'w')
            #write header
            dd =  out_q.get(timeout=0.1)
            for key in sorted(dd.keys()):
                if isinstance(dd[key], collections.Iterable):
                    continue
                else:
                    f.write(str(key) + ", ")
                    f.write("\n")
            save_dict_scalar(dd, f, ", ")
            dictout.append(dd)
            for i in range(len(inds)):
                try:
                    dd =  out_q.get(timeout=0.1)
                    if dd['id'] == 0:
                        continue
                    save_dict_scalar(dd, f, ", ")
                    dictout.append(dd)
                except:
                    continue        
            f.close()    
            print("Text file written")
        except:
            print("No filename is given.\n ")
        
        fcat = out_dir_cat +"catalog" + snout + cat_suffix +".pickle"
        with open(fcat, 'wb') as f:
            pickle.dump(dictout, f)
            
        s = 0
        # minimum stellar mass check only for the final snapshot galaxies,
        # No more mstar_min test.
        print("------------------ \n")
        #print("main loop took ", time.time() - t0, "seconds")
    


if __name__ == "main":
    main()
