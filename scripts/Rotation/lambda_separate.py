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
 

"""
import matplotlib
matplotlib.use("Agg")

import numpy as np
import collections
from galaxy import galaxy
import pickle
import load

from analysis.cal_lambda import *

def worker(gals, hals, out_q, info, inds,
           dump_gal=False,
           reorient=False,
           galaxy_plot=False,
           galaxy_plot_dir='galaxy_plot/',
           region_plot=False,
           wdir='./',
           with_DM=True,
           with_cell=True,
           mk_gal_params={},
           cal_lambda_params={},
           **kwargs):



    import galaxy.make_gal as mkg
    from galaxy import galaxy

    nout = info.nout
    if type(inds) == int:
        inds = [inds]
    for i in inds:
        #print(mp.current_process())
        print("\n \nThis is {}-th galaxy".format(i))
        print("halo: {}".format(hals.data['id'][i]), )
        print("gal: {}".format(gals.data['id'][i]), )
        
        gg = gals.data[i]
        hh = hals.data[i]

        galid = gg['id']
        halid = hh['id']
        gm = load.rd_GM.rd_gal(info.nout, galid, wdir=wdir)
        # print('!!!! mima vx in gm', min(gm.star['vx']), max(gm.star['vx']))

        if with_cell:
            gm.cell = load.rd_GM.rd_cell(info.nout, galid, wdir=wdir)
        else:
            gm.cell = None
        if with_DM:
            if halid > 0 :
                gm.dm = load.rd_GM.rd_dm(info.nout, halid, wdir=wdir)
            else:
                gm.dm = extract_dm(hh, info.nout, wdir)
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

        gal = galaxy.Galaxy(halo=gg, info=info)
        good_gal = gal.mk_gal(gm.star, gm.dm, gm.cell,
                              unit_conversion="GM",
                              verbose=False)
        if with_cell:
            print("GAS MASS in solar mass", gal.meta.mgas)

        galid = gal.meta.id

        #do other things
        if not good_gal:
            print(galid, " Not a good galaxy")
            #out_q.put(gal.meta)
        else:
            print("galaxy {} is made \n".format(galid))
            lambdas = gal.cal_lambda_r_eps(**cal_lambda_params)
#                           npix_per_reff=npix_lambda,
#                           rscale=rscale_lambda, method='ellip', verbose=False)

            gal.meta.lambda_arr, gal.meta.lambda_arrh, gal.meta.lambda_12kpc= lambdas[0]
            gal.meta.lambda_r,   gal.meta.lambda_rh,   gal.meta.lambda_12kpc = lambdas[1]
            if galaxy_plot:
                gal.plot_gal(fn_save = galaxy_plot_dir + str(nout).zfill(3) \
                                 + "_" + str(galid) + ".png", ioff=True)

            gal.meta.__dict__['idx'] = hals.data['idx'][i]
            gal.meta.__dict__['rhalo'] = hals.data['rvir'][i]

            print("   lambda calculation done. ID, IDx", galid, gal.meta.idx)

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
# Because lambda is measured for ALL stars,
# B/T must also be measured for ALL stars, not only for bound stars.
                gal.meta.lambda_arr2, gal.meta.lambda_arr2h, gal.meta.lambda_arr2q = lambdas[0]
                gal.meta.lambda_r2,   gal.meta.lambda_r2h  , gal.meta.lambda_r2q   = lambdas[1]
            if galaxy_plot:
#                gal.cal_b2t(ptype='star',
#                            bound_only=False,
#                            hist=False,
#                            proj='y')
                gal.plot_gal(fn_save = galaxy_plot_dir + str(nout).zfill(3) \
                                     + "_" + str(galid) + "_reori.png", ioff=True)
        out_q.put(gal.meta.__dict__)
        #print("where are you gone, dict?", out_q.get())


def main(wdir='./',
         ncore=4,
         nout_ini=37,
         nout_end=187,
         read_halo_list=True,
         out_dir="",
         cat_suffix="",
         r_cluster_scale=2.9,
         with_DM=False,
         with_cell=False):

    from tree import treemodule
    import os
    import pickle
    import tree.ctutils as ctu
    import pandas as pd
    
    multi = False # 
    if multi: 
        import multiprocessing as mp
    else:
        import queue
 
    is_gal = True
    dump_gal = True
    nout_complete = 87 # at least complete up to z = 1
    
    masscut_a = 1256366362.16
    masscut_b = -20583566.5218

    # optional parameters ----------------------------------------------------
    lambda_method = 'ellip' 
    galaxy_plot = False
    reorient = True
    verbose=False
    region_plot = False

    verbose=False
    dir_cat = out_dir

    if out_dir == "":
        out_dir = "out_/"

    if out_dir[-1] != '/':
        out_dir = out_dir + '/'
        # need a tailing slash to be joined with sub-directories

    if not os.path.isdir(wdir + out_dir):
        os.mkdir(wdir + out_dir)

    out_base = wdir + out_dir

    # 27 : z=4;  37 : z=3;  20 : ~ z=5
    if nout_ini == "":
        nout_ini = 37
    else:
        nout_ini = int(nout_ini)

    if nout_end == "":
        nout_end = 187
    else:
        nout_end = int(nout_end)

    nout_fi = 187
    nouts = range(nout_fi, nout_ini -1, -1)
    #nouts = [187, 121, 87, 54, 37]
    #nouts_dump_gal = [187, 154, 130, 121, 112,  98,  87,  67,  54,  37,  20]
    nouts_dump_gal = []
    try: nouts_dump_gal
    except NameError: nouts_dump_gal = None
    #----------------------------------------------------------------------

    # Only galaxies above this stellar mass at the final snapshot are considered.
    mstar_min = 5e9
    mk_gal_rscale = 1.1 # unit of Rvir,galaxy
    #r_cluster_scale = 2.9 # maximum radius inside which galaxies are searched for
    npix=800


    # worker options
    worker_params = {"dump_gal":False,
                   "reorient":False,
                   "galaxy_plot":True,
                   "galaxy_plot_dir":out_dir + 'galaxy_plot/',
                   "region_plot":False,
                   "wdir":'./',
                   "with_DM":False,
                   "with_cell":False}

    if worker_params["galaxy_plot"]:
        if not os.path.isdir(worker_params["galaxy_plot_dir"]):
            os.mkdir(worker_params["galaxy_plot_dir"])


    voronoi_dict = dict(targetSN=15, plot=False, quiet=True)
    # voronoi tesselation need many points,
    # always use voronoi with pseudo particle option
    cal_lambda_params = dict(npix_per_reff=5,
                             rscale=3.0, 
                             method='ellip',
                             n_pseudo=60,
                             verbose=False,
                             voronoi=voronoi_dict,
                             mge_interpol = True)


    lmax = 19
    ptypes=["star id pos mass vel time metal", "dm id pos mass vel"]

    ## halo part -----------------------------------------------------------
    dir_out = wdir + dir_cat + '/'


    # Load complete tree -----------------------------------------------------
    if is_gal:
        # Galaxy tree
        tree_path = 'GalaxyMaker/Trees/'
        m_halo_min = 5e9 # minimum galaxy mass above which galaxies are searched for. 
    else:
        # halo tree
        tree_path = 'halo/Trees/'
        m_halo_min = 2e10 # minimum halo mass. 
    try:
        alltrees = pickle.load(open(wdir + tree_path + "extended_tree.pickle", "rb" ))
        print("Loaded an extended tree")
    except:
        alltrees = treemodule.CTree()
        alltrees.load(filename= wdir + tree_path + 'tree_0_0_0.dat')
        # Fix nout -----------------------------------------------------
        nout_max = alltrees.data['nout'].max()
        alltrees.data['nout'] += nout_fi - nout_max
        print("------ NOUT fixed")
        alltrees.data = ctu.augment_tree(alltrees.data, wdir, is_gal=is_gal)
        print("------ tree data extended")
        pickle.dump(alltrees, open(wdir + tree_path + "extended_tree.pickle", "wb" ))

    # Reading tree done
    td = alltrees.data


    if read_halo_list:
        try:
            print("loading pickled halo list done:")
            prg_only_tree = pickle.load(open(wdir + "prg_only_tree.pickle", 'rb'))
            read_halo_list = True
        except:
            read_halo_list = False

    if not read_halo_list:
        info = load.info.Info(nout=nout_fi, base=wdir, load=True)
        prg_only_tree = get_sample_tree(alltrees, info,
                        wdir=wdir,
                        nout_ini=nout_ini,
                        nout_fi =nout_fi,  
                        is_gal = is_gal,
                        r_cluster_scale = r_cluster_scale,
                        m_halo_min = m_halo_min,
                        nout_complete = nout_complete)
        pickle.dump(prg_only_tree, open(wdir + 'prg_only_tree.pickle', 'wb'))


    with open(wdir + 'lambda_mp_status.txt', 'w') as f:
        f.write("mstar_min = " + str(mstar_min) + "\n")
        f.write("Rscale cluster : " + str(r_cluster_scale))
        f.write("Rscale mk_gal : " + str(mk_gal_rscale))
        f.write("npix : " + str(npix))
        save_dict_scalar(cal_lambda_params, f, "\n")
        #f.write("Rscale lambda calculation : " + str(rscale_lambda))
        #f.write("npix per 1Reff for lambda calculation : " + str(npix_lambda))
        f.write("ptypes : \n")
        for i in ptypes:
            f.write("  " + str(i) + "\n")


    for nout in nouts:
        print(nout, nout_fi)

        snout = str(nout)
        if nout > nout_end:
            continue
        if nout in nouts_dump_gal:
            dump_gal = True
            print("DUMP TRUE___________________________________")
        else:
            dump_gal = False

        fcat = dir_out +"catalog" + snout + cat_suffix +".pickle"

        galaxy_plot_dir = out_base + 'galaxy_plot' + snout + '/'
        if not os.path.isdir(galaxy_plot_dir):
            os.mkdir(galaxy_plot_dir)

        info = load.info.Info(nout=nout, base=wdir, load=True)
        if nout < 187:
            mstar_min = 2 * get_mstar_min(info.aexp)
        print("Mstar now, min: {:.2e} {:.2e}".format(2 * \
                    (masscut_a * info.aexp + masscut_b), mstar_min))

        allgal, allhal = get_sample_gal(wdir, nout, info, prg_only_tree, mstar_min)

###
#  required quantities
#  halo.data[] x, y, z, vx, vy, vz, m, mvir, r, rivr, id, idx, 
#  halo.idlists
#
        nh = len(allgal.data)
        print("Total # galaxies to analyze",nh)
        print("# complete-tree galaxies",sum(allhal.data['idx'] > 0))
        print("# non-complete-tree galaxies",sum(allhal.data['idx'] < 0))


        mk_gal_params = dict(verbose=verbose,
                             mstar_min=mstar_min,
                             unit_conversion="GM")

    #   Multiprocessing --------------------------------------------------
        if multi:
            m = mp.Manager()
            out_q = m.Queue()
            print("Analyzing galaxies inside {} halos".format(nh))
         
            # Distribute galaxies among cores
            inds=[]
            [inds.append([]) for i in range(ncore)]
         
            for i in range(nh):
                j = i % ncore
                inds[j].append(i)
         
            processes = [mp.Process(target=worker, args=(allgal, allhal, out_q,
                        info, inds[i],
                        dump_gal,
                        reorient,
                        galaxy_plot,
                        galaxy_plot_dir,
                        region_plot, 
                        wdir,
                        with_DM,
                        with_cell,
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
            worker(allgal, allhal, out_q, info, range(nh),
                   mk_gal_params=mk_gal_params,
                   cal_lambda_params=cal_lambda_params,
                   **worker_params)
        
        print("----------Done---------")
    
        dictout=[]
        try:
            if not os.path.isdir(dir_out):
                os.mkdir(dir_out)
            f = open(dir_out + 'galaxies' + snout + '.txt', 'w')
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
            for i in range(nh-1):
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
        
        with open(fcat, 'wb') as f:
            pickle.dump(dictout, f)
            
        s = 0
        # minimum stellar mass check only for the final snapshot galaxies,
        # No more mstar_min test.
        print("------------------ \n")
        #print("main loop took ", time.time() - t0, "seconds")
    
