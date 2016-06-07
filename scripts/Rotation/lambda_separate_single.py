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

"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt 

import numpy as np
import utils.sampling as smp

import collections
from galaxy import galaxy
import utils.match as mtc
import tree.ctutils as ctu
import pickle
import load
import tree.halomodule as hmo 
import general


from analysis.cal_lambda import *


#%%
    
def worker(gals, hals, out_all, out_young, out_old, info,
           dump_gal=False,
           reorient=False,
           lambda_method='ellip',
           rscale_lambda=3,
           npix_lambda=10,
           galaxy_plot=False,
           galaxy_plot_dir='galaxy_plot/',
           region_plot=False,
           wdir='./',
           with_DM=True,
           with_cell=True,
           age_cut=10.0,
           **kwargs):

    import galaxy.make_gal as mkg
    from galaxy import galaxy
    import copy
    import utils.cosmology

    nout = info.nout
    #if type(inds) == int:
    #    inds = [inds]
    inds=[0]
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
            if halid > 0:
                gm.dm = load.rd_GM.rd_dm(info.nout, halid, wdir=wdir)
            else:
                gm.dm = extract_dm(hh, info.nout, wdir)
                # These are not in the final units that I want,
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

        # copy - galaxy will be reoriented.
        gm_young = copy.copy(gm)
        gm_old = copy.copy(gm)
 
        gm.star['time'] = utils.cosmology.time2gyr(gm.star['time'],
                                     z_now = info.zred,
                                     info=info)


        print(min(gm.star['time']), max(gm.star['time']))

        gm_young.star = gm.star[gm.star['time'] < age_cut]
        gm_old.star = gm.star[gm.star['time'] > age_cut]
 
        print("all/ young/ old stelar mass \n {:.2e}  {:.2e}  {:.2e}".format(
               sum(gm.star['m']),
               sum(gm_young.star['m']),
               sum(gm_old.star['m'])))

        gal = _make_n_run(gg, gm, info, hals, i=i,
           with_DM=with_DM, with_cell=with_cell,
           reorient=True,
           npix_lambda=npix_lambda,
           rscale_lambda=rscale_lambda,
           fn_save=galaxy_plot_dir + str(nout).zfill(3)+ "_all_" + str(galid),
           galaxy_plot=galaxy_plot)

        out_all.put(gal.meta.__dict__)

        try:
            print(" \n Now young stellar particles ({})".format(len(gm_young.star)))
            gal = _make_n_run(gg, gm_young, info, hals, i=i,
               with_DM=with_DM, with_cell=with_cell,
               reorient=False,
               npix_lambda=npix_lambda,
               rscale_lambda=rscale_lambda,
               fn_save=galaxy_plot_dir + str(nout).zfill(3)+ "_young_" + str(galid),
               galaxy_plot=galaxy_plot)

            out_young.put(gal.meta.__dict__)
        except:
            pass


#        if len(gm_old.star) > 1e3:
        try:
            print(" \n And old stellar particles ({})".format(len(gm_old.star)))
            gal = _make_n_run(gg, gm_old, info, hals, i=i,
               with_DM=with_DM, with_cell=with_cell,
               reorient=False,
               npix_lambda=npix_lambda,
               rscale_lambda=rscale_lambda,
               fn_save=galaxy_plot_dir + str(nout).zfill(3)+ "_old_" + str(galid),
               galaxy_plot=galaxy_plot)
         
            out_old.put(gal.meta.__dict__)
        except:
            pass

        gal = 0
        gm = 0
        gm_young = 0
        gm_old = 0


def _make_n_run(gg, gm, info, hals,
                i=0,
                with_DM=True,
                with_cell=True,
                reorient=True,
                npix_lambda=5,
                rscale_lambda=3.0,
                fn_save='pic.png',
                galaxy_plot=True):

    gal = galaxy.Galaxy(halo = gg,
                        radius_method='eff',
                        info=info)
    #print('!!!! mima vx gm.dm', min(gm.dm['vx']), max(gm.dm['vx']))
    good_gal = gal.mk_gal(gm.star, gm.dm, gm.cell,
                          unit_conversion="GM",
                          verbose=False)

    #print("GAS MASS in solar mass", gal.meta.mgas)

    galid = gal.meta.id
    gal.meta.nout = info.nout

    #do other things
    if not good_gal:
        print(galid, " Not a good galaxy")
        #out_q.put(gal.meta)
    else:
        print("galaxy {} is made \n".format(galid))
        lambdas = gal.cal_lambda_r_eps(npix_per_reff=npix_lambda,
                       rscale=rscale_lambda, method='ellip', verbose=False) # calculate within 1.0 * reff
        gal.meta.lambda_arr, gal.meta.lambda_arrh, gal.meta.lambda_12kpc= lambdas[0]
        gal.meta.lambda_r,   gal.meta.lambda_rh,   gal.meta.lambda_12kpc = lambdas[1]
        if galaxy_plot:
            gal.plot_gal(fn_save = fn_save, ioff=True)

        gal.meta.__dict__['idx'] = hals.data['idx'][i]
        gal.meta.__dict__['rhalo'] = hals.data['rvir'][i]

        print("   lambda calculation done. ID, IDx", galid, gal.meta.idx)

        gal.cal_norm_vec(["star"], dest=[0.,1.,0.])
        if reorient:
            reori_pops = ["star"]
            if with_cell is None:
                reori_pops.append("cell")
            if with_DM is None:
                reori_pops.append("dm")

            gal.cal_rotation_matrix(dest = [0., 1., 0.])
            gal.reorient(dest=[0., 1., 0.], pop_nvec = ['star'],
                         pops=reori_pops, verbose=True)
            lambdas = gal.cal_lambda_r_eps(npix_per_reff=npix_lambda,
                            rscale=rscale_lambda, method='ellip')
            # Because lambda is measured for ALL stars,
            # B/T must also be measured for ALL stars, not only for bound stars.
            gal.meta.lambda_arr2, gal.meta.lambda_arr2h, gal.meta.lambda_arr2q = lambdas[0]
            gal.meta.lambda_r2,   gal.meta.lambda_r2h  , gal.meta.lambda_r2q   = lambdas[1]
#        if galaxy_plot:
#            gal.cal_b2t(ptype='star',
#                        bound_only=False,
#                        hist=False,
#                        proj='y')
#            gal.plot_gal(fn_save = fn_save + "_reori.png", ioff=True)

    return gal



def dump_queue(qq, fname, safe=True):
    dictout=[]
    f = open(fname + '.txt', 'w')
    #write header
    dd =  qq.get(timeout=0.1)
    for key in sorted(dd.keys()):
        if isinstance(dd[key], collections.Iterable):
            continue
        else:
            f.write(str(key) + ", ")
            f.write("\n")
    save_dict_scalar(dd, f, ", ")
    dictout.append(dd)
    while not qq.empty():
        dd =  qq.get(timeout=0.1)
        if dd['id'] == 0:
            continue
        save_dict_scalar(dd, f, ", ")
        dictout.append(dd)
    f.close()  

    if safe:
        catalog = dictout
    else:
        import pandas as pd
        catalog = pd.DataFrame(dictout).to_records() 

    with open(fname + '.pickle', 'wb') as f:
        pickle.dump(catalog, f)



def main(galidx,
         wdir='./',
         ncore=1,
         nout_ini=37,
         nout_end=187,
         read_halo_list=True,
         cat_suffix="",
         r_cluster_scale=2.9,
         with_cell = False,
         with_DM = False):

    from tree import treemodule
    import os
    import pickle
    import tree.ctutils as ctu
    from queue import Queue
    #import multiprocessing as mp
 
    
    #multi = 1 # 
    is_gal = True
    dump_gal = True
    nout_complete = 87 # at least complete up to z = 1

    masscut_a = 1256366362.16
    masscut_b = -20583566.5218
    
    age_cut = 9
    mstar_min = 1e20 # no more galaxies than I chose.

    # optional parameters ----------------------------------------------------
    lambda_method = 'ellip' 
    galaxy_plot = True
    reorient = True
    verbose=False
    region_plot = False
    dump_gal = False # No need to save galaxy anymore.
                     # Now star, dm, cell are stored in separate files. 

    verbose=False
    out_dir = "gal_" + str(galidx) + '/'
    dir_cat = out_dir

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
    mstar_min = 5e9
    # Only galaxies above this stellar mass at the final snapshot are considered.
    mk_gal_rscale = 1.1 # unit of Rvir,galaxy
    npix=800
    rscale_lambda = 3.0 # Reff unit.
    npix_lambda = 5 # per 1Reff

    ## halo part -----------------------------------------------------------
    dir_out = wdir + out_dir + '/'

    # Load complete tree -----------------------------------------------------
    # Galaxy tree
    tree_path = 'GalaxyMaker/Trees/'
    m_halo_min = 5e9 # minimum galaxy mass above which galaxies are searched for. 


    if read_halo_list:
        try:
            print("loading pickled halo list done:")
            prg_only_tree = pickle.load(open(wdir + "prg_only_tree.pickle", 'rb'))
            read_halo_list = True
        except:
            read_halo_list = False

    if not read_halo_list:
        # Load tree 
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
 
        # build progenitor-only tree
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


    # Tree of only one galaxy.
    prg_only_tree = ctu.extract_main_tree(prg_only_tree, idx=galidx)
  
    out_all = Queue()
    out_young = Queue()
    out_old = Queue()    


    for nout in nouts:
        print(nout, nout_fi)

        snout = str(nout)
        if nout > nout_end:
            continue

        galaxy_plot_dir = out_base
        info = load.info.Info(nout=nout, base=wdir, load=True)
        
        # What to do with Phantom galaxy??
        this_gal = prg_only_tree[prg_only_tree['nout'] == nout]
        if this_gal['phantom'] == 1:
            continue

        allgal, allhal = get_sample_gal(wdir,
                             nout, info, prg_only_tree, mstar_min)

 
        keywords = dict(rscale = mk_gal_rscale,
                    verbose=verbose,
                    mstar_min=mstar_min)#, dump_gal = dump_gal,

        worker(allgal,allhal, out_all, out_young, out_old,
               info,
               dump_gal=dump_gal,
               reorient=reorient,
               lambda_method=lambda_method,
               rscale_lambda=rscale_lambda,
               npix_lambda=npix_lambda,
               galaxy_plot=galaxy_plot,
               galaxy_plot_dir=galaxy_plot_dir,
               region_plot=region_plot,
               wdir=wdir,
               with_DM=with_DM,
               with_cell=with_cell,
               age_cut=age_cut,
               **keywords)
        
    print("-------- All nouts Done---------")

    
    if not os.path.isdir(dir_out):
        os.mkdir(dir_out)
    
    fname = dir_out + 'galaxy_' + str(galidx)

    dump_queue(out_all, dir_out + 'galaxy_' + str(galidx))
    dump_queue(out_young, dir_out + 'galaxy_young' + str(galidx))
    dump_queue(out_old, dir_out + 'galaxy_old' + str(galidx))
  
    m=0
