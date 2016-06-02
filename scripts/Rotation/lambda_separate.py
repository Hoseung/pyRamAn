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

#%%
def extract_halos_within(halos, i_center, info, dist_in_mpc=1.0):
    xc = halos['x'][i_center]
    yc = halos['y'][i_center]
    zc = halos['z'][i_center]

    xx = halos['x']
    yy = halos['y']
    zz = halos['z']

    dd = np.multiply(distance_to([xc,yc,zc], [xx,yy,zz]), info.pboxsize)

    return (dd < (dist_in_mpc))


def distance_to(xc, xx):
    return np.sqrt([(xc[0] - xx[0])**2 + (xc[1] - xx[1])**2 + (xc[2] - xx[2])**2])[0]


def all_gals_org(treedata, final_gals, nout_ini=None, nout_fi=None):
    """
       build a list of all progenitors of the final_gals from nout_fi up to nout_ini
       [ [final_gals (at nout = nout_fi) ], 
         [all progenitors of final gals at nout = nout_fi -1],
         [ '' at the at nout = nout_fi -2 ], ...]
    """
    if nout_ini == None:
        nout_ini = min(treedata['nout'])
    if nout_fi == None:
        nout_fi = max(treedata['nout'])

    all_gals_at_nouts = []
    for inout, nout in enumerate(range(nout_ini, nout_fi+1)):
        all_gals_this_nout = []
        tree_now = treedata[np.where(treedata['nout'] == nout)]

        for finalgal in final_gals:
            i_gals_include = np.where(tree_now['tree_root_id'] == finalgal)[0]
            [all_gals_this_nout.append(gal) for gal in tree_now['id'][i_gals_include]]

        all_gals_at_nouts.append(all_gals_this_nout)

    return all_gals_at_nouts


def all_gals(treedata, final_gals, nout_ini=None, nout_fi=None):
    """
       build a list of all progenitors of the final_gals from nout_fi up to nout_ini
       [ [final_gals (at nout = nout_fi) ], 
         [all progenitors of final gals at nout = nout_fi -1],
         [ '' at the at nout = nout_fi -2 ], ...]
    """
    if nout_ini == None:
        nout_ini = min(treedata['nout'])
    if nout_fi == None:
        nout_fi = max(treedata['nout'])

    inds=[]
    for finalgal in final_gals:
        inds.extend(np.where(treedata['tree_root_id'] == finalgal)[0])
        # Add main progenitor tag
    
    return treedata[inds]


def halo_from_tree(tree_element, info):
    import tree.halomodule as hmo

    dtype_halo = [('id', '<i4'), ('idx', '<i4'), ('m', '<f4'), ('mvir', '<f4'),
              ('r', '<f4'), ('rvir', '<f4'),
              ('x', '<f4'), ('y', '<f4'), ('z', '<f4'),
              ('vx', '<f4'), ('vy', '<f4'), ('vz', '<f4')]

    cboxsize = 200.

    nout = tree_element['nout']
    h = hmo.Halo(nout=nout, halofinder='HM', info=info, is_gal=True)

    h.data = np.recarray(len(tree_element), dtype=dtype_halo)
    h.nout = nout

    h.data['m'] = tree_element['m']
    h.data['mvir'] = tree_element['mvir']
    h.data['x'] = tree_element['x'] / cboxsize
    h.data['y'] = tree_element['y'] / cboxsize
    h.data['z'] = tree_element['z'] / cboxsize # Mpc/h -> code unit
    h.data['vx'] = tree_element['vx']
    h.data['vy'] = tree_element['vy']
    h.data['vz'] = tree_element['vz']
    h.data['r'] = tree_element['r'] # already in code unit
    h.data['rvir'] = tree_element['rvir'] / (cboxsize * 1000) # kpc/h -> code unit
    h.data['id'] = tree_element['Orig_halo_id']
    h.data['idx'] = tree_element['id']
    h.aexp = tree_element['aexp']

    return h

def dist(data, center):
    return np.sqrt(np.square(center['x'] - data['x']) +
                np.square(center['y'] - data['y']) +
                np.square(center['z'] - data['z']))

def dist2(data, center):
    return (np.square(center['x'] - data['x']) +
         np.square(center['y'] - data['y']) +
         np.square(center['z'] - data['z']))


def find_halo_near(data, center, rscale=1.0):
    import numpy as np
    i_gal_ok = dist2(data, center) <  np.square(rscale * center['rvir'])
    return i_gal_ok


def unique(a,b):
    a = np.concatenate((a,b))
    a = np.sort(a)
    b = np.diff(a)
    b = np.r_[1, b]
    return a[b != 0]


def get_mstar_min(aexp):
    masscut_a = 1256366362.16
    masscut_b = -20583566.5218
    
    return masscut_a * aexp + masscut_b



def associate_gal_hal(allgal, allhal, plot_check=False, dir_out=""):
    """ 
    associate halos with galaxies. 
    Arbitrary maching parameters are used.
    """
    import numpy as np
    
    def dist(data, center):
        return np.sqrt(np.square(center['x'] - data['x']) +
                np.square(center['y'] - data['y']) +
                np.square(center['z'] - data['z']))
    
    def distv(halo, center):
        norm = np.sqrt(np.square(center['vx'] - halo.vx) + 
                       np.square(center['vy'] - halo.vy) + 
                       np.square(center['vz'] - halo.vz))
    
        return norm

    i0=[] # matched list
    i1=[] # unmatched list

    newhals = np.recarray(len(allgal.data), dtype=allhal.data.dtype) # equivalent with allhal.data
    
    for i, gal in enumerate(allgal.data):
        dd = dist(allhal.data, gal)# 3d distance. simply sqrt((x-xc)**2 + ...)
        d_sort = np.argsort(dd)
        # if closest halo is closer by 0.1 than the second and close than 10kpc/h, good match.
        if (dd[d_sort[0]] < 0.1 * dd[d_sort[1]]) and (dd[d_sort[0]] < 5e-5):
            gal['hosthalo'] = allhal.data['id'][d_sort[0]]
            i0.append(i)
            newhals[i] = allhal.data[d_sort[0]]
    
        # if closest halo is closer by 0.3 and than the second and closer than 5kpc/h, good match.
        elif (dd[d_sort[0]] < 0.3 * dd[d_sort[1]]) and (dd[d_sort[0]] < 2.5e-5):
            gal['hosthalo'] = allhal.data['id'][d_sort[0]]
            i0.append(i)
            newhals[i] = allhal.data[d_sort[0]]
        else:
            # if closet than 40kpc/h and has similar velocity.
            dv = distv(allhal.data, gal)
            d_nomi = dd < 2e-4 # within 40kpc!! 
            v_nomi = dv < 150 
            dvnomi = d_nomi * v_nomi
            if sum(dvnomi) > 1:
                dddv = dd[dvnomi]/2e-4 + dv[dvnomi]/150
                dddv_sort = np.argsort(dddv)
                gal['hosthalo'] = allhal.data['id'][dddv_sort[0]]
                i0.append(i)
                newhals[i] = allhal.data[dddv_sort[0]]
            else:
                # otherwise non-matched.
                gal['hosthalo'] = -1
                i1.append(i)
                # copy galaxy catalog data to new 'fake' halo catalog. 
                newhals[i]['x'] = gal['x']
                newhals[i]['y'] = gal['y']
                newhals[i]['z'] = gal['z']
                newhals[i]['vx'] = gal['vx']
                newhals[i]['vy'] = gal['vy']
                newhals[i]['vz'] = gal['vz']
                newhals[i]['r'] = gal['r'] * 2 # arbitrary
                newhals[i]['m'] = gal['m'] * 2 # arbitrary
                newhals[i]['rvir'] = gal['rvir'] * 2 # arbitrary
                newhals[i]['mvir'] = gal['mvir'] * 2 # arbitrary
                # small radius assumes most of the DM, gas have been lost.
                # star is the only major component.
                newhals[i]['id'] = -1 * gal['id'] # negative ID = fake halo.

    allhal.data = newhals
    
    if plot_check:
        from draw import pp
        import matplotlib.pyplot as plt
        #plt.ioff()
        fig, ax = plt.subplots(1)
        pp.pp_halo(allhal, 400, axes= ax, verbose=False, colors=10, rscale=0.3)
        pp.pp_halo(allgal, 400, axes= ax, verbose=False, colors=200, rscale=0.3)
        plt.savefig(dir_out + "associate_gal_hal.png")
        plt.close()
    
    
    return allhal

    
    #%%
    #import numpy as np

def get_sample_tree(alltrees,
                    info,
                    wdir = './',
                    nout_ini=57,
                    nout_fi = 18,  
                    is_gal = True,
                    r_cluster_scale = 2.9,
                    m_halo_min = 5e9,
                    nout_complete = 87):
                        
    import load
#    import pickle
    import tree.ctutils as ctu
    
#    gd = general.defaults.Default()
    
#    info = load.info.Info(nout=nout_fi, base=wdir, load=True)
    # load galaxy tree 
    
#    alltrees = pickle.load(open(wdir + \
#                    gd.dir_galaxy_tree + "extended_tree.pickle", "rb"))
    
    td = alltrees.data
    
    # halo catalog
    hhal = hmo.Halo(base=wdir, nout=nout_fi, halofinder='HM', load=True, is_gal=False)
    # cluster radius 
    
    i_center = np.where(hhal.data['np'] == max(hhal.data['np']))[0]
    r_cluster = hhal.data['rvir'][i_center].squeeze() * hhal.info.pboxsize
    
    # galaxies
    hh = hmo.Halo(base=wdir, nout=nout_fi, halofinder='HM', info=info, load=True, is_gal=is_gal)
    i_center = np.where(hh.data['np'] == max(hh.data['np']))[0]
    # galaxies that are within r_cluster_scale * the cluster halo from the BCG 
    # (not the cluster halo center, although two must be virtually identical)
    
    # All galaxies inside the cluster radius * r_scale
    i_satellites = extract_halos_within(hh.data, i_center, info, dist_in_mpc = r_cluster * r_cluster_scale)
    print("Total {} galaxies \n {} galaxies are within {} times the cluster virial radius, {} Mpc".format(
          len(i_satellites),sum(i_satellites), r_cluster_scale, r_cluster_scale * r_cluster))
    
    # Above a mass cut at nout_fi
    # halos found inside the cluster and have complete tree back to nout_ini
    large_enough = hh.data['mvir'] > m_halo_min
    halo_list = hh.data['id'][i_satellites * large_enough]
    final_ids = ctu.check_tree_complete(td, nout_complete, nout_fi, halo_list, idx=False) # 87: z = 1
    
    
    # build list of progenitors (all mass)
    tt_final = td[td['nout'] == nout_fi]
    final_gals_idx = [tt_final['id'][tt_final['Orig_halo_id'] == final_gal] for final_gal in final_ids]
    ngals = len(final_gals_idx)
    print(" {} galaxies have complete tree up to nout = {}".format(ngals, nout_complete))
    # Search for all galaxies that are listed in the trees of final_gals
    prg_only_tree = all_gals(td, final_gals_idx, nout_fi=nout_fi, nout_ini=nout_ini)
    #all_gals_in_trees = prg_only_tree['']

    # idx list of galaxies at the current nout (from all_gals_in_trees list)
    #t_now = prg_only_tree[prg_only_tree['nout'] == nout]
    #idxs_tree_now = all_gals_in_trees[inout]

    return prg_only_tree

#%%
def get_sample_gal(wdir, nout, info, prg_only_tree, mstar_min): 
    import utils.match as mtc
    #gals_in_tree_now = halo_from_tree(t_now[mtc.match_list_ind(t_now['id'], np.array(idxs_tree_now))], info)
    gals_in_tree_now = prg_only_tree[prg_only_tree['nout'] == nout]
    id_now = gals_in_tree_now['Orig_halo_id'] # this is Orig_halo_id
    
    # Galaxies near the cluster
    allhal = hmo.Halo(base=wdir, nout=nout, is_gal=False, halofinder='HM', return_id=False, load=True)
    cluster_now = allhal.data[allhal.data.np.argmax()]
    
    max_dist_prg = max(np.sqrt(np.square(cluster_now['x'] * 200 - gals_in_tree_now['x']) + 
                               np.square(cluster_now['y'] * 200 - gals_in_tree_now['y']) + 
                               np.square(cluster_now['z'] * 200 - gals_in_tree_now['z']))) # in Mpc/h
                        
    print(nout, max_dist_prg)
    
    # unless you have "perfect" trees, do not rely on what tree gives you.
    #mstar_min = min(gals_in_tree_now['m'][(gals_in_tree_now['mmp'] == 1) * (gals_in_tree_now['phantom'] == 0)])
    
    # Galaxy with enough stellar mass
    # -> there is no stellar mass cut when building all_gals_in_trees list.
    allgal = hmo.Halo(base=wdir, nout=nout, is_gal=True, halofinder='HM',
                      return_id=False, load=True)
    
    dd_cat = np.sqrt(np.square(cluster_now['x'] - allgal.data['x']) + 
                     np.square(cluster_now['y'] - allgal.data['y']) + 
                     np.square(cluster_now['z'] - allgal.data['z'])) * 200
    
    igal_ok_cat = (dd_cat < max_dist_prg) * (allgal.data['m'] > mstar_min)
    print("(catalogue) # galaxies more massive than {:.2e} at nout ="
    " {}".format(mstar_min, nout, sum(igal_ok_cat)))
    
    final_sample_galaxies = \
        np.unique(np.concatenate((allgal.data[igal_ok_cat]['id'], id_now)))
    print(" Originally the tree selected:", len(id_now))
    print("Total set length:", len(final_sample_galaxies))
    
    i_gal_ok_final = mtc.match_list_ind(allgal.data['id'], final_sample_galaxies)

    # load GalaxyMaker output again, with return_id this time.
    allgal.data = allgal.data[i_gal_ok_final]
#    allgal = hmo.Halo(base=wdir, nout=nout, is_gal=True, halofinder='HM',
#                      return_id=final_sample_galaxies, load=True)
    print("length of the original allgal: {}".format(len(allgal.data)))

    
    # Associate galaxies with halos. 
    # Mark non-matching galaxies.
    print("Before associate_gal_hal," 
          "length of the original allhal: {}".format(len(allhal.data)))
    allhal = associate_gal_hal(allgal, allhal, plot_check=True, dir_out=wdir)
    print("Now, {}".format(len(allhal.data)))
    

    # load DM ids
#    allhal = hmo.Halo(base=wdir, nout=nout, is_gal=False, halofinder='HM', return_id=allhal.data['id'], load=True)
    # Associate again.
#    allhal = associate_gal_hal(allgal, allhal, plot_check=True, dir_out=wdir)

    # Add idx to catalog.
#    print(gals_in_tree_now['Orig_halo_id'])
#    print(allgal.data['id'])
    for gal in allgal.data:
        ii = np.where(gals_in_tree_now['Orig_halo_id'] == gal['id'])[0]
        if len(ii) > 0:
            gal['idx'] = gals_in_tree_now['id'][ii]
        else:
            gal['idx'] = -1
#    print(gal['idx'])
    allhal.data['idx'] = allgal.data['idx']

    return allgal, allhal
    


def extract_dict_scalar(a):
    import collections
    d = dict()
    for key, val in a.items():
        if isinstance(a[key], collections.Iterable):
            continue
        else:
            d[key] = val
    return d


def save_dict_scalar(cc,f, delim="   "):
    import collections
    keys=[]
    data=[]
    for key in sorted(cc.keys()):
        if isinstance(cc[key], collections.Iterable):
            continue
        else:
            f.write(str(cc[key]) + delim)
    f.write("\n")


    
def worker(gals, hals, out_q, info, inds,
           dump_gal=False,
           reorient=False,
           lambda_method='ellip',
           rscale_lambda=3,
           npix_lambda=10,
           galaxy_plot=False,
           galaxy_plot_dir='galaxy_plot/',
           region_plot=False,
           hydro=False,
           wdir='./',
           **kwargs):
    import galaxy.make_gal as mkg
    from galaxy import galaxy

    nout = info.nout
#    worker_q = Queue()
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
        if hydro:
            gm.cell = load.rd_GM.rd_cell(info.nout, galid, wdir=wdir)
        else:
            gm.cell = None
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


        gal = galaxy.Galaxy(halo = gg,
                            radius_method='eff',
                            info=info)
        good_gal = gal.mk_gal(gm.star, gm.dm, gm.cell,
                              unit_conversion="GM",
                              verbose=True)

        print("GAS MASS in solar mass", gal.meta.mgas)

        galid = gal.meta.id
        print("galaxy {} is made \n".format(galid))

        #do other things
        if not good_gal:
            print(galid, " Not a good galaxy")
            #out_q.put(gal.meta)
        else:
            lambdas = gal.cal_lambda_r_eps(npix_per_reff=npix_lambda,
                           rscale=rscale_lambda, method='ellip', verbose=False) # calculate within 1.0 * reff
            gal.meta.lambda_arr, gal.meta.lambda_arrh, gal.meta.lambda_12kpc= lambdas[0]
            gal.meta.lambda_r,   gal.meta.lambda_rh,   gal.meta.lambda_12kpc = lambdas[1]
            if galaxy_plot:
                gal.plot_gal(fn_save = galaxy_plot_dir + str(nout).zfill(3) \
                                 + "_" + str(galid) + ".png", ioff=True)

            gal.meta.__dict__['idx'] = hals.data['idx'][i]
            gal.meta.__dict__['rhalo'] = hals.data['rvir'][i]

            print("\n \n Good galaxy. ID, IDx", galid, gal.meta.idx)
#            if dump_gal:
#                # save gal before reorientation
#                print(galaxy_plot_dir + str(nout).zfill(3) + "_" + str(galid) + ".h5")
#                mkg.save_gal(gal, galaxy_plot_dir + str(nout).zfill(3) + "_" + str(galid) + ".h5")

            if reorient:
                gal.cal_norm_vec(["star"], dest=[0.,1.,0.])
                gal.cal_rotation_matrix(dest = [0., 1., 0.])
                gal.reorient(dest=[0., 1., 0.], pop_nvec = ['star'],
                             pops=['star', 'dm'], verbose=True)
                lambdas = gal.cal_lambda_r_eps(npix_per_reff=npix_lambda,
                                rscale=rscale_lambda, method='ellip')
# Because lambda is measured for ALL stars,
# B/T must also be measured for ALL stars, not only for bound stars.
#                print(lambdas[0],"d", lambdas[1])
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



def not_all_galaxies_have_hosthalo(allhal, allgal):
    """
       But I don't need this.... 
       No need to know the ids of particles when I aleardy have the particle data.
    """
    newlist = [[]] * len(allhal.data)
    
    for i, this_hal in enumerate(allhal.data['id']):
        ind_idlist = np.where(this_hal == allhal.hal_idlists)[0]
        if len(ind_idlist) > 0:
            newlist[i] = allhal.idlists[ind_idlist]
        else:
            newlist[i] = np.array([])
    	
    allhal.hal_idlists = allhal.data['id']
    # substitute idlists array
    allhal.idlists = newlist


def extract_dm(halo, nout, wdir):
    import utils.sampling as smp
    region = smp.set_region(xc=halo['x'],
                            yc=halo['y'],
                            zc=halo['z'],
                            radius=halo['r'])	
    s = load.sim.Sim(nout=nout, base=wdir, setup=True)
    s.set_ranges(region['ranges'])

    s.add_part(ptypes=["dm id mass pos vel"],
           load=True)

    return s.part.dm

def main(wdir='./',
         ncore=4,
         nout_ini=37,
         nout_end=187,
         read_halo_list=True,
         out_dir="",
         cat_suffix="",
         r_cluster_scale=2.9,
         hydro = False):

    from tree import treemodule
    import os
    import pickle
    import tree.ctutils as ctu
    from queue import Queue
    import pandas as pd
    import multiprocessing as mp
 
    
    #multi = 1 # 
    is_gal = True
    dump_gal = True
    nout_complete = 87 # at least complete up to z = 1
    
    # optional parameters ----------------------------------------------------
    lambda_method = 'ellip' 
    galaxy_plot = True
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
    mstar_min = 5e9
    # Only galaxies above this stellar mass at the final snapshot are considered.
    mstar_min_plot = 5e9
    mk_gal_rscale = 1.1 # unit of Rvir,galaxy
    #r_cluster_scale = 2.9 # maximum radius inside which galaxies are searched for
    npix=800
    rscale_lambda = 3.0 # Reff unit.
    npix_lambda = 5 # per 1Reff
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
        f.write("Rscale lambda calculation : " + str(rscale_lambda))
        f.write("npix per 1Reff for lambda calculation : " + str(npix_lambda))
        f.write("ptypes : \n")
        for i in ptypes:
            f.write("  " + str(i) + "\n")

    #print("------ target galaxy list built")

    masscut_a = 1256366362.16
    masscut_b = -20583566.5218

#%%
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
            mstar_min = get_mstar_min(info.aexp)
        print("Mstar now, min: {:.2e} {:.2e}".format(2 * \
                    (masscut_a * info.aexp + masscut_b), mstar_min))

        # Fix it, no need to load idlist.
        allgal, allhal = get_sample_gal(wdir, nout, info, prg_only_tree, mstar_min)

        # not_all_galaxies_have_hosthalo(allhal, allgal)
        # print("allgal and allgal",allgal.data['id'], allhal.data['id'])

        
###
#  required quantities
#  halo.data[] x, y, z, vx, vy, vz, m, mvir, r, rivr, id, idx, 
#  halo.idlists
#
#
#
#
        nh = len(allgal.data)
        print("Total # galaxies to analyze",nh)
        print("# complete-tree galaxies",sum(allhal.data['idx'] > 0))
        print("# non-complete-tree galaxies",sum(allhal.data['idx'] < 0))

        keywords = dict(rscale = mk_gal_rscale,
                    verbose=verbose,
                    mstar_min=mstar_min
                    ,hydro=hydro)#, dump_gal = dump_gal,
    #                reorient=reorient)

    #   Multiprocessing -----------------------------------------------------------
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
                    lambda_method,
                    rscale_lambda,
                    npix_lambda,
                    galaxy_plot,
                    galaxy_plot_dir,
                    region_plot), kwargs=keywords) for i in range(ncore)]

    
        keywords = dict(rscale = mk_gal_rscale,
                        verbose=verbose,
                      mstar_min=mstar_min)#, dump_gal = dump_gal,
    #                reorient=reorient)
    
        # Run processes
        for p in processes:
            p.start()
        
        # Exit the completed processes
        for p in processes:
            p.join()
                
        
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
        
        catalog = pd.DataFrame(dictout).to_records() 
        with open(fcat, 'wb') as f:
            pickle.dump(catalog, f)
            
        s = 0
        # minimum stellar mass check only for the final snapshot galaxies,
        # No more mstar_min test.
        print("------------------ \n")
        #print("main loop took ", time.time() - t0, "seconds")
    
