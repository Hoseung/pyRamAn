# -*- coding: utf-8 -*-
"""
    Created on Mon May 30 15:37:03 2016

    @author: hoseung

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
#import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from ..utils import sampling as smp
from . import galaxy
from ..utils import match as mtc
from ..tree import ctutils as ctu
from ..tree import halomodule as hmo

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
        norm = np.sqrt(np.square(center['vx'] - halo["vx"]) +
                       np.square(center['vy'] - halo["vy"]) +
                       np.square(center['vz'] - halo["vz"]))

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
        fig.set_size_inches(10,10)
        npix = 1200
        pp.pp_halo(allhal, npix, axes=ax, verbose=False, colors=10, rscale=0.3, name=True)
        pp.pp_halo(allgal, npix, axes=ax, verbose=False, colors=200, rscale=0.3)
        ax.set_xlim(-50, npix+50)
        ax.set_ylim(-50, npix+50)
        ax.set_aspect('equal')
        plt.savefig(dir_out + "associate_gal_hal.pdf")
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
    print("Total {} galaxies \n {} galaxies are within {} times the cluster virial radius, {:.2f} Mpc (No mass cut)".format(
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
    print(" {} galaxies have complete tree up to nout = {}, above mass cut {:.2e}".format(ngals, nout_complete, m_halo_min))
    # Search for all galaxies that are listed in the trees of final_gals
    prg_only_tree = all_gals(td, final_gals_idx, nout_fi=nout_fi, nout_ini=nout_ini)
    #all_gals_in_trees = prg_only_tree['']

    # idx list of galaxies at the current nout (from all_gals_in_trees list)
    #t_now = prg_only_tree[prg_only_tree['nout'] == nout]
    #idxs_tree_now = all_gals_in_trees[inout]

    return prg_only_tree

#%%
def get_sample_gal(wdir, nout, info, prg_only_tree, mstar_min, verbose=False):
    #gals_in_tree_now = halo_from_tree(t_now[mtc.match_list_ind(t_now['id'], np.array(idxs_tree_now))], info)
    gals_in_tree_now = prg_only_tree[prg_only_tree['nout'] == nout]
    id_now = gals_in_tree_now['Orig_halo_id'] # this is Orig_halo_id

    # Galaxies near the cluster
    allhal = hmo.Halo(base=wdir, nout=nout, is_gal=False, halofinder='HM', return_id=False, load=True)
    cluster_now = allhal.data[allhal.data.np.argmax()]

    max_dist_prg = max(np.sqrt(np.square(cluster_now['x'] * 200 - gals_in_tree_now['x']) +
                               np.square(cluster_now['y'] * 200 - gals_in_tree_now['y']) +
                               np.square(cluster_now['z'] * 200 - gals_in_tree_now['z']))) # in Mpc/h

    if verbose: print(nout, max_dist_prg)

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
    if verbose:
        print("(catalogue) # galaxies more massive than {:.2e} at nout ="
        " {}".format(mstar_min, nout, sum(igal_ok_cat)))

    final_sample_galaxies = \
        np.unique(np.concatenate((allgal.data[igal_ok_cat]['id'], id_now)))
    if verbose:
        print(" Originally the tree selected:", len(id_now))
        print("Total set length:", len(final_sample_galaxies))

    # load GalaxyMaker output again, with return_id this time.
    allgal = hmo.Halo(base=wdir, nout=nout, is_gal=True, halofinder='HM',
                      return_id=final_sample_galaxies, load=True)
    if verbose: print("length of the original allgal: {}".format(len(allgal.data)))


    # Associate galaxies with halos.
    # Mark non-matching galaxies.
    if verbose:
        print("Before associate_gal_hal,"
            "length of the original allhal: {}".format(len(allhal.data)))
    allhal = associate_gal_hal(allgal, allhal, plot_check=True, dir_out=wdir)
    # there

    if verbose: print("Now, {}".format(len(allhal.data)))

    #print()
    # load DM ids
    allhal = hmo.Halo(base=wdir, nout=nout, is_gal=False, halofinder='HM',
                      return_id=allhal.data['id'], load=True)
    # Associate again.
    allhal = associate_gal_hal(allgal, allhal, plot_check=True, dir_out=wdir)

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


##########################################################################
def mk_gal(galdata, halodata, info, star_all, cell_all, dm_all,
           star_id=None, dm_id=None,
           save=False, rscale=1.1, verbose=False, galaxy_plot_dir='./',
           npix=400,
           region_plot=False,
           method_com=2, mstar_min=5e9, dump_gal=False,
           min_gas_density = 1e-1):
    """
        Direct plot,
        Create galaxy,
        Calculate lambda_r (using Cappellari 2003)
        Draw ed map of galaxy.

        minimum gas densiy = mean gas density of the universe?
        Galaxy center = galdata['x'], ['y'], ['z']. Not halo center.

        To do:
        extract cells by density threshold, not radial cut.
        Gas extent varies a lot, and there no good single choice of radius cut.

    """
    print("mk_gal ID", galdata['id'])
    plt.close()
    import galaxy.make_gal as mkg

    if star_id is None:
        star, dm, cell = mkg.extract_data_old(halodata, rscale=rscale)
        # Only when not using GalaxyMaker output, "not enough stars" error can occur.
        if sum(star['m']) * info.msun < mstar_min:
            print("(1) Not enough stars: {:.2f} Msun".format(sum(star['m']) * info.msun))
            print("Aborting... \n", " Not a good galaxy")
            return -1, False
        return
    else:
        xc_tmp = galdata['x'] #
        yc_tmp = galdata['y']
        zc_tmp = galdata['z']
        rr_tmp0 = galdata['r']
# If memory access is the bottle neck,
# There is no reason the following part need to be called
# by individual process.
# put it in the main body.

        if cell_all is not None:
            ind_c = ((cell_all['x'] - xc_tmp)**2 + (cell_all['y'] - yc_tmp)**2 +
                    (cell_all['z'] - zc_tmp)**2) < rr_tmp0**2
            ind_c = ind_c * (cell_all['var0'] > min_gas_density)

            ind_c = np.where((cell_all['x'] - xc_tmp)**2 + (cell_all['y'] - yc_tmp)**2
                            + (cell_all['z'] - zc_tmp)**2 < rr_tmp0**2)[0]
            cell = cell_all[ind_c]
        else:
            cell = None
            #return None
        star = star_all[mtc.match_list_ind(star_all['id'], star_id)]
        # There are a few fake halos.
        # Assume a halo should have at least 64 dm particles.
        if len(dm_id) > 64:
            dm = dm_all[mtc.match_list_ind(dm_all['id'], dm_id)]
        else:
            dm = dm_all[np.where((dm_all['x'] - xc_tmp)**2 +
                             (dm_all['y'] - yc_tmp)**2 +
                             (dm_all['z'] - zc_tmp)**2 < halodata['r']**2)[0]]


    # Direct plot ---------------------------------------------------------
    if region_plot:
        import draw
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

        rgal = region.radius * info.pboxsize * 1e3

        ax = plt.gca()
        ax.set_xlabel("position [kpc]")
        ax.set_xticks(np.linspace(0,npix,5))
        xticks = ["{:.2f}".format(x) for x in np.linspace(-rgal, rgal, num=5)]
        ax.set_xticklabels(xticks)
        ax.set_ylabel("position [kpc]")
        ax.set_yticks(np.linspace(0,npix,5))
        yticks = ["{:.2f}".format(y) for y in np.linspace(-rgal, rgal, num=5)]
        ax.set_yticklabels(yticks)

        plt.savefig(galaxy_plot_dir+"2dmap_"+str(halodata['id']).zfill(5)+'.png', dpi=144)
        plt.close()


    #Create galaxy ----------------------------------------------
    gal = galaxy.Galaxy(galdata, radius_method='eff', info=info)
    print("gal instance", gal.meta)
    good_gal = gal.mk_gal(star, dm, cell,
                        mstar_min=mstar_min,
               rscale=rscale, verbose=verbose, method_com=method_com)
    return gal, good_gal


def extract_dict_scalar(a):
    import collections
    d = dict()
    for key, val in a.items():
        if isinstance(a[key], collections.Iterable):
            continue
        else:
            d[key] = val
    return d


def save_dict_scalar(cc, f, delim="   "):
    import collections
    for key in sorted(cc.keys()):
        if isinstance(cc[key], collections.Iterable):
            continue
        else:
            f.write(str(cc[key]) + delim)
    f.write("\n")



def worker(gals, hals, out_q, info, inds,
           star_all, cell_all, dm_all,
           dump_gal=False,
           reorient=False,
           lambda_method='ellip', rscale_lambda=3,npix_lambda=10,
           galaxy_plot=False,
           galaxy_plot_dir='galaxy_plot/',
           region_plot=False,
           **kwargs):
    import galaxy.make_gal as mkg

    if type(inds) == int:
        inds = [inds]
    for i in inds:
        print("This is {}-th galaxy".format(i))
        print("halo: {}".format(hals.data['id'][i]), )
        print("gal: {}".format(gals.data['id'][i]), )
#        print("number of particles", len(gals.idlists[i]), gals.data['np'][i])
        gal, good_gal = mk_gal(gals.data[i],
                               hals.data[i],
                               info,
                               star_all,
                               cell_all,
                               dm_all,
                               star_id=gals.idlists[i],
                               dm_id=hals.idlists[i],
                               region_plot=region_plot,
                               **kwargs)
        nout = info.nout
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
            if dump_gal:
                # save gal before reorientation
                print(galaxy_plot_dir + str(nout).zfill(3) + "_" + str(galid) + ".h5")
                mkg.save_gal(gal, galaxy_plot_dir + str(nout).zfill(3) + "_" + str(galid) + ".h5")

            if reorient:
                gal.cal_norm_vec(["star"], dest=[0.,1.,0.])
                gal.cal_rotation_matrix(dest = [0., 1., 0.])
                gal.reorient(dest=[0., 1., 0.], pop_nvec = ['star'],
                             pops=['star', 'dm'], verbose=False)
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
