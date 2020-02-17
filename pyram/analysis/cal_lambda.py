import numpy as np
import utils.sampling as smp
import collections
from galaxymodule import galaxy
import utils.match as mtc
import tree.ctutils as ctu
import pickle
import load
import tree.halomodule as hmo
import general
from ..load.part import Part
from ..load.hydro import Hydro


def get_mstar_min(aexp):
    masscut_a = 1256366362.16
    masscut_b = -20583566.5218

    return masscut_a * aexp + masscut_b


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
                dddv = dd/2e-4 + dv/150
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


def get_sample_tree(alltrees,
                    info,
                    wdir = './',
                    nout_ini=57,
                    nout_fi = 18,
                    is_gal = True,
                    r_cluster_scale = 2.9,
                    m_halo_min = 5e9,
                    nout_complete = 87):
    """
    Returns a list of galaxies of interest over snapshots

    Todo
    ----
    Add host halo property by runing associate_gal_hal() routine,
    so that all information needed to make a Galaxy instance is
    contained. (2016/06/06)
    """
    import load
    import tree.ctutils as ctu

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


def get_sample_gal(wdir, nout, info, prg_only_tree, mstar_min):
    """
        return a list of galaxies of interest a the given nout.

        Based on the progenitor-only-tree, add galaxies
        inside the zoom region above the mass cut.
    """
    import utils.match as mtc
    import numpy.lib.recfunctions as rf

    gals_in_tree_now = prg_only_tree[prg_only_tree['nout'] == nout]
    id_now = gals_in_tree_now['Orig_halo_id'] # this is Orig_halo_id

    # Galaxies near the cluster
    allhal = hmo.Halo(base=wdir, nout=nout,
                      is_gal=False,
                      halofinder='HM',
                      return_id=False,
                      load=True)
    cluster_now = allhal.data[allhal.data.np.argmax()]

    max_dist_prg = max(np.sqrt(np.square(cluster_now['x'] * 200 - gals_in_tree_now['x']) +
                               np.square(cluster_now['y'] * 200 - gals_in_tree_now['y']) +
                               np.square(cluster_now['z'] * 200 - gals_in_tree_now['z']))) # in Mpc/h

    # unless you have "perfect" trees, do not take the mass of the least massive progentor
    # as the non-tree galaxy mass cut. Sometimes tree marks a ridiculously small galaxy
    # as a progenitor.
    # mstar_min = min(gals_in_tree_now['m'][(gals_in_tree_now['mmp'] == 1) * (gals_in_tree_now['phantom'] == 0)])

    # Galaxy with enough stellar mass
    # -> there is no stellar mass cut when building all_gals_in_trees list.
    allgal = hmo.Halo(base=wdir, nout=nout, is_gal=True, halofinder='HM',
                      return_id=False, load=True)
    allgal.data = rf.append_fields(allgal.data, "tree_root_id", dtypes='i8', data = -1 * np.ones(len(allgal.data)))

    dd_cat = np.sqrt(np.square(cluster_now['x'] - allgal.data['x']) +
                     np.square(cluster_now['y'] - allgal.data['y']) +
                     np.square(cluster_now['z'] - allgal.data['z'])) * 200

    igal_ok_cat = (dd_cat < max_dist_prg) * (allgal.data['m'] > mstar_min)
    print("[get_sample_gal] (catalogue) # galaxies more massive than {:.2e} at nout ="
    " {}".format(mstar_min, nout, sum(igal_ok_cat)))

    final_sample_galaxies = \
        np.unique(np.concatenate((allgal.data[igal_ok_cat]['id'], id_now)))
    print("[get_sample_gal] Originally the tree selected:", len(id_now))
    print("[get_sample_gal] Total set length:", len(final_sample_galaxies))

    i_gal_ok_final = mtc.match_list_ind(allgal.data['id'], final_sample_galaxies)

    # load GalaxyMaker output again, with return_id this time.
    allgal.data = allgal.data[i_gal_ok_final]
    print("[get_sample_gal] length of the original allgal: {}".format(len(allgal.data)))


    # Associate galaxies with halos.
    # Mark non-matching galaxies.
    print("[get_sample_gal] Before associate_gal_hal,"
          "[get_sample_gal] length of the original allhal: {}".format(len(allhal.data)))
    allhal = associate_gal_hal(allgal, allhal, plot_check=False, dir_out=wdir)
    print("[get_sample_gal] Now, {}".format(len(allhal.data)))


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
            gal['idx'] = np.int32(gals_in_tree_now['id'][ii])
            gal["tree_root_id"] = gals_in_tree_now["tree_root_id"][ii]
        else:
            gal['idx'] = np.int32(-1)
            gal["tree_root_id"] = -1
#    print(gal['idx'])
    allhal.data['idx'] = allgal.data['idx']

    return allgal, allhal

def find_host_single(galdata, allhal):
    """
        find the host halo of the given galaxy.

        galdata should have ['x', 'y', 'z'. 'vx', 'vy', 'vz', 'm', 'r']

        INCOMPLETE
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


def load_dm_direct(halo, info, wdir, region=None):
    """
    Load and return DM paticles at the location of the halo.
    But unlike HaloMaker dump files, subhalos may be included.
    """
    if region is None:
        import utils.sampling as smp
        region = smp.set_region(xc=halo['x'],
                                yc=halo['y'],
                                zc=halo['z'],
                                radius=halo['r'])
    part = Part(info, region=region,
                 base=wdir,
                 ptypes=["dm id pos vel"],
                 load=True)
    return part.dm


def load_cell_direct(halo, info, wdir, region=None, rscale=1.0):
    """
    Actually, it load the data from the binary.
    Does not extract from a larger chunk of memory.
    """
    if region is None:
        import utils.sampling as smp
        region = smp.set_region(xc=halo['x'],
                                yc=halo['y'],
                                zc=halo['z'],
                                radius=halo['r']*rscale)

    hh = Hydro(info=info, region=region,
#               base=wdir,
               load=True)

    return hh.cell
