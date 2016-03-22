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

"""
import matplotlib
matplotlib.use("Agg")

import numpy as np
import utils.sampling as smp
import matplotlib.pyplot as plt 

from galaxy import galaxy
import utils.match as mtc
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
    

def associate_gal_hal(allgal, allhal):
    # associate halos with galaxies. 
    def distv(halo, center):
        norm = np.sqrt(np.square(center['vx'] - halo.vx) + 
                       np.square(center['vy'] - halo.vy) + 
                       np.square(center['vz'] - halo.vz))
        
        return norm

    i0=[]
    i1=[]

    newhals = np.recarray(len(allgal.data), dtype=allhal.data.dtype) # equivalent with allhal.data
    
    for i, gal in enumerate(allgal.data):
        dd = dist(allhal.data, gal)
        d_sort = np.argsort(dd)
        if (dd[d_sort[0]] < 0.1 * dd[d_sort[1]]) and (dd[d_sort[0]] < 5e-5):
            gal['hosthalo'] = allhal.data['id'][d_sort[0]]
            i0.append(i)
            newhals[i] = allhal.data[d_sort[0]]
        elif (dd[d_sort[0]] < 0.3 * dd[d_sort[1]]) and (dd[d_sort[0]] < 2.5e-5):
            gal['hosthalo'] = allhal.data['id'][d_sort[0]]
            i0.append(i)
            newhals[i] = allhal.data[d_sort[0]]
        else:
            dv = distv(allhal.data, gal)
            d_nomi = dd < 2e-4 # within 40kpc!! 
            v_nomi = dv < 150
            
            if sum(d_nomi * v_nomi) > 1:
                dddv = dd[d_nomi]/2e-4 + dv[v_nomi]/150
                dddv_sort = np.argsort(dddv)
                gal['hosthalo'] = allhal.data['id'][dddv_sort[0]]
                i0.append(i)
                newhals[i] = allhal.data[dddv_sort[0]]
            else:
                gal['hosthalo'] = -1
                i1.append(i)
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
    return allhal


def join_struct_arrays(arrays):
    sizes = np.array([a.itemsize for a in arrays])
    offsets = np.r_[0, sizes.cumsum()]
    n = len(arrays[0])
    joint = np.empty((n, offsets[-1]), dtype=np.uint8)
    for a, size, offset in zip(arrays, sizes, offsets):
        joint[:,offset:offset+size] = a.view(np.uint8).reshape(n,size)
    dtype = sum((a.dtype.descr for a in arrays), [])
    return joint.ravel().view(dtype)

def augment_tree(treedata, base, is_gal=False):
    """
        Add more quantities to existing tree data. 
        
        Consistent tree (built with HM/GM output) does not provide much detail of halos/galaxies.
        I need to add some more information from original HM/GM output.
    """
    
    dtype_new_quantities = [('np', '<i4'), ('id', '<i4'), ('m', '<f4'), ('mvir', '<f4'),
                            ('r', '<f4'), ('rvir', '<f4'), ('tvir', '<f4'), ('cvel', '<f4'),
                            ('x', '<f4'), ('y', '<f4'), ('z', '<f4'),
                            ('vx', '<f4'), ('vy', '<f4'), ('vz', '<f4'),
                            ('ax', '<f4'), ('ay', '<f4'), ('az', '<f4'),
                            ('sp', '<f4')]
    if is_gal:
        [dtype_new_quantities.append(i) for i in [('sig', '<f4'), ('sigbulge', '<f4'), ('mbulge', '<f4')]]
           
    New_arr = np.zeros(len(treedata), dtype=dtype_new_quantities)
    import tree.halomodule as hmo
    for nout in np.unique(treedata['nout']):
        # nout and Orig_halo_id are required.
        gal_org = hmo.Halo(base=wdir, nout=nout, halofinder='HM', load=True, is_gal=is_gal)
        # Before we start, remove unnecessary coulmns
        dtype_names = [field[0] for field in dtype_new_quantities]
        gal_org = gal_org.data[dtype_names]
        
        ind_tree_this_nout = np.where(treedata['nout'] == nout)[0]
        ok_gals = treedata['Orig_halo_id'][ind_tree_this_nout]
        
        # Galaxies are from a snapshot. Galaxy ID list must be a unique set.
        assert len(ok_gals) == len(np.unique(ok_gals))
        
        ind_org_gals = [np.where(gal_org['id'] == gal)[0] for gal in ok_gals]
        
        for i, ind in enumerate(ind_org_gals):
            assert sum(New_arr[ind_tree_this_nout[i]]) == 0. # array must be empty
            New_arr[ind_tree_this_nout[i]] = gal_org[ind]
 
    # Drop duplicate fields
    #["id", "mvir", "rvir", "x", "y", "z", "vx", "vy", "vz"]
    keep_fields = ["np", "m", "r", "tvir", "cvel"]
    if is_gal:
        [keep_fields.append(i) for i in ['sig', 'sigbulge', 'mbulge']]
        
    return join_struct_arrays([treedata, New_arr[keep_fields]])


def distance_to(xc, xx):
    import numpy as np
    return np.sqrt([(xc[0] - xx[0])**2 + (xc[1] - xx[1])**2 + (xc[2] - xx[2])**2])[0]

def extract_halos_within(halos, i_center, info, dist_in_mpc=1.0):

    xc = halos['x'][i_center]
    yc = halos['y'][i_center]
    zc = halos['z'][i_center]

    xx = halos['x']
    yy = halos['y']
    zz = halos['z']

    dd = np.multiply(distance_to([xc,yc,zc], [xx,yy,zz]), info.pboxsize)

    return (dd < (dist_in_mpc))

def all_gals(treedata, final_gals, nout_ini=None, nout_fi=None):
    if nout_ini == None:
        nout_ini = min(treedata['nout'])
    if nout_fi == None:
        nout_fi = max(treedata['nout'])
    
    all_gals_at_nouts = []
    for inout, nout in enumerate(nouts):
        all_gals_this_nout = []
        tree_now = treedata[np.where(treedata['nout'] == nout)]

        for finalgal in final_gals:
            i_gals_include = np.where(tree_now['tree_root_id'] == finalgal)[0]
            [all_gals_this_nout.append(gal) for gal in tree_now['id'][i_gals_include]]
            
        all_gals_at_nouts.append(all_gals_this_nout)
        
    return all_gals_at_nouts



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
    h.data['x'] = tree_element['x']/cboxsize
    h.data['y'] = tree_element['y'] /cboxsize
    h.data['z'] = tree_element['z'] /cboxsize # Mpc/h -> code unit
    h.data['vx'] = tree_element['vx']
    h.data['vy'] = tree_element['vy']
    h.data['vz'] = tree_element['vz']
    h.data['r'] = tree_element['r'] # already in code unit
    h.data['rvir'] = tree_element['rvir'] / (cboxsize * 1000) # kpc/h -> code unit
    h.data['id'] = tree_element['Orig_halo_id']
    h.data['idx'] = tree_element['id']
    h.aexp = tree_element['aexp']
    
    return h


def show_target_summary(allgal):
    print("OK galaxies")
    ok = allgal.data['idx'] > 0
    print("IDX", allgal.data['idx'][ok])
    print("mass", allgal.data['m'][ok])
    
    print("/n bad galaxies")
    print("IDX", allgal.data['idx'][~ok])
    print("mass", allgal.data['m'][~ok])




##########################################################################
def mk_gal(galdata, halodata, info, star_id=None, dm_id=None,
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
    
    """
    plt.close()
    import make_gal as mkg

    if star_id is None:
        star, dm, cell = mkg.extract_data_old(halodata, rscale=rscale)
        # Only when not using GalaxyMaker output, "not enough stars" error can occur.
        if sum(star['m']) * info.msun < mstar_min:
            print("(1)Not enough stars: {:.2f} Msun".format(sum(star['m']) * info.msun))
            print("Aborting... \n", " Not a good galaxy")
            out_q.put(gal_out)
        return 
    else:
        xc_tmp0 = halodata['x']
        yc_tmp0 = halodata['y']
        zc_tmp0 = halodata['z']
 
        rr_tmp0 = min([halodata['r'] * rscale, 0.0002])
        rr_tmp0 = max([rr_tmp0, 0.000025])
 
        if cell_all is not None:
            ind_c = ((cell_all['x'] - xc_tmp0)**2 + (cell_all['y'] - yc_tmp0)**2 +
                    (cell_all['z'] - zc_tmp0)**2) < rr_tmp0**2
            ind_c = ind_c * (cell_all['var0'] > min_gas_density)
 
            ind_c = np.where((cell_all['x'] - xc_tmp0)**2 + (cell_all['y'] - yc_tmp0)**2
                            + (cell_all['z'] - zc_tmp0)**2 < rr_tmp0**2)[0]
        else:
            return None
 
        cell = cell_all[ind_c]
        star = star_all[mtc.match_list_ind(star_all['id'], star_id)]
        # There are a few fake halos.
        # Assume a halo should have at least 64 dm particles.
        if len(dm_id) > 64:
            dm = dm_all[mtc.match_list_ind(dm_all['id'], dm_id)]
        else:
            dm = dm_all[np.where((dm_all['x'] - halodata['x'])**2 +
                             (dm_all['y'] - halodata['y'])**2 +
                             (dm_all['z'] - halodata['z'])**2 < halodata['r']**2)[0]]

               
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
    
        rgal = region['radius'] * s.info.pboxsize * 1e3

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
    good_gal = gal.mk_gal(star, dm, cell,
                        mstar_min=mstar_min,
               rscale=rscale, verbose=verbose, method_com=method_com)
    return gal, good_gal

    
def set_affinity_on_worker():
    import os
    """When a new worker process is created, the affinity is set to all CPUs"""
    print("I'm the process %d, setting affinity to all CPUs." % os.getpid())
    os.system("taskset -p 0xff %d" % os.getpid())    
    

def worker(gals, hals, out_q, info, inds,
           dump_gal=False,
           reorient=False,
           lambda_method='ellip', rscale_lambda=3,npix_lambda=10,
           galaxy_plot=False, galaxy_plot_dir='galaxy_plot/',
           **kwargs):
    import make_gal as mkg

#    worker_q = Queue()
    if type(inds) == int:
        inds = [inds]
    for i in inds:
        print("This is {}-th galaxy".format(i))
        gal, good_gal = mk_gal(gals.data[i],
                               hals.data[i],
                               s.info,
                               star_id=gals.idlists[i],
                               dm_id=hals.idlists[i],
                               **kwargs)
        print("galaxy is made \n")
        gal.meta = {"id":0, "idx":hals.data['idx'][i], "xc":0.0, "yc":0.0, "zc":0.0,
                    "vx":0.0, "vy":0.0, "vz":0.0,
                    "mstar":0.0, "nstar":0.0, "mgas":0.0,
                    "lambda_arr":[], "lambda_r":0, "rgal":0,
                    "rhalo":hals.data['rvir'][i], "boxtokpc":info.pboxsize*1e3,
                    "b2t":0.0, "eps":0.0, "mjr":0.0}

        #do other things
        if not good_gal:
            print(gal.id, " Not a good galaxy")
            #out_q.put(gal.meta)
        else:
            if reorient:
                gal.cal_norm_vec(["star"])
                gal.cal_rotation_matrix(dest = [0., 0., 1.]
                gal.reorient(dest=[0., 0., 1], pop_nvec = ['star'], verbose=False)
 
            # Calculate lambda_r ---------------------------------------------------
            gal.cal_lambda_r_eps(npix=npix_lambda, method=lambda_method,
                                 rscale=rscale_lambda,
                                 galaxy_plot_dir=galaxy_plot_dir) 
 
            print("\n \n Good galaxy. ID, IDx", gal.id, gal.meta['idx'])
            # Save to catalog ----------------------------------------------------
            gal.meta['mstar'] = gal.mstar
            gal.meta['mgas'] = gal.mgas
            gal.meta['nstar'] = gal.nstar
            gal.meta['id'] = gal.id
            gal.meta['xc'] = gal.xc * info.pboxsize
            gal.meta['yc'] = gal.yc * info.pboxsize
            gal.meta['zc'] = gal.zc * info.pboxsize
            gal.meta['vx'] = gal.vxc * info.kms
            gal.meta['vy'] = gal.vyc * info.kms
            gal.meta['vz'] = gal.vzc * info.kms        
            gal.meta['lambda_arr'] = gal.lambda_arr
            gal.meta['lambda_r'] = gal.lambda_r
            gal.meta['rgal'] = gal.reff# * info.pboxsize * 1000.0 # in kpc  
            gal.meta['eps'] = gal.eps # Eps = 1-b/a
            gal.meta['pa'] = gal.pa
            gal.meta['mjr'] = gal.mrj
            out_q.put(gal.meta)
 
            if dump_gal:
                print(galaxy_plot_dir + str(nout).zfill(3) + "_" + str(gal.id) + ".h5")
                mkg.save_gal(gal, galaxy_plot_dir + str(nout).zfill(3) + "_" + str(gal.id) + ".h5")
 
            if galaxy_plot:
                gal.plot_gal(fn_save = galaxy_plot_dir + str(nout).zfill(3) \
                              + "_" + str(gal.id) + ".png", ioff=True)
#    for i in inds:
#        out_q.put(worker_q.get())

###########################################################################
"""
The processing pool needs to be instantiated in the main 
thread of execution. 
"""
from queue import Queue
import pandas as pd
import pickle
import multiprocessing as mp    
import load
from tree import treemodule
import tree.ctutils as ctu
import tree.halomodule as hmo 
import os

multi = 1 # 
hydro = True
is_gal = True
dump_gal = True

# optional parameters ----------------------------------------------------
lambda_plot = False 
lambda_method = 'ellip' 
galaxy_plot = False
reorient = False

wdir = input("Working directory \n")
#wdir = './'
#ncore=1
#nout_ini=186
#nout_end=186

if multi: ncore = int(input("How many cores? \n"))
nout_ini = input("First snapshot: (default = 37 , z=3) \n")
nout_end = input("Last snapshot: (default = 187, z=0) \n")
read_halo_list = input("Load existing galaxy list? (y/n)")
out_dir = input("output directory: (default = out_001++)")
cat_suffix = input("suffix to the catalog file? (default = '')")
dir_cat = out_dir

if read_halo_list == 'y':
    read_halo_list = True
else:
    read_halo_list = False

if out_dir == "":
    dir_increment = 0
    out_dir = "out_" + str(dir_increment).zfill(3) + '/'
    while not os.path.isdir(out_dir):
        dir_increment += 1
        out_dir = "out_" + str(dir_increment).zfill(3) + '/'

if out_dir[-1] != '/':
    out_dir = out_dir + '/' 
    # need a tailing slash to be joined with sub-directories

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

out_base = wdir + out_dir

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
#nouts = range(nout_fi, nout_ini -1, -1) 
nouts = [187, 121, 87, 54, 37]
nouts_dump_gal = [187, 154, 130, 121, 112,  98,  87,  67,  54,  37,  20]
try: nouts_dump_gal
except NameError: nouts_dump_gal = None
    

#----------------------------------------------------------------------
mstar_min = 5e9
# Only galaxies above this stellar mass at the final snapshot are considered.
mstar_min_plot = 5e9
mk_gal_rscale = 1.1 # unit of Rvir,galaxy
r_cluster_scale = 2.5 # maximum radius inside which galaxies are searched for
npix=800
rscale_lambda = 3.0 # Reff unit.
npix_lambda = 10 # per 1Reff
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

# pickle tree data -----------------------------------------------------
tt = alltrees.data
tt_final = tt[tt['nout'] == nout_fi]

if read_halo_list:
    try:
        all_gals_in_trees = pickle.load(open(wdir + "all_gals_in_trees.pickle", 'rb'))
        final_gals_idx = all_gals_in_trees[-1]
        print("loading pickled halo list done:")
        print(len(final_gals_idx), "halos left")
        ngals = len(final_gals_idx)
        read_halo_list = True
    except:
        read_halo_list = False

if not read_halo_list:
    info = load.info.Info(nout=nout_fi, base=wdir, load=True)
    print(" It's new run starting from the final snapshot.")
    print(" Loading halo list and selecting only massive halo/galaxies"
            " close enough to the central cluster")
    hhal = hmo.Halo(base=wdir, nout=nout_fi, halofinder='HM', info=info, load=True, is_gal=False)
    i_center = np.where(hhal.data['np'] == max(hhal.data['np']))[0]
    r_cluster = hhal.data['rvir'][i_center] * info.pboxsize

    hh = hmo.Halo(base=wdir, nout=nout_fi, halofinder='HM', info=info, load=True, is_gal=is_gal)
    i_center = np.where(hh.data['np'] == max(hh.data['np']))[0]
    i_satellites = extract_halos_within(hh.data, i_center, info, dist_in_mpc = r_cluster * r_cluster_scale)
    print("Total {0} galaxies \n{1} galaxies are selected".format(
          len(i_satellites),sum(i_satellites)))
    
    # halos found inside the cluster and have complete tree back to nout_ini
    large_enugh = hh.data['mvir'] > m_halo_min
    halo_list = hh.data['id'][i_satellites * large_enugh]
    final_ids = ctu.check_tree_complete(tt, 87, nout_fi, halo_list, idx=False) # 87: z = 1

    final_gals_idx = [tt_final['id'][tt_final['Orig_halo_id'] == final_gal] for final_gal in final_ids]
    print(len(final_gals_idx), "halos left")
    ngals = len(final_gals_idx)
    # Search for all galaxies that listed in the trees of final_gals
    all_gals_in_trees = all_gals(tt, final_gals_idx)

if not read_halo_list:    
    print(" Pickling halo list...")
    pickle.dump(all_gals_in_trees, open( wdir + "all_gals_in_trees.pickle", "wb" ))
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


#import time
#os.system("taskset -p 0xff %d" % os.getpid())
# http://stackoverflow.com/questions/15414027/multiprocessing-pool-makes-numpy-matrix-multiplication-slower

#print(all_gals_in_trees)
masscut_a = 1256366362.16
masscut_b = -20583566.5218

#import utils.match as mtc

for inout, nout in enumerate(nouts):
    print(inout, nout, nout_fi)

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
    
#    mstar_min = min([masscut_a * info.aexp + masscut_b, mstar_min]) * 2
    print("Mstar min: {:.2e}".format(mstar_min))

    t_now = alltrees.data[alltrees.data['nout'] == nout]
    idxs_now = all_gals_in_trees[inout]
#    print(idxs_now)

    h = halo_from_tree(t_now[mtc.match_list_ind(t_now['id'], np.array(idxs_now))], info)
#    print(h.data,h.data.dtype)


    allhal = hmo.Halo(base=wdir, nout=nout, is_gal=False, halofinder='HM', return_id=False, load=True)
    clu = allhal.data[allhal.data.np.argmax()]

    # Galaxy with enough stellar mass
    allgal = hmo.Halo(base=wdir, nout=nout, is_gal=True, halofinder='HM', return_id=False, load=True)
    i_gal_ok = allgal.data['m'] > mstar_min
#    print("i_gal_ok", i_gal_ok)
    i_gal_near = find_halo_near(allgal.data, clu, rscale=2.5)
    id_gal_ok = allgal.data['id'][i_gal_ok * i_gal_near]
    print("# id_gal_ok", len(id_gal_ok))


    # get the stellar particle IDs of OK galaxies.
    allgal = hmo.Halo(base=wdir, nout=nout, is_gal=True, halofinder='HM', return_id=True)
    allgal.set_return_id_list(id_gal_ok)
    allgal.load()
    allgal.data = allgal.data[i_gal_ok * i_gal_near] # only good gals

    # associate halos with galaxies. #########################################
    # Virial radius is required to measure gas mass associated with the galaxy

    allhal.data = allhal.data[find_halo_near(allhal.data, clu, rscale=2.5)] # only good hals

    # new halo catalog with 1:1 correlation to gals.
    allhal = associate_gal_hal(allgal, allhal) 
#    print("allhal.data", allhal.data['id'])
    
    tmp = hmo.Halo(base=wdir, nout=nout, is_gal=False, halofinder='HM', return_id=True, load=False)
    tmp.set_return_id_list(allhal.data['id'])
    tmp.load()

    # move tmp.idlist to allhal with re-ordering.
    # Now idlist of allhal and allgal are corresponding to each other.
    allhal.idlists=[]
    for thisid in allhal.data['id']:
        if thisid < 0:
            allhal.idlists.append([])
        else:
            allhal.idlists.append(tmp.idlists[np.where(tmp.hal_idlists == thisid)[0]])
    
    for gal in allgal.data:
        ii =  np.where(h.data['id'] == gal['id'])[0]
        if len(ii) > 0:
            gal['idx'] = h.data['idx'][ii]
        else:
            gal['idx'] = -1

    allhal.data['idx'] = allgal.data['idx']

    print("idxxxxxx", allgal.data['idx'], allgal.data['id'])
#    show_target_summary(allgal)
    
    nh = len(allgal.data)
    print("Total # galaxies to analyze",nh)
    print("# complete-tree galaxies",sum(allhal.data['idx'] > 0))
    print("# non-complete-tree galaxies",sum(allhal.data['idx'] < 0))
    
    #Load simulation data##########################################################
    s = load.sim.Sim(nout=nout, base=wdir, setup=True)
    s.add_part(ptypes, load=True, fortran=True)
    #assert s.part.nstar > 0, "Not enough stellar particles in given cpus"
    if hydro:
        s.add_hydro(load=True, lmax=lmax)
        cell_all = s.hydro.cell
    else:
        cell_all = None

    star_all = s.part.star
    dm_all = s.part.dm   

    keywords = dict(rscale = mk_gal_rscale,
                verbose=False,
                mstar_min=mstar_min)#, dump_gal = dump_gal,
#                reorient=reorient)

    if multi == 1:
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
                    s.info, inds[i],
                    dump_gal,
                    reorient,
                    lambda_method,
                    rscale_lambda,
                    npix_lambda,
                    galaxy_plot,
                    galaxy_plot_dir), kwargs=keywords) for i in range(ncore)]

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
            pool.apply_async(mk_gal, args=(allhal.data[i], out_q,
                              s.info, i), kwds=keywords)
        pool.close()
        pool.join()
    else:
        for i in range(nh):
            out_q = Queue()
            mk_gal(allhal.data[i], out_q, s.info, i, **keywords)
    
    print("----------Done---------")

#%%    
    dictout=[]
    try:
        if not os.path.isdir(dir_out):
            os.mkdir(dir_out)
        f = open(dir_out + 'galaxies' + snout + '.txt', 'w')
    except:
        print("No filename is given.\n ")
    
    f.write(" nout    ID     IDx      x          y       z[Mpc]    vx      vy     vz[km/s]")
    f.write("    Reff[kpc]     Mstar    Mgas[Msun]  Rhalo[kpc]  boxtokpc  final_ID \n")
    
    for i in range(nh):
        try:
            dd =  out_q.get(timeout=0.1)
            if dd['id'] == 0:
                continue
            f.write("{:<4}   {:<4}  {:<6}  {:.5f}  {:.5f}  {:.5f}".format(nout,
                    dd['id'], dd['idx'], dd['xc'],dd['yc'],dd['zc']))
            f.write("  {:.3f}  {:.3f}  {:.3f}".format(dd['vx'],dd['vy'],dd['vz']))
            f.write("  {:.6f}  {:.0f}  {:.0f}".format(dd['rgal'],dd['mstar'], dd['mgas']))
            f.write("  {:.5f}  {:.5f}  \n".format(dd['rhalo'],dd['boxtokpc']))
            dictout.append(dd)
        except:
            continue
    
    f.close()    
    print("Text file written")

    catalog = pd.DataFrame(dictout).to_records()    

    with open(fcat, 'wb') as f:
        pickle.dump(catalog, f)
        
    star_all = 0
    dm_all = 0
    cell_all = 0
    s = 0
    # minimum stellar mass check only for the final snapshot galaxies,
    # No more mstar_min test.
    print("------------------")
    #print("main loop took ", time.time() - t0, "seconds")

