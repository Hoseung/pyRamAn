"""
    MODIFICATIONS
    
    2015.12.08
        Without tree version.
        Only use mass cut.
    
"""

import numpy as np
import utils.sampling as smp

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


def extract_data(halo, rscale=1.25):
    xc_tmp0 = halo['x']
    yc_tmp0 = halo['y']
    zc_tmp0 = halo['z']
    
    rr_tmp0 = min([halo['r'] * rscale, 0.0002]) 
    # it's galaxy, not halo. 
    # 'r' is the furthest stellar particle
    # 'rvir' doesn't have robust physical meaning 
    # because galaxies are hardly virialized systems.
    
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


# In[6]:

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


def mk_gal(halodata, out_q, info, i,
           save=False, rscale=1.1, verbose=False, galaxy_plot_dir='./',
           rscale_lambda=2.0, npix_lambda=50, npix=400, galaxy_plot=False,
           method_com=2, mstar_min=5e9, dump_gal=False):
    """
    Direct plot,
    Create galaxy, 
    Calculate lambda_r (using Cappellari 2003)
    Draw ed map of galaxy.
    
    """

    print("This is {}-th halo".format(i))
    from galaxy import galaxy


    gal_out = {"id":0, "xc":0.0, "yc":0.0, "zc":0.0,
               "vx":0.0, "vy":0.0, "vz":0.0,
               "mstar":0.0, "nstar":0.0, "mgas":0.0,
               "lambda_arr":[], "lambda_r":0, "rgal":0,
               "rhalo":halodata['rvir'], "boxtokpc":info.pboxsize*1000}
               
    star, dm, cell = extract_data(h.data[i], rscale=rscale)

    if sum(star['m']) * info.msun < mstar_min:
        print("(1)Not enough stars: {:.2f} Msun".format(sum(star['m']) * info.msun))
        print("Aborting... \n")
        print(" Not a good galaxy")
        out_q.put(gal_out)
        return
               
               
    # Direct plot ---------------------------------------------------------                                
    if galaxy_plot:
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
        xticks = ["{:.2f}".format(x) for x in np.linspace(-rgal, rgal, num=5)]
        ax.set_xticklabels(xticks)
        ax.set_ylabel("position [kpc]")
        ax.set_yticks(np.linspace(0,npix,5))
        yticks = ["{:.2f}".format(y) for y in np.linspace(-rgal, rgal, num=5)]
        ax.set_yticklabels(yticks)
        
        plt.savefig(galaxy_plot_dir+"2dmap_"+str(halodata['id']).zfill(5)+'.png', dpi=144)
        plt.close()

    #Create galaxy ----------------------------------------------
    gal = galaxy.Galaxy(halodata, radius_method='eff', info=info)
    good_gal = gal.mk_gal(star, dm, cell,
                        mstar_min=mstar_min,
               rscale=rscale, verbose=verbose, method_com=method_com)

    #-----------------------------------------------------------------------    
    def get_metadata(clazz):
            """
                Out of all attributes of a galaxy instance, leave only data.
            """
            return {name: attr for name, attr in clazz.__dict__.items()
                    if not name.startswith("__") 
                    and not callable(attr)
                    and not type(attr) is staticmethod}
        
    def get_metadata2(adict):
        return {name: attr for name, attr in adict.items()
                if not isinstance(attr, (np.ndarray, np.recarray, dict, list, str))}
    
    def save_gal(galaxy, filename, verbose=False):    
        import h5py as hdf
        
        outfile = hdf.File(filename, 'w')
        
        # Store metadata using HDF5 attributes
        attrs = get_metadata(galaxy)
        attrs = get_metadata2(attrs)
        #outfile.attrs.create("all", attrs)
        for name, atr in attrs.items():
            if atr != None:
                outfile.attrs.create(name, atr)
            #outfile.attrs[name] = atr
        
        # Store data under /galaxy with direct assignment
        if hasattr(galaxy, 'star'):
            if verbose: print("Saving star")
            star = outfile.create_group("star")
            for field in galaxy.star.dtype.names:
                star.create_dataset(field, data=galaxy.star[field])
        
        if hasattr(galaxy, 'dm'):
            if verbose: print("Saving DM")
            dm = outfile.create_group("dm")
            for field in galaxy.dm.dtype.names:
                dm.create_dataset(field, data=galaxy.dm[field])
            
        if hasattr(galaxy, 'cell'):
            if verbose: print("Saving cell")
            gas = outfile.create_group("cell")
            for field in galaxy.cell.dtype.names:
                gas.create_dataset(field, data=galaxy.cell[field])
        
        outfile.close()

    if not good_gal:
        print(gal.id, " Not a good galaxy")
        out_q.put(gal_out)
    else:
        # Save to catalog -------------------------------------------------------
        gal.cal_lambda_r(npix=npix_lambda, method=1, rscale=rscale_lambda) # calculate within 1.0 * reff    
        
        if dump_gal:
            #gal.save_gal(galaxy_plot_dir + str(nout).zfill(3) \
            save_gal(gal, galaxy_plot_dir + str(nout).zfill(3) + "_" + str(gal.id) + ".h5")        
        
        # Calculate lambda_r ---------------------------------------------------

        gal.plot_gal(fn_save = galaxy_plot_dir + str(nout).zfill(3) \
                             + "_" + str(gal.id) + ".png", ioff=True)


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

    
def set_affinity_on_worker():
    import os
    """When a new worker process is created, the affinity is set to all CPUs"""
    print("I'm the process %d, setting affinity to all CPUs." % os.getpid())
    os.system("taskset -p 0xff %d" % os.getpid())    
    
    

def worker(halodata, out_q, info, inds, **kwargs):
    worker_q = Queue()
    if type(inds) == int:
        inds = [inds]
    for i in inds:
        mk_gal(h.data[i], worker_q, s.info, i, **kwargs)
    for i in inds:
        out_q.put(worker_q.get())


def halo_from_tree(tree_element, info):
    dtype_halo = [('id', '<i4'), ('m', '<f4'), ('mvir', '<f4'),
              ('r', '<f4'), ('rvir', '<f4'), 
              ('x', '<f4'), ('y', '<f4'), ('z', '<f4'),
              ('vx', '<f4'), ('vy', '<f4'), ('vz', '<f4')]
    
    cboxsize = 200.
    
    nout = tree_element['nout']
    h = hmo.Halo(nout=nout, halofinder='HM', info=info, is_gal=True)
    
    h.data = np.recarray(1, dtype=dtype_halo)
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
    h.aexp = tree_element['aexp']
    
    return h
#%%

"""
The processing pool needs to be instantiated in the main 
thread of execution. 
"""
import multiprocessing as mp    
import load
from tree import treemodule
import tree.ctutils as ctu
import tree.halomodule as hmo 
import os
import pickle
multi = True # 
hydro = True
is_gal = True
dump_gal = True

# optional parameters ----------------------------------------------------
lambda_plot = False 


#wdir = input("Working directory \n")
wdir = './'
#ncore=1
#nout_ini=186
#nout_end=186

if multi: ncore = int(input("How many cores? \n"))
nout_ini = input("First snapshot: (default = 37 , z=3) \n")
nout_fi = input("Last snapshot: (default = 187, z=0) \n")
out_dir = input("output directory: (default = out_001++)")
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

if nout_fi == "":
    nout_fi = 187
else:
    nout_fi = int(nout_fi)

nouts = range(nout_fi, nout_ini -1, -1) 

#zz = 1/nn.aout[1:] - 1 # 0-th aexp = 0. causes divide by 0 error.
#inds = mtc.find_closest(zz[::-1], [0,0.2,0.4,0.6,0.8,1.0, 1.5, 2.0, 3.0, 5.0])
#nouts = 187 - inds, where nn is the namelist.

# close to zred = [0,0.2,0.4,0.6,0.8,1.0, 1.5, 2.0, 3.0, 5.0]
#nouts = [187, 154, 130, 112,  98,  87,  67,  54,  37,  20]
nouts_dump_gal = [187, 154, 130, 112,  98,  87,  67,  54,  37,  20]
try: nouts_dump_gal
except NameError: nouts_dump_gal = None


#----------------------------------------------------------------------
m_min = 5e9
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
dir_out = out_base + 'catalog_GM/'
if not os.path.isdir(dir_out):
        os.mkdir(dir_out)


# In[71]:
from queue import Queue
import pandas as pd
import pickle
#import time
#os.system("taskset -p 0xff %d" % os.getpid())
# http://stackoverflow.com/questions/15414027/multiprocessing-pool-makes-numpy-matrix-multiplication-slower

"""
    at nout = 187, zoom-in region = 3 * Rvir_cluster.
    But at earlier nouts, there is no relation between clsuter radius and zoom-in region.
    Thus I need to load refine_params.txt
    Although r_refine determines the zoom-in region exactly, 
    let's give some buffer region. 
    For exmaple at nout = 187, Rvir_scale = 3.0, but I use only r < 2.5*Rvir_cluster
    which is ~ 83%.
    So I consider galaxies only within 85% of r_refine at each nout.

"""

rr = load.runparam.RefineParam()
rr.loadRegion(wdir + 'refine_params.txt')

nn = load.runparam.Nml()
nn.loadNml(wdir + 'cosmo_200.nml')

from utils import match

masscut_a = 1256366362.16
masscut_b = -20583566.5218

for inout, nout in enumerate(nouts):
    print(inout, nout, nout_fi)


    if nout in nouts_dump_gal:
        dump_gal = True
    else:
        dump_gal = False

    snout = str(nout)
    fcat = dir_out +"catalog" + snout + ".pickle"    
    galaxy_plot_dir = out_base + 'galaxy_plot' + snout + '/'
    if not os.path.isdir(galaxy_plot_dir):
        os.mkdir(galaxy_plot_dir)

    info = load.info.Info(nout=nout, base=wdir, load=True)
    
    m_min = masscut_a * info.aexp + masscut_b
    print("Mstar min:", m_min)
    # Load all halo
    hh = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info, load=True, is_gal=is_gal)
    i_center = np.where(hh.data['np'] == max(hh.data['np']))[0]
    
    i_aexp = match.closest(nn.aout[nout-1], rr.aexp) # or nout - 1?

    x_refine = rr.x_refine[i_aexp]
    y_refine = rr.y_refine[i_aexp]
    z_refine = rr.z_refine[i_aexp]
    r_refine = rr.r_refine[i_aexp] * 0.5
    
    i_satellites = smp.distance_to([x_refine, y_refine,z_refine],
                                   [hh.data['x'], hh.data['y'], hh.data['z']]) < 0.85 * r_refine
    
#    i_satellites = extract_halos_within(hh.data, i_center, info, dist_in_mpc = r_cluster * r_cluster_scale)
    print("Total {0} galaxies \n{1} galaxies are selected".format(
          len(i_satellites),sum(i_satellites)))
    
    large_enugh = hh.data['mvir'] > m_min
    hind = i_satellites * large_enugh
    h = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info, is_gal=is_gal)
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
                mstar_min=m_min, dump_gal = dump_gal)

    if multi == 1:
#   Multiprocessing -----------------------------------------------------------
        m = mp.Manager()
        out_q = m.Queue()
        print("Analyzing galaxies inside {} halos".format(nh))
        inds=[]
        [inds.append([]) for i in range(ncore)]
        
        for i in range(nh):
            j = i % ncore
            inds[j].append(i)

#        print(inds)
        processes = [mp.Process(target=worker, args=(h.data, out_q,
                    s.info, inds[i]), kwargs=keywords) for i in range(ncore)]

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
                              s.info, i), kwds=keywords)
        pool.close()
        pool.join()
    else:
        for i in range(nh):
            out_q = Queue()
            mk_gal(h.data[i], out_q, s.info, i, **keywords)
    
    print("----------Done---------")

#%%    
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
            dd =  out_q.get(timeout=0.1)
            if dd['id'] == 0:
                continue
            f.write("{:<4}   {:<4}   {:.5f}  {:.5f}  {:.5f}".format(nout,
                    dd['id'],dd['xc'],dd['yc'],dd['zc']))
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

