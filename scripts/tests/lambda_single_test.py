# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 12:48:45 2015

@author: hoseung
""" 

def plot_region(star, region, haloid, pboxsize, npix=300, ids=[], color=None):
    import matplotlib.pyplot as plt
    import numpy as np
    import draw
    
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

    rgal = region['radius'] * pboxsize * 1000

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
    
    plt.savefig(galaxy_plot_dir+"2dmap_"+str(haloid).zfill(5)+'.png', dpi=144)
    plt.close()


def extract_data(halo, rscale=1.25, radius='rvir'):
    xc_tmp0 = halo['x']
    yc_tmp0 = halo['y']
    zc_tmp0 = halo['z']
       
    #rr_tmp0 = min([halo[radius] * rscale, 0.0002]) 
    rr_tmp0 = max([halo[radius] * rscale, 0.0001])
    # arbitrary! < 20kpc/h
    
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
        
    return star_all[ind_s], dm_all[ind_d], cell_all[ind_c]

def mk_gal(halodata, out_q, info, i, idx,
           save=False, rscale=1.5, verbose=False, galaxy_plot_dir='./',
           rscale_lambda=3.0, npix_lambda=50, npix=400, galaxy_plot=False,
           method_com=2, mstar_min=5e9):
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
               "lambda_arr":[], "lambda_r":0, "rgal":0, "idx":idx,
               "rhalo":halodata['rvir'], "boxtokpc":info.pboxsize*1000}
    
    rscale_extract = 1.1
    star, dm, cell = extract_data(h.data[i], rscale=rscale_extract, radius='rvir')
    if mstar_min > 0.:
        if sum(star['m']) * info.msun < mstar_min:
            print("(1)Not enough stars: {:.2f} Msun".format(sum(star['m']) * info.msun))
            print("Aborting... \n")
            print(" Not a good galaxy")
            out_q.put(gal_out)
            return
               
    # Direct plot ---------------------------------------------------------                                
    if galaxy_plot:
        import matplotlib.pyplot as plt
        plt.ioff()
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
        xticks = ["{:.2f}".format(x) for x in np.linspace(-rgal, rgal, num=5)]
        ax.set_xticklabels(xticks)
        ax.set_ylabel("position [kpc]")
        ax.set_yticks(np.linspace(0,npix,5))
        yticks = ["{:.2f}".format(y) for y in np.linspace(-rgal, rgal, num=5)]
        ax.set_yticklabels(yticks)
        
        plt.savefig(galaxy_plot_dir+"2dmap_"+str(halodata['id']).zfill(7)+'.png', dpi=144)
        plt.close()

    #Create galaxy ----------------------------------------------
    gal = galaxy.Galaxy(halodata, radius_method='eff', info=info)
    good_gal = gal.mk_gal_from_gal(star, dm, cell,
                        mstar_min=mstar_min,
               rscale=min([rscale,rscale_extract]), verbose=verbose, method_com=method_com)

    #-----------------------------------------------------------------------    
    if not good_gal:
        print(gal.id, " Not a good galaxy")
        out_q.put(gal_out)
    else:
        # Save to catalog -------------------------------------------------------
        gal.cal_lambda_r(npix=npix_lambda, method=1, rscale=rscale_lambda) # calculate within 1.0 * reff    
        # Calculate lambda_r ---------------------------------------------------
        gal.plot_gal(fn_save = galaxy_plot_dir + str(nout).zfill(3) \
                             + "_" + str(idx).zfill(7) + "_"  \
                             + str(gal.id) + ".png", ioff=True)

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
    return gal


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
import load
from tree import treemodule
import tree.ctutils as ctu
import numpy as np
import utils.sampling as smp
import tree.halomodule as hmo 
import os
import pickle

#wdir = input("Working directory \n")
#wdir = '/home/hoseung/Work/data/05427/'

wdir = './'
nout_fi = 187 # fixed. Tree depends on this value. 
#final_gal = int(input("idgal: \n"))
#nout_ini = input("First snapshot: (default = 37 , z=3) \n")
#nout_end = input("Last snapshot: (default = 187, z=0) \n")

final_gal = 40
nout_ini=37
nout_end=187


nout_ini=37
nout_fi=187
if final_gal == "":
    import sys
    sys.exit(0)

if nout_ini == "":
    nout_ini = 37
else:
    nout_ini = int(nout_ini)

if nout_end == "":
    nout_end = 187
else:
    nout_end = int(nout_end)

read_halo_list = False

#nout_end = 131
#----------------------------------------------------------------------
# 27 : z=4;  37 : z=3;  20 : ~ z=5
nouts = range(nout_fi, nout_ini -1, -1) 
#----------------------------------------------------------------------

mstar_min_plot = 1e9
mk_gal_rscale = 2.0 # unit of Rvir,halo
rscale = 1.5
# should show fly-bys..

npix=800
rscale_lambda = 3.0 # Reff unit.
npix_lambda = int(10 * rscale_lambda)
lmax = 19
ptypes=["star id pos mass vel time metal", "dm id pos mass vel"]

# optional parameters ----------------------------------------------------
lambda_plot = False 
do_galaxy_plot=True
## halo -----------------------------------------------------------
mstar_min = 5e8 # minimum halo mass above which galaxies are searched for. 
dir_out = wdir + 'catalog_GM/'

is_gal = True

#%%

if is_gal:
    # Galaxy tree
    tree_path = 'GalaxyMaker/Trees/'
else:
    # halo tree
    tree_path = 'halo/Trees/'
    
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

tt = alltrees.data
tt_final = tt[tt['nout'] == nout_fi]

info = load.info.Info(nout=nout_fi, base=wdir, load=True)

# list of all galaxt targets from nout_ini to nout_fi --------------------
# This is the sum of all galaxies of tree of each final galaxy.
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
    print(" It's new run starting from the final snapshot.")
    print(" Loading halo list and selecting only massive halo/galaxies"
            " close enough to the central cluster")           

    idx = ctu.id_to_idx(tt, final_gal, nout_fi)[0]
    all_gals_in_trees = ctu.extract_main_tree(tt, idx)

if not read_halo_list:    
    print(" Pickling halo list...")
    pickle.dump(all_gals_in_trees, open( wdir + "all_gals_in_trees.pickle", "wb" ))
    with open(wdir + 'lambda_mp_status.txt', 'w') as f:
        f.write("mstar_min = " + str(mstar_min) + "\n")
        f.write("ptypes : \n")
        for i in ptypes:
            f.write("  " + str(i) + "\n")
        
print("------ target galaxy list built")


#%%
#import galaxy
import utils.util
utils.util.reimport(load)
utils.util.reimport(load.utils)
try:
    if not os.path.isdir(dir_out):
        os.mkdir(dir_out)
    f = open(dir_out + 'galaxies' +str(final_gal) + '.txt', 'w')
except:
    print("No filename is given.\n ")

f.write(" Final ID:" + str(final_gal) + "\n")
f.write(" nout     ID        x          y       z[Mpc]       vx      vy     vz[km/s]")
f.write("    Reff[kpc]     Mstar    Mgas[Msun]  Rhalo[kpc]  boxtokpc  \n")    

print("!!!!!!!!!!!!!!!!!!")
import queue
out_q = queue.Queue()
dictout=[]
dd = None
galaxy_plot_dir = wdir + 'galaxy_' + str(final_gal).zfill(5) + '/'
if os.path.isdir(galaxy_plot_dir) is False:
    os.mkdir(galaxy_plot_dir)

#import matplotlib.pyplot as plt
#%%
#for inout, nout in enumerate(nouts):
inout, nout = 0, 187
#if nout > nout_end:
#    continue
print(inout, nout, nout_end)
snout = str(nout)

info = load.info.Info(nout=nout, base=wdir, load=True)
       
# Load all halo
hh = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info, load=True, is_gal=True)
h = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info, is_gal=True)
h.derive_from(hh, hh.data['id'] == all_gals_in_trees['Orig_halo_id'][all_gals_in_trees['nout'] == nout])
halodata = h.data
region = smp.set_region(xc=halodata['x'],
                    yc=halodata['y'],
                    zc=halodata['z'],
                    radius = halodata['rvir'])


s = load.sim.Sim(nout=nout, base=wdir, setup=True, region = region)
s.add_part(ptypes, load=True, fortran=True)

assert s.part.nstar > 0, "Not enough stellar particles in given cpus"
    
s.add_hydro(load=True, lmax=lmax)

star_all = s.part.star
dm_all = s.part.dm
cell_all = s.hydro.cell

keywords = dict(galaxy_plot_dir=galaxy_plot_dir,
            rscale = mk_gal_rscale,
            verbose=False, rscale_lambda=rscale_lambda,
            npix_lambda=npix_lambda, galaxy_plot = do_galaxy_plot,
            mstar_min=mstar_min)
# To keep the consistency, use queue to get output data. 

gal = mk_gal(h.data[0], out_q, info, 0, final_gal, **keywords)
#gal.save_gal('./test.h5')
  
tmp =  out_q.get()
#if tmp['id'] == 0:
#    continue
dd = tmp
f.write("{:<4}   {:<4}   {:.5f}  {:.5f}  {:.5f}".format(nout,
        dd['id'],dd['xc'],dd['yc'],dd['zc']))
f.write("  {:.3f}  {:.3f}  {:.3f}".format(dd['vx'],dd['vy'],dd['vz']))
f.write("  {:.6f}  {:.0f}  {:.0f}".format(dd['rgal'],dd['mstar'], dd['mgas']))
f.write("  {:.5f}  {:.5f}    \n".format(dd['rhalo'],dd['boxtokpc']))
dictout.append(dd)
print("  ")

print("----------Done---------")
f.close()    
print("Text file written")

import pandas as pd
catalog = pd.DataFrame(dictout).to_records()
# pickling also works!
# The problem was memory being locked!
import pickle
fcat = dir_out +"catalog" + str(nout_fi) + "_" + str(final_gal) + ".pickle"
with open(fcat, 'wb') as f:
    pickle.dump(catalog, f)  

#%%
import matplotlib.pyplot as plt 
import load
info = load.info.Info(nout=187, base='./')

#%%

if lambda_plot:
    i_early = np.where(catalog['mstar'] > mstar_min_plot)[0]
    i_late = []
    i_bad = np.where(catalog.id == 0)[0]    
    plot_lambda(catalog, i_early, i_late, i_bad,
                fn_out = snout + "_lambdar_disk.png")






#%%
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
            if not isinstance(attr, (np.ndarray, np.recarray, dict, list))}

def save_gal(gal, filename):    
    import h5py as hdf
    
    outfile = hdf.File(filename, 'w')
    
    # Store metadata using HDF5 attributes
    attrs = get_metadata(gal)
    attrs = get_metadata2(attrs)
    #outfile.attrs.create("all", attrs)
    print(attrs)
    for name, atr in attrs.items():
        outfile.attrs[name] = atr
    
    # Store data under /galaxy with direct assignment
    if hasattr(gal, 'star'):
        print("Saving star")
        star = outfile.create_group("star")
        for field in gal.star.dtype.names:
            star.create_dataset(field, data=gal.star[field])
    
    if hasattr(gal, 'dm'):
        print("Saving DM")
        dm = outfile.create_group("dm")
        for field in gal.dm.dtype.names:
            dm.create_dataset(field, data=gal.dm[field])
        
    if hasattr(gal, 'cell'):
        print("Saving gas")
        gas = outfile.create_group("gas")
        for field in gal.cell.dtype.names:
            gas.create_dataset(field, data=gal.cell[field])
    
    outfile.close()


#new_gal = hdf.File(filename, 'r')
#%%
save_gal(gal, './test.h5')
        
#%%
import h5py as hdf
infile = hdf.File('./galaxy_plot187/187_40.h5', 'r')

#gal = galaxy.galaxy.Galaxy()
#%%
# Store metadata using HDF5 attributes



#%%
def load_gal(filename):    
    import h5py as hdf
    import galaxy
    import numpy as np
    
    def get_dtype(file, grp):
        dtype =[]
        for field in file[grp]:
            dtype.append( (field, file[grp][field].dtype))
        
        return dtype
    

    def assign_data(target, grp):
        for field in grp:
            target[field] = grp[field]
            
    
    infile = hdf.File(filename, 'r')
    
    gal = galaxy.galaxy.Galaxy()
    for name, val in infile.attrs.items():
        setattr(gal, name, val)

    # initialize
    gal.star = np.recarray(infile['star']['x'].shape[0], dtype=get_dtype(infile, 'star'))
    gal.dm = np.recarray(infile['dm']['x'].shape[0], dtype=get_dtype(infile, 'dm'))
    gal.gas = np.recarray(infile['gas']['x'].shape[0], dtype=get_dtype(infile, 'gas'))
   
    # assign
    assign_data(gal.star, infile['star'])
    assign_data(gal.dm, infile['dm'])
    assign_data(gal.gas, infile['gas'])
    infile.close()
    
    return gal
    
#%%
gal = load_gal('./galaxy_plot187/187_40.h5')
#%%


# Now calculate bulge mass (and compare with GalaxyMaker result)












