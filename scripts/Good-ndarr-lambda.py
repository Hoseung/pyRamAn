
# coding: utf-8

# In[131]:

def mk_gal(halodata, out_q, info, i,
            save=False, rscale=0.3, verbose=False, galaxy_plot_dir='./'):
    from galaxymodule import galaxy
#    print("IDs:", id(star), id(dm), id(cell))
    rscale_lambda = 2
    print(halodata['id'])
    print("{}-th **********************************************".format(i))

    gal = galaxy.Galaxy(halodata, radius_method='simple', info=info)
    is_gal = gal.mk_gal(star=star, dm=dm, cell=cell,
               rscale=rscale, verbose=verbose)
    gal_out = {"id":0, "xc":0.0, "yc":0.0, "zc":0.0,
               "vx":0.0, "vy":0.0, "vz":0.0,
               "mstar":0.0, "nstar":0.0, "mgas":0.0,
               "lambda_arr":[], "lambda_r":0, "rgal":0}
    print("IS_GAL",is_gal)
    if not is_gal:
        print(gal.id, " Not a good galaxy")
        out_q.put(gal_out)
    else:
        print("Good galaxy, R eff:", gal.reff * info.pboxsize)
        gal.cal_lambda_r(npix=npix_lambda, method=1, rscale=rscale_lambda) # calculate within 1.0 * reff    
        gal.plot_gal(save_dir= galaxy_plot_dir, ioff=True)
        print("plot done")
    #            gal.save_gal(base=wdir)

        # Instead of galaxy class, save them in a dict.    
        gal_out['mstar'] = gal.mstar
        gal_out['mgas'] = gal.mgas
        gal_out['nstar'] = gal.nstar
        gal_out['id'] = gal.id
        gal_out['xc'] = gal.xc * info.pboxsize
        gal_out['yc'] = gal.yc * info.pboxsize
        gal_out['zc'] = gal.zc * info.pboxsize
        gal_out['vx'] = gal.vxc
        gal_out['vy'] = gal.vyc
        gal_out['vz'] = gal.vzc        
        gal_out['lambda_arr'] = gal.lambda_arr
        gal_out['lambda_r'] = gal.lambda_r
        gal_out['rgal'] = gal.reff * info.pboxsize * 1000.0 # in kpc
        print(gal_out)
        out_q.put(gal_out)
    return
    
def extract_halos_within(halos, ind_center, scale=1.0, Mcut=1e5):
    import numpy as np
    import utils.sampling as smp
    '''
    Returns indices of halos within SCALE * Rvir of the central halo.

    def extract_halos_within(halos, ind_center, scale=1.0)
    halos : halo finder output (single snapshot)
    ind_center : index of central halo
    scale : multiplying factor to the Rvir of the central halo
    '''
    xc = halos['x'][i_center]
    yc = halos['y'][i_center]
    zc = halos['z'][i_center]
    rvir= halos['rvir'][i_center]

    xx = halos['x']
    yy = halos['y']
    zz = halos['z']
    m = np.array(halos['m'])

    dd = smp.distance_to([xc,yc,zc],[xx,yy,zz])

    i_m = m > Mcut
    i_ok = np.logical_and(dd < (rvir * scale), i_m)

    return i_ok
#%%


# In[145]:

"""
The processing pool needs to be instantiated in the main 
thread of execution. 
"""
import multiprocessing as mp    
import load
from tree import tmtree
import numpy as np
import utils.sampling as smp
#import ctypes
import tree.halomodule as hmo 
from utils import match
#ncore = int(input("How man cores? \n"))
#wdir = input("Working directory \n")
#nout = int(input("nout? \n"))
wdir = './'
#wdir = '/home/hoseung/Work/data/AGN2/'
ncore = 16

rscale = 0.8
npix=800
rscale_lambda = 3.0
npix_lambda = 30
lmax = 19
ptypes=["star id pos mass vel time metal", "dm id pos mass vel"]
frefine= 'refine_params.txt'

## halo part ###
m_halo_min = 1e10 # minimum halo mass above which galaxies are searched for. 
fixed_position = True
nout_ini = 184
# 27 : z=4;  37 : z=3;  20 : ~ z=5
nout_fi = 187
rvir=3.0
# nout_halo = 122 == nout 10, nout_halo = 0   == nout 132
nouts = range(nout_fi, nout_ini -1, -1) 
Nnouts = len(nouts)

tt = tmtree.load(work_dir=wdir, filename="halo/TMtree.fits")
tfin = tt[np.where(tt['NOUT'] == 0)]
tini = tt[np.where(tt['NOUT'] == nout_fi - nout_ini)]
#%%
info = load.info.Info()
info.setup(nout=nout_fi, base=wdir)
info.load()
hh = hmo.Halo(base=wdir, nout=nout_fi, halofinder='HM', info=info)
hh.load()
halo = hh.data
#halo = hmu.load_data(nout_fi, work_dir=work_dir, normalize=True)
i_center = np.where(halo['np'] == max(halo['np']))
i_satellites = extract_halos_within(halo, i_center, scale=3.0)[0]
print("Total {0} halos \n{1} halos are selected".format(
      len(i_satellites),sum(i_satellites)))
#%%
# halos found inside the cluster and has tree back to nout_ini
halo_list = halo['id'][i_satellites]
#print(halo_list)
h_ind_ok, halo_ok = tmtree.check_tree_complete(tt, 0, nout_fi - nout_ini, halo_list)

print(len(halo_ok), "halos left")


# In[146]:

import time
import os
nouts = range(nout_fi, nout_fi - 1, -1) 
dir_out = wdir + 'catalog/'

fnml = 'cosmo_200.nml'
for inout, nout in enumerate(nouts):
    print(inout, nout)
    snout = str(nout)

    fcat = dir_out +"catalog" + snout + ".pickle"
    galaxy_plot_dir = wdir + 'galaxy_plot' + snout
    if os.path.isdir(galaxy_plot_dir) is False:
        os.mkdir(galaxy_plot_dir)
    
    info = load.info.Info()
    info.setup(nout=nout, base=wdir)
    info.load()
           
    # Load all halo
    hh = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info)
    hh.load()  

    hind = match.match_list_ind(hh.data['id'], halo_ok[:,inout])
    
    h = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info)
    h.derive_from(hh, hind)
    region = smp.set_region_multi(h.data.x, yc=h.data.y, zc=h.data.z, radius = h.data.rvir * rscale)

    #%%
    s = load.sim.Sim()
    s.setup(nout, wdir)    
    s.set_ranges(region["ranges"])
    s.show_cpus()
    s.add_part(ptypes)    
    s.add_hydro()
    t0 = time.time()
    s.part.load(fortran=True)
    t1 = time.time()
    s.hydro.amr2cell(lmax=19)
    t2 = time.time()
    print("Loading particle took {}, \n and loading hydro took {}".format(t1-t0, t2-t1))    
    star = s.part.star
    dm = s.part.dm
    cell = s.hydro.cell


# In[147]:

m = mp.Manager()
out_q = m.Queue()
ncore = 1
nh = len(h.data)
#    nh = 100
t3 = time.time()

import galaxy
import draw
utils.util.reimport(galaxy)
utils.util.reimport(galaxy.galaxy)
utils.util.reimport(draw.pp)
#    keywords = dict(galaxy_plot_dir=galaxy_plot_dir, verbose=False)
keywords = {"galaxy_plot_dir":'./', "verbose":False}
print(nh)
#    pool = mp.Pool(processes=ncore)
for i in range(nh):
    mk_gal(h.data[i], out_q, s.info, i, keywords)
#        pool.apply_async(mk_gal, args=(h.data[i], out_q, s.info, i, keywords))
#        print(i,"-th started")
    # Pass optional parameters in dict!s

#    pool.close()
#    pool.join()
print("----------Done---------")
print(" Took", time.time() - t3, "Seconds")

dictout=[]
try:
    f = open(dir_out + 'galaxies' + snout + '.txt', 'w')
except:
    print("No filename is given.\n ")

f.write(" #      ID        x          y       z[Mpc]       vx      vy     vz[km/s]")
f.write("    Reff""[kpc]      Mstar    Mgas[Msun]\n")    

for i in range(nh):
    print("dict out", i)
    dd = out_q.get(timeout=2)
    f.write("{:<4}   {:<4}   {:.5f}  {:.5f}  {:.5f}".format(i,
            dd['id'],dd['xc'],dd['yc'],dd['zc']))
    f.write("  {:.3f}  {:.3f}  {:.3f}".format(dd['vx'],dd['vy'],dd['vz']))
    f.write("  {:.6f}  {:.0f}  {:.0f}     \n".format(dd['rgal'],dd['mstar'], dd['mgas']))
    dictout.append(dd)

f.close()    
print("Text file written")
import gc
gc.collect()
import pandas as pd
catalog = pd.DataFrame(dictout).to_records()    
# pickling also works!
# The problem was memory being locked!
import pickle
with open(fcat, 'wb') as f:
    pickle.dump(catalog, f)  

i_early = np.where(catalog.id != 0)[0]
i_late = []
i_bad = np.where(catalog.id == 0)[0]
if not os.path.isdir(wdir + snout + '/'):
    os.mkdir(wdir + snout + '/')
plot_lambda(catalog, i_early, i_late, i_bad, out_dir = wdir + snout + '/')



# In[ ]:

def plot_lambda(catalog, i_early, i_late, i_bad, out_dir='./'):
    import matplotlib.pyplot as plt
    plt.ioff()
    f = plt.figure()
    ax = f.add_subplot(111)
    #for i, val in enumerate(lambdar_arr):
    cnt = 0
    for i in i_early:
        print(i)
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
    plt.savefig(out_dir + "lambdar_disk.png")
    plt.close()
    
#plot_lambda(catalog, i_early, i_late, i_bad, out_dir = wdir + snout + '/')


# In[139]:

a = catalog['lambda_arr'][0]


# In[141]:




# In[ ]:



