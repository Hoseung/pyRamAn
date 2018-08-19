# coding: utf-8


"""
Created on Mon Jun  1 00:15:25 2015

@author: hoseung
"""
from galaxymodule import galaxy
import load
import tree
import numpy as np
import utils.sampling as smp
from utils import util

#ncore = int(input("How man cores? \n"))
#wdir = input("Working directory \n")
wdir = '/home/hoseung/Work/data/AGN2/'
ncore = 2
nout=132
snout = str(nout)
rscale = 0.8
npix=800

info = load.info.Info(nout=nout, base=wdir)

frefine= 'refine_params.txt'
fnml = 'cosmo_200.nml'

ptypes=["star id pos mass vel", "dm id pos mass vel"]

# Load all halo
hall = tree.halomodule.Halo(nout=nout, base=wdir, halofinder="RS", info=info)
#hall = tree.halomodule.Halo(nout=nout, base=wdir, halofinder="HM", info=info)

hall.load()
# convert to code unit. - done by default
#hall.normalize()

# subset of halos ONLY inside zoom-in region
i_center = np.where(hall.data['np'] == max(hall.data['np']))[0]
h_ind = smp.extract_halos_within(hall, i_center, scale=2.0)
#%%
h = tree.halomodule.Halo()
h.derive_from(hall, h_ind)
#h.derive_from(hall, h_ind)#,  [4921, 5281, 5343, 5365, 5375, 5412, 5415], 5442, 5639, 5665, 6095])

region = smp.set_region_multi(h.data.x, yc=h.data.y, zc=h.data.z, radius = h.data.rvir * rscale)
print(region)

hind = np.where((hall.data.x > region["xr"][0]) & (hall.data.x < region["xr"][1]) &
                (hall.data.y > region["yr"][0]) & (hall.data.y < region["yr"][1]) &
                (hall.data.z > region["zr"][0]) & (hall.data.z < region["zr"][1]) &
                (hall.data.mvir > 1e11))[0]
h.derive_from(hall, hind[5:7])
region = smp.set_region_multi(h.data.x, yc=h.data.y, zc=h.data.z, radius = h.data.rvir * rscale)

#%%
s = load.sim.Sim()
s.setup(nout, wdir)

s.set_ranges(region["ranges"])
s.show_cpus()

#%%
s.add_part(ptypes)
s.part.load()
# convert to km/s
s.part.star['vx'] *= s.info.kms
s.part.star['vy'] *= s.info.kms
s.part.star['vz'] *= s.info.kms
s.part.star["m"]  *= s.info.msun

if 'dm' in s.part.pt:
    s.part.dm['vx'] *= s.info.kms
    s.part.dm['vy'] *= s.info.kms
    s.part.dm['vz'] *= s.info.kms
    s.part.dm["m"]  *= s.info.msun

#%%
s.add_hydro()
s.hydro.amr2cell(lmax=19)

#%%
print("Go!")
# Now, make galaxy data set
import numpy as np
dtype_catalog = [('id', int), ('mtot', '<f8'), ('mgas', '<f8'),
                 ('mstar', '<f8'), ('mdm', '<f8'), ('mhal', '<f8'),
                 ('nstar', int), ('ndm', int), ('nsink', int),
                 ('pos', '<f8', 3), ('vel', '<f8', 3), ('lx', '<f8', 3),
                 ('dcluster', '<f8'), ('b2t', float), ('mag', float, 5),
                 ('sfr', '<f8'), ('lambdar', '<f8'), ('lambdar_arr', '<f8', 30), 
                 ("nfit",float), ("morph_vi",str),
                 ("morph_b2t",str), ("add2",float), ("add3",float)]


x_clu = hall.data['x'][i_center]
y_clu = hall.data['y'][i_center]
z_clu = hall.data['z'][i_center]

util.reimport(galaxy)
#util.reimport(draw)
#util.reimport(draw.pp)
#util.reimport(draw.img_obj)

#plt.ioff()
import multiprocessing as mp
import queue
import matplotlib.pyplot as plt
# First halo is the cen

def mkgal(s, halodata, info, i_queue, out_q, star=None, dm=None, cell=None, rscale=0.4):
    while True:
        try:
            i = tasks.get(block=False)

            print(i,'/', halodata['id'][i])    
            gal = galaxy.Galaxy(halodata[i], radius_method='simple', info=info)
            gal.mk_gal(star=s.part.star, dm=None, cell=None,
                       rscale=rscale)
            gal_out = {"id":0, "xc":0, "yc":0, "zc":0, "mstar":0.0, "nstar":0,
                       "lambda_arr":[], "lambda_r":0}
                       
            if gal.star is False:
                print(gal.id, " Not a good galaxy")
                out_q.put(gal_out)
                continue
            else:
                print("R eff:", gal.reff * info.pboxsize)
                gal.cal_lambda_r(npix=20, method=1, rscale=1.5) # calculate within 1.0 * reff    
    #            gal.plot_gal(base=wdir + 'galaxy_plot4/')
    #            gal.save_gal(base=wdir)
            
                # Instead of galaxy class, save them in a dict.    
                gal_out['mstar'] = gal.mstar
                gal_out['nstar'] = gal.nstar
                gal_out['id'] = gal.id
                gal_out['xc'] = gal.xc
                gal_out['yc'] = gal.yc
                gal_out['zc'] = gal.zc
                gal_out['lambda_arr'] = gal.lambda_arr
                gal_out['lambda_r'] = gal.lambda_r
                out_q.put(gal_out)

        except queue.Empty:
            print("Queue closed. Exiting.")
            return
        except:
            continue

nh = len(h.data)
out_q = mp.Queue()
tasks = mp.Queue()

for i in range(nh):
    tasks.put(i) # send tasks to workers

processes = [mp.Process(target=mkgal, args=(s, h.data, s.info, tasks, out_q)) 
             for i in range(ncore)]

# run processes
for p in processes:
    p.start()

# exit completed processes   
print("QUEUE Done")
#print(out_q.empty())
for p in processes:
    p.join()
    print(p.name,"has terminated and joined.")

print("-------------------- Done --------------------")
#%%
for i in range(nh):
    print("dict out", i)
    dictout = out_q.get(timeout=2)

#%%
import pandas as pd
catalog = pd.DataFrame(dictout).to_records()

import pickle
with open(wdir + "catalog.pickle", 'wb') as f:
    pickle.dump(catalog, f)
#%%
#with open(wdir + 'catalog.pickle', 'rb') as f:
#    catalog = pickle.load(f)a
    
#dir+'000014gal.hdf5', "r")
#sx = infile['galaxy/star/x']
#infile.close()
#util.reimport(utils)

#%%    
import utils

# only massive galaxies
lambdar_arr = np.asarray(catalog["lambda_arr"])
i_cat = np.where(catalog['id'] != 0)
catalog = catalog[i_cat]
#%%

i_truegal = catalog['mstar'] > 3e9

disk_list=[]

# Exclude seemingly interacting galaxies
id_nogood = [ ]#, 6033, 6179]
#id_interacting=[5886,8909,6158,6226,]
#id_nogood += disk_list
i_ng = utils.match.match_list_ind(catalog['id'], id_nogood)
tflist=np.full(len(catalog), True, dtype=bool)
tflist[i_ng] = False

#i_voffset = np.wehre(lambdar_arr)

# intersection of two criteria
i_ok = np.logical_and(i_truegal, tflist)
#%%
f = plt.figure()
ax = f.add_subplot(111)
#for i, val in enumerate(lambdar_arr):
cnt = 0
for i, ok in enumerate(i_ok):
    if ok:
        cnt += 1
        if catalog[i]['id'] in disk_list :
            print("disk")
            ax.plot(lambdar_arr[i], 'r-', alpha=0.5) # up to 1Reff
            pass
#                ax.plot(val, 'r-') # up to 1Reff
        else:
            ax.plot(lambdar_arr[i], 'b-', alpha=0.3) # up to 1Reff
            
print(cnt)
#plt.xlabel() # in the unit of Reff
ax.set_title(r"$\lambda _{R}$") 
ax.set_ylabel(r"$\lambda _{R}$") 
ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
ax.set_xlim(right=9)
ax.set_xticks([0, 4.5, 9])
ax.set_xticklabels(["0", "0.5", "1"])
plt.savefig(wdir + "lambda_disk.png")
plt.close()

#%%
ll = catalog['lambda_r'][i_ok]

"""
f = plt.figure()
ax = f.add_subplot(111)
ax.scatter(catalog['dcluster'][i_ok],ll)
#ax.scatter(catalog['dcluster'][i_disk],lld, color='r')
#           catalog['lambdar_arr'][i_ok])# catalog['lambdar'][i_ok])
ax.set_xlim(0,1.1)
ax.set_ylim(0,1)
ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
ax.set_ylabel(r"$\lambda_{R}$")
plt.savefig(wdir + "hdclustervslambdar.png")
plt.close()
"""

f = plt.figure()
ax2 = f.add_subplot(111)
ax2.scatter(np.log10(catalog['mstar'][i_ok]), ll)
#ax2.scatter(np.log10(catalog['mstar'][i_disk]),lld, color='r')
# catalog['lambdar'][i_ok])#catalog['lambdar'][i_ok] )
ax2.set_xlim([9,11])
ax2.set_ylim(0,1)
ax2.set_xlabel("Stellar mass " + r"$[10^{10} M_{\odot}]$")
ax2.set_ylabel(r"$\lambda_{R}$")
plt.savefig(wdir + "msvslambdar.png")
plt.close()

