# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 00:15:25 2015

@author: hoseung
"""
from galaxymodule import galaxy


def get_center(x, y, z, m, method='default'):
    """
    Determines the center of mass of given particles.
    IDL version was to iteratively search for mass weighted mean position of
    particles within gradually decreasing search radius. 
    Let's do the same thing here.
    """
    import numpy as np
    r = 10. / 200 * 1000  # 10kpc/h
    # Stops if r < 10pc
    # or, if dx_3d < 5 pc
    msum = sum(m)
    xc_old = sum(x * m)/msum
    yc_old = sum(y * m)/msum
    zc_old = sum(z * m)/msum

    i = 0    
    while i < 10 and r > 1e-7 :
        
        msum = sum(m)
        xc = sum(x * m)/msum
        yc = sum(y * m)/msum
        zc = sum(z * m)/msum
        
        dx3 = np.sqrt((xc_old - xc)**2 + (yc_old - yc)**2 + (zc_old - zc)**2)
        if dx3 < 5e-8:
            break
        xc_old = xc
        yc_old = yc
        zc_old = zc

        # shrink the search radius by half.
        ind = np.where( (x-xc)**2 + (x-xc)**2 + (x-xc)**2 < r**2)
        if len(ind) < 100:
            break
        x = x[ind]
        y = y[ind]
        z = z[ind]

        r = 0.2 * max([x.ptp(), y.ptp(), z.ptp()])

        i +=1
    
    return xc, yc, zc


# set region to encircle multiple halos. 

def set_region(xc, yc, zc, radius):
    import utils.sampling as smp
    import numpy as np
    if len(xc) > 1:
        ranges = set_region_multi(xc=xc, yc=yc, zc=zc, radius=radius)
        rr = max([np.ptp(ranges[0:2]),
                  np.ptp(ranges[2:4]),
                  np.ptp(ranges[4:6])]) * 0.5
        
        region = smp.set_region(xc=0.5 * sum(ranges[0:2]),
                                yc=0.5 * sum(ranges[2:4]),
                                zc=0.5 * sum(ranges[4:6]), radius = rr)
    else:
        region = smp.set_region(xc=xc, yc=yc, zc=zc, radius = radius)
    
    return region


def set_region_multi(xc=[], yc=[], zc=[], radius = []):
    xmi=[]
    xma=[]
    ymi=[]
    yma=[]
    zmi=[]
    zma=[]
    for i in range(len(xc)):
        xmi.append(xc[i] - radius[i])
        xma.append(xc[i] + radius[i])
        ymi.append(yc[i] - radius[i])
        yma.append(yc[i] + radius[i])
        zmi.append(zc[i] - radius[i])
        zma.append(zc[i] + radius[i])
    
    return (min(xmi), max(xma), min(ymi), max(yma), min(zmi), max(zma))

import draw 
import matplotlib.pyplot as plt

def part_2_den(part, info, region=None, proj='z', npix=800, ptype=None):
    """
    creates img object, calculate 2dprojection map and return the img object.
    given part, info, region, it passes x,y,z,m,npix,info arrays to pp.den2d.
    """
    import numpy as np
    from draw import img_obj, pp
    if region is None:
        region={"ranges":[[0,1]]*3, "center":[0.5]*3, "radius":0.5, \
        "xr":[0, 1], "yr":[0, 1], "zr":[0, 1]}

    ind_ok = np.where((part.x > region["xr"][0]) & (part.x < region["xr"][1]) & \
                        (part.y > region["yr"][0]) & (part.y < region["yr"][1]) & \
                        (part.z > region["zr"][0]) & (part.z < region["zr"][1]))

    if len(ind_ok[0]) > 0:
        img = img_obj.MapImg(info=info, proj='z', npix=npix, ptype=ptype)
        img.set_region(region)
        img.set_data(pp.den2d(part["x"][ind_ok], part["y"][ind_ok], part["z"][ind_ok], \
                    part["m"][ind_ok], npix, info, cic=True, norm_integer=True, region=region))

        return img
    else:
        return False


#%%
import load
import tree
import numpy as np
import utils.sampling as smp
from utils import util

nout=132
snout = str(nout)
rscale = 0.8
npix=800
wdir = '/home/hoseung/Work/data/AGN2/'
info = load.info.Info(nout=nout, base=wdir)

frefine= 'refine_params.txt'
fnml = 'cosmo_200.nml'

ptypes=["star id pos mass vel time metal", "dm id pos mass vel"]

# Load all halo
#hall = tree.halomodule.Halo(nout=nout, base=wdir, halofinder="RS", info=info)
hall = tree.halomodule.Halo(nout=nout, base=wdir, halofinder="HM", info=info)

hall.load()
# convert to code unit. - done by default
#hall.normalize()

# subset of halos ONLY inside zoom-in region
i_center = np.where(hall.data['np'] == max(hall.data['np']))[0]
h_ind = smp.extract_halos_within(hall, i_center, scale=1.0)
#%%
h = tree.halomodule.Halo()
#h.derive_from(hall, h_ind[5:35])
h.derive_from(hall, h_ind)#,  [4921, 5281, 5343, 5365, 5375, 5412, 5415], 5442, 5639, 5665, 6095])

region = set_region(h.data.x, yc=h.data.y, zc=h.data.z, radius = h.data.rvir * rscale)

hind = np.where( (hall.data.x > region["xr"][0]) & (hall.data.x < region["xr"][1]) &
                (hall.data.y > region["yr"][0]) & (hall.data.y < region["yr"][1]) &
                (hall.data.z > region["zr"][0]) & (hall.data.z < region["zr"][1]) &
                (hall.data.mvir > 1e11))[0]
h.derive_from(hall, hind[10:20])

region = set_region(h.data.x, yc=h.data.y, zc=h.data.z, radius = h.data.rvir * rscale)
#%%
s = load.sim.Sim(nout, wdir)

s.set_ranges(region["ranges"])
s.show_cpus()
s.add_part(ptypes)
s.part.load()
# convert to km/s
s.part.star['vx'] *= s.info.kms
s.part.star['vy'] *= s.info.kms
s.part.star['vz'] *= s.info.kms
s.part.star["m"]  *= s.info.msun

from utils import cosmology
s.part.star['time'] = cosmology.time2gyr(s.part.star['time'], z_now=info.zred, info=info)

if 'dm' in s.part.pt:
    s.part.dm['vx'] *= s.info.kms
    s.part.dm['vy'] *= s.info.kms
    s.part.dm['vz'] *= s.info.kms
    s.part.dm["m"]  *= s.info.msun

s.add_hydro()
s.hydro.amr2cell(lmax=19)

# set range directly from a halo instance
#region = smp.set_region(xc=h.data.x, yc=h.data.y, zc=h.data.z, radius = h.data.rvir * rscale)
#util.reimport(draw)
#util.reimport(draw.pp)


#img = part_2_den(s.part.dm, s.info, region=region, npix=npix)
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#img.plot_2d_den(axes=ax1, show=False, vmax=1e10, dpi=400, cname='brg', zposition=False)
#draw.pp.pp_halo(hall, npix, ind=hind, region=region, rscale=30, \
#            axes=ax1, facecolor='none', edgecolor='b', name=True)
#plt.savefig(wdir + snout + "dmmap_w_halo.png")
#plt.close()

#%%
# Now, make galaxy data set
import numpy as np
dtype_catalog = [('id', int), ('mtot', '<f8'), ('mgas', '<f8'),
                 ('mstar', '<f8'), ('mdm', '<f8'), ('mhal', '<f8'), 
                 ('nstar', int), ('ndm', int), ('nsink', int),
                 ('pos', '<f8', 3), ('vel', '<f8', 3), ('lx', '<f8', 3),
                 ('dcluster', '<f8'), ('b2t', float), ('mag', float, 5),
                 ('sfr', '<f8'), ('lambdar', '<f8'), ('lambdar_arr', '<f8', 30), 
                 ("nfit",float), ("morph_vi",str),
                 ("morph_b2t",str), ("nvec",float), ("reff",float)]


x_clu = hall.data['x'][i_center]
y_clu = hall.data['y'][i_center]
z_clu = hall.data['z'][i_center]


#%%
util.reimport(galaxy)
util.reimport(draw)
util.reimport(draw.pp)
util.reimport(draw.img_obj)

plt.ioff()
nh = len(h.data.id)

import matplotlib.pyplot as plt
#from draw import pp

#fig, axs = plt.subplots(np.ceil(nh/5), 5)
#ax = axs[0,0]
# First halo is the cen

lambdar_arr=[]
catalog = np.zeros(nh, dtype=dtype_catalog)
for icat, i in enumerate(range(nh)):
    gal = galaxy.Galaxy(h.data[i], radius_method='simple', info=info)
    gal.mk_gal(star=s.part.star, dm=s.part.dm, cell=s.hydro.cell, rscale=0.4)
    print(gal.id)
    if gal.star is False:
        print(" Not a good galaxy")
        continue
#    gal.get_radius("eff")
    gal.cal_lambda_r(npix=20, method=1, rscale=1.0)
    gal.plot_gal(base=wdir + 'galaxy_plot4/')
    gal.reorient(gal.star)
    gal.cal_b2t()
    isort = np.argsort(gal.dist1d)
    print(" DONE \n")
#    gal.save_gal(base = wdir)
    lambdar_arr.append(gal.lambda_arr)
    catalog[i]['mstar'] = gal.mstar
    catalog[i]['reff'] = gal.reff
    catalog[i]['id'] = gal.id
    catalog[i]['pos'][0] = gal.xc
    catalog[i]['pos'][1] = gal.yc
    catalog[i]['pos'][2] = gal.zc
    catalog[i]['vel'][0] = gal.vxc
    catalog[i]['vel'][1] = gal.vyc
    catalog[i]['vel'][2] = gal.vzc
    catalog[i]['lambdar_arr'][0:len(gal.lambda_arr)] = gal.lambda_arr
    catalog[i]['lambdar'] = gal.lambda_r
    catalog[i]['nvec'] = gal.nvec
    catalog[i]['dcluster'] = np.sqrt((gal.xc - x_clu)**2 + 
                                    (gal.yc - y_clu)**2 + 
                                    (gal.zc - z_clu)**2)/hall.data['rvir'][i_center]

#    field = pp.den2d(gal.star["x"], gal.star["y"], gal.star["z"],
#                     gal.star["m"], 600, s.info, region=None,

#                    ngp=False, cic=True, tsc=False)
print(catalog)
#%%
import pickle
with open("catalog.pickle", 'wb') as f:
    pickle.dump(catalog, f)

import pyfits
f = "/home/hoseung/Work/data/AGN2/catalog132.fits"
hdu = pyfits.open(f)
cat = hdu[1].data
#%%
for i, idgal in enumerate(catalog['id'][1:]):
    print(i,idgal)
    ind = np.where(cat['hnu'] == idgal)
    try:
        catalog["b2t"][i] = 1 - cat['d2t'][ind]
    except:
        print(idgal, "is missing")
        
    try:
        catalog["nfit"][i] = cat['nfit'][ind]
    except:
        print(idgal, "is missing")
#    morp.append()
#%%
def save_catalog(catalog, base='./', fname='catalog.hdf5'):
    import h5py as hdf
    # Save data into a hdf5 file
    filename = base + fname
    
    # Create and open the file with the given name
    with hdf.File(filename, 'w', libver='latest') as outfile:
    
    # Make a data set called catalog
#        outfile.create_dataset("galaxy")
    
    # Create a group for stroing data
        cat = outfile.create_group(snout)
        # Store data under /galaxy with direct assignment        
        cat.create_dataset("galaxy/id", data=catalog['id'])
        cat.create_dataset("galaxy/mstar", data=catalog['mstar'])        
        cat.create_dataset("galaxy/reff", data=catalog['reff'])        

        cat.create_dataset("galaxy/pos", data=catalog['pos'])
#        cat.create_dataset("galaxy/y", data=catalog['y'])
#        cat.create_dataset("galaxy/z", data=catalog['z'])
        cat.create_dataset("galaxy/mhal", data=catalog['mhal'])
        cat.create_dataset("galaxy/vel", data=catalog['vel'])
#        cat.create_dataset("galaxy/vy", data=catalog['vyc'])
#        cat.create_dataset("galaxy/vz", data=catalog['vz'])

#        cat.create_dataset("galaxy/time", data=catalog.galaxy['time'])
        cat.create_dataset("galaxy/mag_r", data=catalog['mag'])
        cat.create_dataset("galaxy/b2t", data=catalog['b2t'])
        cat.create_dataset("galaxy/nfit", data=catalog['nfit'])
        cat.create_dataset("galaxy/nvec", data=catalog['nvec'])

save_catalog(catalog, base=wdir)        
#%%    

        
        
        
#%%
import utils

# only massive galaxies
ll = np.asarray(lambdar_arr)
i_cat = np.where(catalog['id'] != 0)
catalog = catalog[i_cat]
#%%

i_truegal = catalog['mstar'] > 3e9

# Exclude seemingly interacting galaxies
disk_list = [4332, 4711, 4824, 4925, 5131, 5331, 5385, 5609, 5718, 6033, 6147, 6179, 6202]
i_disk = utils.match.match_list_ind(catalog['id'], disk_list)

id_nogood = [4922, 5122, 5203, 5208, 5305, 5388, 5631, 
             5776, 5886, 5909, 5964, 6158, 6226, 1106,
             4382, 4926, 5759, 5845, 5860 ]#, 6033, 6179]
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
plt.savefig(wdir + "lambdar_disk.png")
plt.close()

#%%
ll = 0.5 * (catalog['lambdar_arr'][i_ok,4] + catalog['lambdar_arr'][i_ok,9])
lld = 0.5 * (catalog['lambdar_arr'][i_disk,4] + catalog['lambdar_arr'][i_disk,9])

f = plt.figure()
ax = f.add_subplot(111)
ax.scatter(catalog['dcluster'][i_ok],ll)
ax.scatter(catalog['dcluster'][i_disk],lld, color='r')
#           catalog['lambdar_arr'][i_ok])# catalog['lambdar'][i_ok])
ax.set_xlim(0,1.1)
ax.set_ylim(0,1)
ax.set_xlabel("["+ r'$R/R_{eff}$'+"]")
ax.set_ylabel(r"$\lambda_{R}$")
plt.savefig(wdir + "hdclustervslambdar.png")
plt.close()


f = plt.figure()
ax2 = f.add_subplot(111)
ax2.scatter(np.log10(catalog['mstar'][i_ok]), ll)
ax2.scatter(np.log10(catalog['mstar'][i_disk]),lld, color='r')
# catalog['lambdar'][i_ok])#catalog['lambdar'][i_ok] )
ax2.set_xlim([9,11])
ax2.set_ylim(0,1)
ax2.set_xlabel("Stellar mass " + r"$[10^{10} M_{\odot}]$")
ax2.set_ylabel(r"$\lambda_{R}$")
plt.savefig(wdir + "msvslambdar.png")
plt.close()
