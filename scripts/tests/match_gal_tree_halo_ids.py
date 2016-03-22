# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 17:48:11 2016

@author: hoseung
"""

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
    # equivalent with allhal.data
    newhals = np.recarray(len(allgal.data), dtype=allhal.data.dtype) 
    
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


import load
import numpy as np
from tree import treemodule
import tree.ctutils as ctu
import tree.halomodule as hmo 
import os
import pickle
import utils.match as mtc
multi = 1 # 
hydro = True
is_gal = True
dump_gal = True

# optional parameters ----------------------------------------------------
lambda_plot = False 
galaxy_plot = True

wdir = './'
nout_ini=37
nout_end=187

nout_ini0 = 37 #intrinsic
nout_fi = 187 #intrinsic


read_halo_list = True


dir_increment = 0
out_dir = "out_" + str(dir_increment).zfill(3) + '/'
while os.path.isdir(out_dir):
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

#nouts = range(nout_fi, nout_ini -1, -1) 
nouts = [187, 121, 87, 54, 37]
nouts_dump_gal = [187, 154, 130, 121, 112,  98,  87,  67,  54,  37,  20]
try: nouts_dump_gal
except NameError: nouts_dump_gal = None
    

#----------------------------------------------------------------------
mstar_min = 5e9
# Only galaxies above this stellar mass at the final snapshot are considered.
mstar_min_plot = 5e9 # What is different from mstar_min?  
mk_gal_rscale = 1.1 # unit of Rvir,galaxy
r_cluster_scale = 2.5 # maximum radius inside which galaxies are searched for
npix=800
rscale_lambda = 3.0 # Reff unit.
npix_lambda = 10 # per 1Reff
lmax = 19
ptypes=["star id pos mass vel time metal", "dm id pos mass vel"]


## halo part -----------------------------------------------------------
dir_out = wdir + 'catalog_reori/'




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
    print("Cannot load an extended tree")
    alltrees = treemodule.CTree()
    alltrees.load(filename= wdir + tree_path + 'tree_0_0_0.dat')
    # Fix nout -----------------------------------------------------
    nout_max = alltrees.data['nout'].max()
    alltrees.data['nout'] += nout_fi - nout_max
    print("------ NOUT fixed")
    alltrees.data = ctu.augment_tree(alltrees.data, wdir, is_gal=is_gal)
    print("------ tree data extended")
    pickle.dump(alltrees, open(wdir + tree_path + "extended_tree.pickle", "wb" ))

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
             





#for inout, nout in enumerate(nouts):
inout = 100
nout = 87

print(inout, nout, nout_fi)

snout = str(nout)
dump_gal = False

masscut_a = 1256366362.16
masscut_b = -20583566.5218  


fcat = dir_out +"catalog" + snout + ".pickle"
galaxy_plot_dir = out_base + 'galaxy_plot' + snout + '/'
if not os.path.isdir(galaxy_plot_dir):
    os.mkdir(galaxy_plot_dir)

info = load.info.Info(nout=nout, base=wdir, load=True)

mstar_min = min([masscut_a * info.aexp + masscut_b, mstar_min])
print("Mstar min: {:.2e}".format(mstar_min))

t_now = alltrees.data[alltrees.data['nout'] == nout]
idxs_now = all_gals_in_trees[inout]

h = halo_from_tree(t_now[mtc.match_list_ind(t_now['id'], np.array(idxs_now))], info)


allhal = hmo.Halo(base=wdir, nout=nout, is_gal=False, halofinder='HM',
                  return_id=False, load=True)
clu = allhal.data[allhal.data.np.argmax()]

# Galaxy with enough stellar mass
allgal = hmo.Halo(base=wdir, nout=nout, is_gal=True, halofinder='HM',
                  return_id=False, load=True)
i_gal_ok = allgal.data['m'] > mstar_min

i_gal_near = find_halo_near(allgal.data, clu, rscale=2.5)
id_gal_ok = allgal.data['id'][i_gal_ok * i_gal_near]


# get the stellar particle IDs of OK galaxies.
allgal = hmo.Halo(base=wdir, nout=nout, is_gal=True, halofinder='HM',
                  return_id=True)
allgal.set_return_id_list(id_gal_ok)
allgal.load()
allgal.data = allgal.data[i_gal_ok * i_gal_near] # only good gals

allgal.data['idx'][:] = -1

# associate halos with galaxies. #########################################
# Virial radius is required to measure gas mass associated with the galaxy
allhal.data = allhal.data[find_halo_near(allhal.data, clu, rscale=2.5)] # only good hals

# new halo catalog with 1:1 correlation to gals.
allhal = associate_gal_hal(allgal, allhal) 

# get the DM particle IDs of OK galaxies.
halo_tmp = hmo.Halo(base=wdir, nout=nout, is_gal=False, halofinder='HM', return_id=True, load=False)
halo_tmp.set_return_id_list(allhal.data['id'])
halo_tmp.load()

# move tmp.idlist to allhal with re-ordering.
# Now idlist of allhal and allgal are corresponding to each other.
allhal.idlists=[]
for thisid in allhal.data['id']:
    if thisid < 0: # galaxy with host halo.
        allhal.idlists.append([])
    else:
        allhal.idlists.append(halo_tmp.idlists[np.where(halo_tmp.hal_idlists == thisid)[0]])

##########################################################################
# Only the galaxies in all_gals_in_tree
for gal in allgal.data:
    ii =  np.where(h.data['id'] == gal['id'])[0]
    if len(ii) > 0:
        gal['idx'] = h.data['idx'][ii]
        print(ii, "idxxxxxx", gal['id'], gal['idx'])

allhal.data['idx'] = allgal.data['idx']
print("idxxxxxx", allgal.data['idx'])