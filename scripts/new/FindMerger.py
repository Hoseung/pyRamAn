
# coding: utf-8

# ###  find merger event
# 
# 
#To derive merger mass ratio, I need mass of galaxy as the sum of particle mass. 
#1) add mass2 to consistent tree parameter set.
#2) retrieve data from original GalaxyMaker output. 

#I prefer 2). 
# In[149]:

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
        Add more values to existing tree data.
        
        1) construct a index list so that tree[i] = gal_org[ind[i]] for a given nout.
        2) 
    """
    

    dtype_new_quantities = [('np', '<i4'), ('id', '<i4'), ('m', '<f4'), ('mvir', '<f4'),
                            ('r', '<f4'), ('rvir', '<f4'), ('tvir', '<f4'), ('cvel', '<f4'),
                            ('x', '<f4'), ('y', '<f4'), ('z', '<f4'),
                            ('vx', '<f4'), ('vy', '<f4'), ('vz', '<f4'),
                            ('ax', '<f4'), ('ay', '<f4'), ('az', '<f4'),
                            ('sp', '<f4')]
    if is_gal:
        dtype_new_quantities.append([('sig', '<f4'), ('sigbulge', '<f4'), ('mbulge', '<f4')])
    
    
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
        keep_fields.append(["sig", "sigbulge", "mbulge"])
        
    return join_struct_arrays([treedata, New_arr[keep_fields]])

##########################################################


# In[166]:

from tree import treemodule
from tree import treeutils
import numpy as np

is_gal = False

alltrees = treemodule.CTree()
wdir = '/home/hoseung/Work/data/05427/'

if is_gal:
    # Galaxy tree
    tree_path = 'rhalo/Trees/tree_0_0_0.dat'
else:
    # halo tree
    tree_path = 'halo/Trees/tree_0_0_0.dat'

alltrees.load(filename= wdir + tree_path)
    
# Fix nout
nout_max = alltrees.data['nout'].max()
nout_fi = 187
alltrees.data['nout'] += nout_fi - nout_max

alltrees.data = augment_tree(alltrees.data, wdir, is_gal=is_gal)


# In[262]:

def get_npr(treedata, idx):
    ind = np.where(treedata['id'] == idx)[0]
    return treedata['nprog'][ind][0] # I want int, not a int array

def idx_to_ind(treedata, idx):
    return np.where(treedata['id'] == idx)[0]

def extract_a_tree(alltrees, idx_last):
    """
        Returns one full tree.
    """
    return alltrees[np.where(alltrees['tree_root_id'] == idx_last)]

#def extract_main_tree(alltrees, idx_last):
#    return alltrees[np.where((alltrees['tree_root_id'] == idx_last) & (alltrees['mmp'] == 1))]

def get_progenitors(treedata, idx, main=False):
    """
        Returns progenitors of a given halo/galaxy. 
        (from only one previous snapshot)
    """    
    if main:
        iprgs = np.where((treedata['desc_id'] == idx) & (treedata['mmp'] == 1))
    else:
        iprgs = np.where(treedata['desc_id'] == idx)

    return treedata['id'][iprgs]

def extract_main_tree(treedata, idx=None):
    """
        Returns a single branch/trunk of tree following only the main progenitors.
        Works whether the treedata is alltrees or atree.
        Search until no progenitor is found. Doesn't matter how long the given tree is. 
    """
    if idx == None:
        print("No idx is given")
        idx = treedata['id'][0]
        print("idx = ", idx)
    
    nprg = 1
    idx_list=[]
    ind_list=[]
    
    while nprg > 0:
        idx = get_progenitors(treedata, idx, main=True)[0]
        ind_list.append(np.where(treedata['id'] == idx)[0])
        nprg = get_npr(treedata, idx)
        
    return treedata[ind_list]

def main_thread(atree, idx):
    main_prg = extract_main_tree(atree, idx)
    
    ind = np.where(atree['id'] == idx)[0]
    nout_now = atree['nout'][ind]
    nout_fi = atree['nout'].max()
    if nout_now < nout_fi:
        desc = idx
        ind_desc_list=[]
        while desc > 0:
            ind = np.where(atree['id'] == desc)[0]
            ind_desc_list.insert(0,ind)
            desc = atree['desc_id'][ind]

    return np.concatenate((atree[ind_desc_list],main_prg))#,axis=1)

def last_halos(treedata, return_ind=False):
    nout_max = treedata['nout'].max()
    if return_ind:
        return np.where(treedata['nout'] == nout_max)[0]
    else:
        return treedata['id'][np.where(treedata['nout'] == nout_max)]


# In[169]:

i_last_halo = last_halos(alltrees.data, return_ind=True)
final_halos = alltrees.data[i_last_halo]
i_cluster = np.argmax(final_halos['np'])
cluster = final_halos['id'][i_cluster]
#print(final_halos['mvir'][i_cluster])


# In[302]:

atree = extract_a_tree(alltrees.data, cluster)
print(atree.dtype)
main = extract_main_tree(alltrees.data, cluster)
x_clu, y_clu = atree['x'][0], atree['y'][0]


# In[198]:

# Calculate merger event parameters
def find_merger(atree, idx=None, aexp_min=0.0):
    """
        find indices of merger event from a tree.
        (Full tree or main progenitor trunk)
    """
    if idx == None:
        idx = atree['id'][0]
        
    nprg = 1
    merger_list=[]

    i = 0
    while nprg > 0:
        idx = get_progenitors(atree, idx, main=True)[0]
        ind = np.where(atree['id'] == idx)[0]
        if atree['aexp'][ind] < aexp_min:
            break
        nprg = get_npr(atree, idx)

        if nprg > 1:
            merger_list.append(i)
        i +=1
    return merger_list


def merger_mass_ratio(atree, idx=None):
    """
    return mass ratio of the given merger event
    """
    if idx == None:
        idx = atree['id'][0]
        
    prgs = get_progenitors(atree, idx)
    
    # only for mergers
    if len(prgs) > 1:
        i_prgs = [np.where(atree['id'] == i)[0] for i in prgs]
        mass = []
        for iprg in i_prgs:
            mass.append(atree['m'])
    else:
        print("This is not a merger")
        return 0
    

def merger_properties_main_prg(atree, idx):
    """
        Calculate merger mass ratio for "one" merger event.

    if idx == None:
        if nout == None:
            print("Both idx and nout are missing")
            return
    else:
        if nout == None:
            nout = np.where(atree['id'] == idx)[0]

    idx = atree['id'][ind]
    """    

    #prgs = get_progenitors(atree, idx)
    #if len(prgs) > 1:
    #    i_prgs = [np.where(atree['id'] == i)[0] for i in prgs]
    
    i_prgs = np.where(atree['desc_id'] == idx)[0]
        
    print(i_prgs)
    id_prgs = atree['id'][i_prgs]
    mass_prgs = atree['m'][i_prgs]
    
    #mass_prgs_norm = mass_prgs / sum(mass_prgs)

    return mass_prgs




# In[303]:

# Create halo evolution animation

import draw.pp as pp
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')

# color some halos
i_z1 = np.where(atree['aexp'] == 0.5)[0]
top5 = np.argsort(atree['mvir'][i_z1])[-5:]
id_top5 = atree['id'][i_z1][top5]

large_branches=[]
for idx in id_top5:
    large_branches.append(main_thread(atree, idx))

x2 = np.zeros((100, len(id_top5)))
y2 = np.zeros((100, len(id_top5)))
r2 = np.zeros((100, len(id_top5)))
for i, large_branch in enumerate(large_branches):
    for j in range(100):
        x2[j,i] = large_branch['x'][j]
        y2[j,i] = large_branch['y'][j]
        r2[j,i] = large_branch['rvir'][j]

for nout in range(11,188):
    plt.cla()
    ax.set_xlim([x_clus-3, x_clus+3])
    ax.set_ylim([y_clus-3, y_clus+3])
    
    halos_this_nout = atree[np.where(atree['nout'] == nout)[0]]

    ax.set_title("aexp = {:<6.3f}".format(halos_this_nout['aexp'][0]))
    x = halos_this_nout['x']
    y = halos_this_nout['y']
    r = halos_this_nout['rvir'] / 1000
    for i in range(len(halos_this_nout)):
        pp.circle_scatter(ax, x[i], y[i], r[i], facecolor='none', edgecolor='b')
    if nout > 90:
        for i in range(len(id_top5)):
            pp.circle_scatter(ax, x2[187-nout,i], y2[187-nout,i], r2[187-nout,i]/1000, facecolor='none', edgecolor='r')
    plt.savefig(wdir + str(nout) + '.png')
#plt.show()


# ##### Merger mass ratio
mergers = find_merger(atree, cluster, aexp_min=0.1)
print(mergers)
main_prg = extract_main_tree(atree)
prg_mass = merger_properties_main_prg(atree, main_prg['id'][87])
MM_ratio= 0.05
for i in range(120,176):
    i_prgs = np.where(atree['desc_id'] == main_prg['id'][i])[0]
    if len(i_prgs) > 1:
        id_prgs = atree['id'][i_prgs]
        mass_prgs = atree['m'][i_prgs]   
        mass_ratios = mass_prgs / max(mass_prgs)
        print(mass_prgs)
        major_mergers = i_prgs[np.where(mass_ratios > MM_ratio)[0]]
        
        if len(major_mergers) > 1:
            print(major_mergers)
            print("Major Merger at nout = {}".format(i), atree['m'][major_mergers])
