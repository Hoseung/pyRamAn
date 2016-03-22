
# coding: utf-8

# In[2]:

import numpy as np
import pickle
from tree import treemodule

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


is_gal=True
wdir = './'

# Load complete tree -----------------------------------------------------
alltrees = treemodule.CTree()
if is_gal:
    # Galaxy tree
    tree_path = 'GalaxyMaker/Trees/'
else:
    # halo tree
    tree_path = 'halo/Trees/'


alltrees.load(filename= wdir + tree_path + 'tree_0_0_0.dat')

# Fix nout -----------------------------------------------------
nout_max = alltrees.data['nout'].max()
nout_fi = 187
alltrees.data['nout'] += nout_fi - nout_max
alltrees.data = augment_tree(alltrees.data, wdir, is_gal=is_gal)


pickle.dump(alltrees, open(wdir + tree_path + "extended_tree.pickle", "wb" ))

