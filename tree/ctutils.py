# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 17:10:10 2015

@author: hoseung
"""
from tree import treemodule
#from tree import treeutils
import numpy as np
#%%

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
        gal_org = hmo.Halo(base=base, nout=nout, halofinder='HM', load=True, is_gal=is_gal)
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
            New_arr[ind_tree_this_nout[i]] = gal_org[ind].copy()
 
    # Drop duplicate fields
    #["id", "mvir", "rvir", "x", "y", "z", "vx", "vy", "vz"]
    keep_fields = ["np", "m", "r", "tvir", "cvel"]
    if is_gal:
        [keep_fields.append(i) for i in ['sig', 'sigbulge', 'mbulge']]
        
    return join_struct_arrays([treedata, New_arr[keep_fields]])


def load_tree(wdir, is_gal=False, no_dump=False, nout_fi=187):
    import pickle
    #import tree.ctutils as ctu

    alltrees = treemodule.CTree()   

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
        print("No extended tree is found")
        alltrees.load(filename= wdir + tree_path + 'tree_0_0_0.dat')
        if not no_dump:
            # Fix nout -----------------------------------------------------
            nout_max = alltrees.data['nout'].max()
            alltrees.data['nout'] += nout_fi - nout_max
            print("nout_fi = {}, fixing nout of the tree".format(nout_fi))
            alltrees.data = augment_tree(alltrees.data, wdir, is_gal=is_gal)
            print("------ tree data extended")
        
    return alltrees


def get_npr(treedata, idx):
    """
        Returns number of progenitors of a given idx
    """
    ind = np.where(treedata['id'] == idx)[0]
    return treedata['nprog'][ind][0] # I want int, not a int array


def get_progenitors(treedata, idx, main=False):
    """
        Returns progenitors of a given halo/galaxy. 
        (from only one previous snapshot)
    """    
    if main:
        iprgs = np.where((treedata['desc_id'] == idx) & (treedata['mmp'] == 1))
        return treedata['id'][iprgs]
    else:
        iprgs = np.where(treedata['desc_id'] == idx)
        return treedata['id'][iprgs]


def idx_to_ind(treedata, idx):
    return np.where(treedata['id'] == idx)[0]

def id_to_idx(treedata, id, nout):
    return treedata['id'][np.where((treedata['Orig_halo_id'] == id) & (treedata['nout'] == nout))]
    

def extract_a_tree(alltrees, idx_last):
    """
        Returns one full tree.
    """
    return alltrees[np.where(alltrees['tree_root_id'] == idx_last)]

def extract_satellite(atree, idx):
    """
        extract a partial tree of satellite.
        Starting from a certain nout, goes back in time. 
    """
    pass
    

def extract_main_tree(treedata, idx=None, no_subset=False):
    """
        Returns a single branch/trunk of tree following only the main progenitors.
        Works with both alltrees or atree.
        Search until no progenitor is found. Doesn't matter how long the given tree is. 
        Only earlier snapshots are searched for.
    """
    
    if idx == None:
        print("No idx is given")
        idx = treedata['id'][0]
        print("idx = ", idx)

    if no_subset:
        smalldata = treedata
    else:
        smalldata = treedata[treedata['tree_root_id'] == idx]

    nprg = 1
    ind_list=[np.where(smalldata['id'] == idx)[0][0]]
      
    while nprg > 0:
        idx = get_progenitors(smalldata, idx, main=True)
        ind_list.append(np.where(smalldata['id'] == idx[0])[0][0])
        nprg = get_npr(smalldata, idx[0])

    return smalldata[ind_list]


def extract_main_tree_full(atree, idx):
    """
        extract a FULL main tree including later snapshots.
        = extract_main_tree + later snapshots
    """
    main_prg = extract_main_tree(atree, idx)
    
    nout_now = atree['nout'][np.where(atree['id'] == idx)[0]]
    nout_fi = atree['nout'].max()
    if nout_now < nout_fi:
        desc = idx
        ind_desc_list=[]
        while desc > 0:
            ind = np.where(atree['id'] == desc)[0]
            ind_desc_list.insert(0,ind)
            desc = atree['desc_id'][ind]
        return np.concatenate((atree[ind_desc_list],main_prg))#,axis=1)
    elif nout_now == nout_fi:
        return main_prg


def last_halos(treedata, return_ind=False):
    nout_max = treedata['nout'].max()
    if return_ind:
        return np.where(treedata['nout'] == nout_max)[0]
    else:
        return treedata['id'][np.where(treedata['nout'] == nout_max)]


def idxs_to_ids(alltree, idxs):
    return [alltree["Orig_halo_id"][np.where(alltree['id'] == idx)][0] for idx in idxs]


def check_tree_complete(treedata, nout_ini, nout_fi, halo_list,\
                        idx=False):
    """
    returns a list of halo IDs at nout_fi that with trees fully linked
    from nout_ini to nout_fi.
    
    parameters
    ----------
    idx : if False, halo_list are ids from halo bricks.
    """    

    if idx:
        idx_list = halo_list
    else:
        halos_now = treedata[np.where(treedata['nout'] == nout_fi)]
        # Not all galaxies are found in the tree.
        idx_list = []
        for hal in halo_list:
            ii = np.where(halos_now['Orig_halo_id'] == hal)[0]
            if len(ii) > 0:
                idx_list.append(halos_now['id'][ii][0])
   
    complete_list=np.zeros(len(halo_list), dtype=bool)
    
    for i, this_hal in enumerate(idx_list):
        complete = False
        try:
            atree = extract_main_tree_full(treedata, this_hal)
            if min(atree['nout']) <= nout_ini and max(atree['nout']) >= nout_fi:
                complete = True
        except:
            pass
        complete_list[i] = complete

    return halo_list[complete_list]
    
    
    
