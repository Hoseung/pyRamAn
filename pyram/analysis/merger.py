# -*- coding: utf-8 -*-
"""
Created on Wed May 11 02:22:44 2016

@author: hoseung
"""

import tree.ctutils as ctu
#from tree import treeutils
import numpy as np
import pickle

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
        idx = ctu.get_progenitors(atree, idx, main=True)[0]
        ind = np.where(atree['id'] == idx)[0]
        if atree['aexp'][ind] < aexp_min:
            break
        nprg = ctu.get_npr(atree, idx)

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
        
    prgs = ctu.get_progenitors(atree, idx)
    
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
    #id_prgs = atree['id'][i_prgs]
    mass_prgs = atree['m'][i_prgs]
    
    #mass_prgs_norm = mass_prgs / sum(mass_prgs)

    return mass_prgs


def load_tree(wdir, is_gal=False, no_dump=False):
    #import pickle
    from tree import treemodule
    import tree.ctutils as ctu

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
        alltrees.load(filename= wdir + tree_path + 'tree_0_0_0.dat')
        if not no_dump:
            # Fix nout -----------------------------------------------------
            nout_max = alltrees.data['nout'].max()
            alltrees.data['nout'] += 187 - nout_max
            print("------ NOUT fixed")
            alltrees.data = ctu.augment_tree(alltrees.data, wdir, is_gal=is_gal)
            print("------ tree data extended")
        
    return alltrees


