# -*- coding: utf-8 -*-

"""

From tree.fits file, select halos by mass cut.

Created on Wed Jan 21 11:45:42 2015

@author: hoseung
"""
import numpy as np
def load(base=None, work_dir=None, filename="halo/TMtree.fits"):
    """
        Load tree.fits file and get the data table.
        use base always. Work_dir is for backward compatibility.
    """
    from astropy.io import fits
    from astropy.table import Table
    if base == None:
        base = work_dir
        
    if base == None:
        raise ValueError ("Please provide the parameter base=")
        
    data = fits.getdata(base + filename, 1)
    return Table(data)


def fix_nout(tt, nout_ini, nout_fi):
    nout_min_org = tt['NOUT'].min()
    nout_max_org = tt['NOUT'].max()

    nnout_org = nout_max_org - nout_min_org + 1
    nnout_new = nout_fi - nout_ini + 1
    
    assert (nnout_org == nnout_new), "number of nouts does not match"
    
    i_z_max = np.where(tt["Z"] == tt["Z"].max())[0]
    assert (tt["NOUT"][np.where(tt["Z"] == tt["Z"].max())[0]][0] == nout_max_org), "The snapshot number of highest redshift snapshot is not 0. Maybe you've already fixed nouts?"
    
    # OK. It's safe.
    tt["NOUT"] = nout_fi - tt["NOUT"]   
    

def get_main_prg(tree, halos, nout_ini=None, nout_fi=None):
    """
    TM version.
    nout in tree =/= nout. correct it.

    parameters
    ----------
    tree 
        tree object (recarray)
    halos
        ID(s) of target halo(s) at the final nout.

    """
    import numpy as np
    from collections import Iterable
    if not isinstance(halos, Iterable):
        halos = [halos]
    i_nout_final = np.where(tree[:]["NOUT"] == 0)[0]
    if nout_ini is None:
        nout_ini = max(tree['NOUT'])
    if nout_fi is None:
        nout_fi = min(tree['NOUT'])
        
    prg_list = np.zeros([len(halos), abs(nout_ini - nout_fi) + 1], dtype=int)
    hnu_list = np.zeros([len(halos), abs(nout_ini - nout_fi) + 1], dtype=int)

    for ihalo, halo in enumerate(halos):
        i_prg = np.where(tree[i_nout_final]["HALNUM"] == halo)[0]
        prg_idx = tree[i_nout_final[i_prg]]['IDX'][0]
        hnu_list[ihalo][0] = halo
        
#        print(prg_idx)
        prg_list[ihalo][0] = prg_idx
        for i in range(nout_fi, nout_ini):
            prg_idx = tree["TREE"][prg_idx][0]
            prg_list[ihalo][i + 1] = prg_idx
            hnu_list[ihalo][i + 1] = tree["HALNUM"][prg_idx]
#            print(prg_idx, tree[prg_idx]["NOUT"])
        # First element of "tree" is the main progenitor.

    return prg_list, hnu_list

def check_tree_complete(tree, nout_fi, nout_ini, halo_list):
    import numpy as np
    '''
    returns a list of halo IDs at nout_fi that with trees fully linked
    from nout_ini to nout_fi.
    
    '''
    # Make sure all parameters are given
    complete_list=np.zeros(len(halo_list), dtype=bool)
    #for ihal,halo in enumerate(halo_list):
    idxlist, idlist=get_main_prg(tree, halo_list, nout_fi=0, nout_ini=nout_ini)
    for i in range(len(halo_list)):
        if len(np.where(idxlist[i] == 0)[0]) > 1:
            complete = False
        elif len(np.where(idxlist[i] == 0)[0]) == 0:
            complete = True
        complete_list[i] = complete
    
    return complete_list, idlist[complete_list]