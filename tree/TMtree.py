# -*- coding: utf-8 -*-

"""


From tree.fits file, select halos by mass cut.


Created on Wed Jan 21 11:45:42 2015

@author: hoseung
"""
def load(work_dir="/home/hoseung/Work/data/", filename="tree.fits"):
    from astropy.io import fits
    from astropy.table import Table
    data = fits.getdata(work_dir + filename, 1)
    return Table(data)


def get_main_prg(tree, halos, nout_ini=None, nout_fi=None):
#
#   nout in tree =/= nout. correct it.
#
#
    import numpy as np
    from collections import Iterable
    if not isinstance(halos, Iterable):
        halos = [halos]
    inout = np.where(tree[:]["NOUT"] == 0)[0]
    prg_list = np.zeros([len(halos), abs(nout_ini - nout_fi) + 1], dtype=int)
    hnu_list = np.zeros([len(halos), abs(nout_ini - nout_fi) + 1], dtype=int)

    for ihalo, halo in enumerate(halos):
        i_prg = np.where(tree[inout]["HALNUM"] == halo)[0]
        prg_idx = tree[inout[i_prg]]['IDX'][0]
        hnu_list[ihalo][0] = halo
#        print(prg_idx)
        prg_list[ihalo][0] = prg_idx
        for i in range(nout_fi, nout_ini):
        # index = idx
            prg_idx = tree[prg_idx]["TREE"][0]
            prg_list[ihalo][i + 1] = prg_idx
            hnu_list[ihalo][i + 1] = tree[prg_idx]["HALNUM"]
#            print(prg_idx, tree[prg_idx]["NOUT"])
        # First element of "tree" is the main progenitor.
        # I assume....

    return prg_list, hnu_list

def check_tree_complete(tree, nout_fi, nout_ini, halo_list):
    import numpy as np
    from tree import TMtree
    '''
    returns a list of halo IDs at nout_fi that with trees fully linked
    from nout_ini to nout_fi.
    '''
    # Make sure all parameters are given
    complete_list=np.zeros(len(halo_list), dtype=bool)
    #for ihal,halo in enumerate(halo_list):
    idxlist, idlist=TMtree.get_main_prg(tree, halo_list, nout_fi=0, nout_ini=nout_ini)
    for i in range(len(halo_list)):
#        print('aa',idxlist[i])
#        print(len(np.where(idxlist[i] == 0)[0]))
        if len(np.where(idxlist[i] == 0)[0]) > 1:
            complete = False
            #print(i,'Broken')
        elif len(np.where(idxlist[i] == 0)[0]) == 0:
            #print(i,'Complete')
            complete = True
        complete_list[i] = complete

    return complete_list, idlist[complete_list]