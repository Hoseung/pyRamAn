# -*- coding: utf-8 -*-
"""
dump id lists of target clusters
ID.dat files are fed into target_clusters_3D.py 
to visualize cluster selection in 3D.

*Complete*

Created on Wed Jun 24 02:38:16 2015

@author: hoseung
"""
#%%
def load_hm_full(nout, halolist, base='./'):
    from load import utils
    import numpy as np

    with open(base + 'halo/tree_bricks080', "rb") as f:
        # header
#        nbodies = utils.read_fortran(f, np.dtype('i4'), 1)
#        massp = utils.read_fortran(f, np.dtype('f4'), 1)
        utils.read_fortran(f, np.dtype('i4'), 1)
        utils.read_fortran(f, np.dtype('f4'), 1)
        utils.read_fortran(f, np.dtype('f4'), 1)
        utils.read_fortran(f, np.dtype('f4'), 1)
        utils.read_fortran(f, np.dtype('f4'), 1)
        halonum, subnum = utils.read_fortran(f, np.dtype('i4'), 2)

        # halo data
        for i in range(halonum + subnum):
            npart = utils.read_fortran(f, np.dtype('i4'), 1)
            allid = utils.read_fortran(f, np.dtype('i4'), npart)
            hnu = utils.read_fortran(f, np.dtype('i4'), 1)
            if hnu in halolist:
                print("write to file", hnu)
                with open(base + 'haloIDs'+ str(hnu[0]).zfill(5)+'.dat', 'wb') as fout:
                    allid.tofile(fout)

            utils.read_fortran(f, np.dtype('i4'), 1)
            utils.read_fortran(f, np.dtype('i4'), 5)
            utils.read_fortran(f, np.dtype('f4'), 1)
            utils.read_fortran(f, np.dtype('f4'), 3)
            utils.read_fortran(f, np.dtype('f4'), 3)
            utils.read_fortran(f, np.dtype('f4'), 3)
            utils.read_fortran(f, np.dtype('f4'), 4)
            utils.read_fortran(f, np.dtype('f4'), 3)
            utils.read_fortran(f, np.dtype('f4'), 1)
            utils.read_fortran(f, np.dtype('f4'), 4)
            utils.read_fortran(f, np.dtype('f8'), 1)

"""
            nstep = utils.read_fortran(f, np.dtype('i4'), 1)
            ha, hb, hc, hd, he = utils.read_fortran(f, np.dtype('i4'), 5)
            mhalo = utils.read_fortran(f, np.dtype('f4'), 1)
            xc, yc ,zc = utils.read_fortran(f, np.dtype('f4'), 3)
            vx, vy, vz = utils.read_fortran(f, np.dtype('f4'), 3)
            jx, jy, jz = utils.read_fortran(f, np.dtype('f4'), 3)
            r,a,b,c = utils.read_fortran(f, np.dtype('f4'), 4)
            ke, pe, etot = utils.read_fortran(f, np.dtype('f4'), 3)
            spin = utils.read_fortran(f, np.dtype('f4'), 1)
            rvir,mvir,tvir,cvel = utils.read_fortran(f, np.dtype('f4'), 4)
            utils.read_fortran(f, np.dtype('f8'), 1)
"""


import tree
import numpy as np
nout=80
wdir = '/home/hoseung/Work/data/DMO/'
# load halo (HM)
hall = tree.halomodule.Halo(nout=nout, base=wdir, halofinder="HM")


#%%
with open(wdir+'allTargets.txt', 'r') as fhlist:
    halolist = [int(i.split('\n')[0]) for i in fhlist.readlines()]
    


#%%
# Read target cluster list (from an ASCII file)

load_hm_full(nout, halolist, base=wdir)
print("done")
