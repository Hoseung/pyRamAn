# -*- coding: utf-8 -*-

"""

From tree.fits file, select halos by mass cut.

Created on Wed Jan 21 11:45:42 2015

@author: hoseung
"""
import numpy as np

from utils.io_module import read_fortran, skip_fortran
import numpy as np
import time

def load_Tree(fname):
    with open(fname, 'rb') as f:
        # Header
        nsteps = read_fortran(f, dtype=np.int32, check=False)
        dummy = read_fortran(f, dtype=np.int32, check=False, n=nsteps)
        nhals_arr = dummy[:nsteps]
        nsubhals_arr = dummy[nsteps:]
        nall = sum(nhals_arr) + sum(nsubhals_arr)

        aexp_arr = read_fortran(f, dtype=np.float32, check=False, n=nsteps)
        omega_t_arr = read_fortran(f, dtype=np.float32, check=False, n=nsteps)
        age_univ_arr = read_fortran(f, dtype=np.float32, check=False, n=nsteps)

        dint = np.int32
        dfloat = np.float32
        ddouble = np.float64
        #ddoublee = np.float64

        dtype_tree = [('zred', '<f8'),
                      ('nstep', '<i4'), ('id', '<i4'),('m', '<f8'),
                      ('macc', '<f8'), ('nsub', '<i4'),# ('nfathers', '<i4'), ('nhalo', '<i4'),
                      ('rvir', '<f8'), ('mvir', '<f8'),('xp', '<f8', (3,)), ('vp', '<f8', (3,)),
                      ('level', '<i4'),#('listid', '<i4', (2,)),
                      ('hosthalo', '<i4'), ('hostsub', '<i4'),
                      ('nextsub', '<i4'), ('idx', '<i4'),
                      ('nprgs', '<i4'),# ('n_mergers', '<i4'),
                      ('flist_index', '<i4'), ('slist_index', '<i4')]

        t=np.zeros(nall, dtype = dtype_tree)
        #startid = np.zeros(nall, dtype=np.int32)
        #endid = np.zeros(nall, dtype=np.int32)

        # Main data
        BIG_RUN = True
        h0 = 70.4
        boxlen_ini = 100/(h0/100)

        t0 = time.time()

        idx = 0
        flist_index = 0
        slist_index = 0

        n_min_list = 1000

        fatherID=[]
        fatherIDx=[]
        fatherMass=[]
        sonID=[]

        cnt_nofather=0

        for i in range(nsteps):
            nhals_now = nhals_arr[i] + nsubhals_arr[i]
            #t[idx:idx + nhals_now]['nhalo'] = nhals_arr[i] + nsubhals_arr[i]
            #t[idx:idx + nhals_now]["boxsize"] = aexp_arr[i] * boxlen_ini
            t[idx:idx + nhals_now]['zred'] = 1/aexp_arr[i] -1
            t[idx:idx + nhals_now]["nstep"] = i
            idx_old = idx
            for j in range(nhals_now):
                #nhalo = nhals_arr[i] + nsubhals_arr[i]

                #halid = read_fortran(f, dtype=dint)
                t[idx]['idx'] = idx
                t[idx]['id'] = read_fortran(f, dtype=dint)
                bushID = read_fortran(f, dtype=dint)
                # What does st stand for?
                st = read_fortran(f, dtype=dint)
                t[idx]["level"], t[idx]["hosthalo"], t[idx]["hostsub"],                t[idx]["nsub"], t[idx]["nextsub"]= read_fortran(f, dtype=dint, n=5)
                # Mass in 1e11 solar mass.
                t[idx]['m'] = read_fortran(f, dtype=dfloat) * 1e11
                t[idx]['macc'] = read_fortran(f, dtype=ddouble)

                t[idx]["xp"] = read_fortran(f, dtype=dfloat, n=3)
                t[idx]["vp"] = read_fortran(f, dtype=dfloat, n=3)
                lp = read_fortran(f, dtype=dfloat, n=3)
                abc = read_fortran(f, dtype=dfloat, n=4)
                ek, ep, et = read_fortran(f, dtype=dfloat, n=3)
                spin = read_fortran(f, dtype=dfloat)
                n_fathers = read_fortran(f, dtype=dint)[0]

                #t[idx]['nhalo'] = nhalo
                #t[idx]['hid'] = halid

                t[idx]['nprgs'] = n_fathers

                #endid[idx] = startid[idx] + n_fathers - 1
                #if idx != nall -1 :
                #    startid[idx + 1] = endid[idx] + 1

                # every halo has at least one father, the background.
                #if n_fathers > 0:
                fid = read_fortran(f, dtype=dint, n=n_fathers)
                fatherMass.append(read_fortran(f, dtype=dfloat, n=n_fathers))
                #if n_fathers > 1:
                fatherID.append(fid.copy())
                # fid[:] copies the value
                # Otherwise, the fatherID entry will be modified afterwards.
                t[idx]["flist_index"] = flist_index
                flist_index += 1

                # Keep fid==0.
                # They are background, but fatherMass has corresponding values.
                # And they are not always the first or last element.
                # So I need to know which values are for background.
                if i > 0 :
                    fid[fid>0] = t_before["idx"][fid[fid > 0]-1]
                fatherIDx.append(fid)

                    #tmp_father_mass = read_fortran(f, dtype=dfloat, n=n_fathers)
                #else:
                #    #cnt_nofather +=1
                #    print("No father")

                nsons = read_fortran(f, dtype=dint)
                if nsons > 0:
                    sonID.append(read_fortran(f, dtype=dint, n=nsons))
                    t[idx]["slist_index"] = slist_index
                    slist_index += 1

                rvir, mvir, tvir, cvel = read_fortran(f, dtype=dfloat, n=4)
                mvir = mvir * 1e11
                t[idx]['mvir'] = mvir
                t[idx]["rvir"] = rvir

                rho_0, rho_c = read_fortran(f, dtype=dfloat, n=2)
                # i-dependent
                if not BIG_RUN:
                    t[idx]["np"] = read_fortran(f, dtype=dint)

                idx = idx + 1

            # Keep tree at the previous snapshot
            # to make fatherIDx list.
            t_before = t[idx_old:idx_old + nhals_now]

    print("Took", time.time() - t0)
    return t, fatherID, fatherIDx, fatherMass



def load_fits(base=None, work_dir=None, filename="halo/TMtree.fits"):
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
