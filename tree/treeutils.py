
"""
Created on ???

This routine makes plots to test CT outputs.
Test halo porperty consistency.

Compatible with Python3

Commonly used routines processing CT output.

@author: hoseung
"""
import numpy as np

def iscomment(s):
	return s.startswith('#')

def load_ascii(filename = None):
    '''
    Loading ~500MB takes ~ 20 seconds, most of which is suspected to be
    processing part.

    '''
    if filename is None :
        from tkinter import Tk
        # Note that the name has changes from Tkinter (python2) to tkinter (python3)
        from tkinter.filedialog import askopenfilename
        Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
        filename = askopenfilename()
# show an "Open" dialog box and return the path to the selected file

    datatype =[ 'f8','i8','i8','i8','i8','i8','i8','i8','i8','f8',\
                'f8','f8','i8','f8','f8','f8','i8','f8','f8','f8',\
                'f8','f8','i8','i8','i8','f8','i8','i8','i8','i8',\
                'i8','i8','i8','f8','i8','i8','i8','i8','i8','i8',\
                'f8','f8','f8','f8','f8','f8','f8','f8','f8','f8',\
                'f8','f8','f8','f8','i8','i8']

    names=['aexp','id','desc_aexp','desc_id','nprog','pid','upid','desc_pid',
           'phantom','sam_mvir','mvir','rvir','rs','vrms','mmp','aexp_last_MM',
           'vmax','x','y','z','vx','vy','vz','jx','jy','jz','spin','b_id',
           'd_id','tree_root_id','Orig_halo_id','nout',
           'next_coprogenitor_d_id','last_progenitor_d_id','rs_Klypin',
           'mvir_all','m200b','m200c','m500c','m2500c','xoff','voff',
           'spin_Bullock','btoc','ctoa','ax','ay','az','btoa500','ctoa500',
           'ax500','ay500','az500','toveru','mpe_b','mpe_d']

    cnt = 0
    with open(filename, 'rb') as f:
        for i in range(180):
            line =f.readline().decode("utf-8")
            cnt += 1
            if not iscomment(line):
                f.readline() # one line
                break
            # data start after the last '#' mark
    data = np.genfromtxt(filename, dtype=datatype, names=names, skip_header=cnt)

    print('DONE')
    return data


def tree2pickle(data):
    import pickle
    with open("tree.pickle", "wb") as f:
        pickle.dump(data, f, protocol = 4)


def show_data(data,ind):
    '''
    This prints list of filed names and values of a numpy array in two column format.
    Sometimes there are too many fields that listing the names in a row makes it difficult
    to match the field name and the field value.
    Example)
    >>>
    '''
    for ii,jj in zip(data.dtype.names,data[ind]):
        print("%s : %f" % (ii,jj))


def get_field_view(arr, fields):
    """
    Return a view to a 2D numpy array with multiple fields.
    By default, veiw = ndarray['field'] returns a view,
    but copeid_array = ndarray['field1','field2'] copies colum data.
    Example:
    >>>
    """
    dtype2 = np.dtype({name:arr.dtype.fields[name] for name in fields})
    return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)


def get_main_prg_tree(all_tree, haloID, nout_ini=None, nout_fi=None):
    """
    Return a sub-tree recarray(or ndarray) of the bigger tree,
    which is a tree of main progeniotrs of the halo, haloID.
    """
    from utils import match

    idx_list = get_main_prg_ID(all_tree, haloID, nout_ini=nout_ini, nout_fi=nout_fi)
    ind = match.match_list_ind(all_tree['id'], idx_list)
    return all_tree[ind]


def richest_halo(all_tree, nout=None):
    """ Returns the ID of the richest halo in a given nout.
        Useful to find the zoomed-in target cluster.
    """
    pass

def get_main_prg_ID(all_tree, haloID, nout_ini=None, nout_fi=None,
                 aexp_ini=None, aexp_fi=None):
    """
    Return main progenitor IDs between nout_ini and nout_fi.

    haloID   -- IDx
    nout_ini -- first snapshot of the tree
    nout_fi  -- last snapshot of the tree

    If not given, min & max nout of all_tree is used.

    """
    import numpy as np

    if (nout_ini is None) & (aexp_ini is None):
        aexp_ini = min(all_tree['aexp'])
        nout_ini = min(all_tree['nout'])

    if (nout_fi is None) & (aexp_fi is None):
        aexp_fi = max(all_tree['aexp'])
        nout_fi = max(all_tree['nout'])

    haloIDx = all_tree[np.where((all_tree['nout'] == nout_fi)
                    & (all_tree['tree_root_id'] == haloID))]['id']

    tree=np.zeros(nout_fi - nout_ini + 1,dtype=np.int64)
    tree.fill(-1)

    if len(haloIDx) == 0:
        print("There is no tree for {0}".format(haloID))
    else:
        this_tree = all_tree[all_tree['tree_root_id'] == haloIDx]

        desc_idx = haloIDx
        tree[0] = desc_idx
        for i in range(nout_ini, nout_fi):
            prog_ind = np.where((this_tree['desc_id'] == desc_idx) & (this_tree['mmp'] == 1))
            if len(prog_ind[0]) == 0:
                print("Tree is broken at {0}-th nout".format(i))
                break
            desc_idx = this_tree['id'][prog_ind]
#            print(prog_ind, len(prog_ind), this_tree['Orig_halo_id'][prog_ind])
            tree[i+1] = desc_idx

    return tree

def get_main_prgTM(tree, halos, nout_ini=None, nout_fi=None):
    """
    TM version.
    nout in tree =/= nout. correct it.


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

    return prg_list, hnu_list


def get_main_prg(trees, haloids=None, haloinds=None, unique_id=True,
                 nout_ini=None, nout_fi=None):
    """
    returns list of main progenitors of the given halo
    and indices of the halos in the tree data.

    parameters
    ----------
    id : int
        list of halos (Tree ID by default)
    	tree_root_id is a unique id. So if a non-unique id is given,
    	convert it to a unique id.
    unique_id : bool
        False means that the input id is halo id, not a unique tree id.
    	If no 'ROOT' galaxy is found, return false


    Examples
    --------
    >>> inds_arr = np.zeros((n_good, 82), dtype=int)
    >>> i_good = 0
    >>> for hid in halo_list:
    >>>     prgs, inds = tru.get_main_prg(trees, haloid=hid, unique_id=False)
    >>>     if prgs is False:
    >>>        print("No tree for {}".format(hid))
    >>>        continue
    >>>     inds_arr[i_good] = inds
    >>>     i_good +=1

    """
    import numpy as np

    if nout_fi is None:
        nout_fi = max(trees.data['nout'])
    if nout_ini is None:
        nout_ini = min(trees.data['nout'])

    prg_arr = np.zeros([len(haloids), abs(nout_ini - nout_fi) + 1], dtype=int)
    ind_arr = np.zeros([len(haloids), abs(nout_ini - nout_fi) + 1], dtype=int)

    for i, haloid in enumerate(haloids):
        # convert to a unique ID
        if haloid is None:
            # select halo by index of the final snapshot tree output.
            haloid = trees.trees_idx[haloinds[i]] # idx is unique id.
        else:
            if unique_id is False:
                haloid = trees.trees_idx[np.where(trees.trees_id == haloid)]

        itt = np.where(trees.data['tree_root_id'] == haloid)[0]

        if len(itt) == 0 :
            return False, False
        tt = trees.data[itt]

        prgs=[]
        prg_inds=[]
        i_hal = [0]
        this_halo = tt['id'][i_hal[0]]
        ind_offset = min(np.where(trees.data['tree_root_id'] == haloid)[0])
        while True:
            prgs.append(this_halo)
            prg_inds.append(i_hal[0])
            i_hal = np.where(tt['desc_id'] == this_halo)[0]
            if len(i_hal) == 0: break
            this_halo = tt['id'][i_hal[0]]
        inds = np.asarray(prg_inds)
        inds += ind_offset
        prg_arr[i,0:len(inds)] = inds
        ind_arr[i,0:len(prgs)] = prgs

    return prg_arr, ind_arr

def check_tree_complete(tree, halo_list, nout_ini=None, nout_fi=None):
    import numpy as np
    '''
    returns a list of halo IDs at nout_fi that with trees fully linked
    from nout_ini to nout_fi.

    example
    -------
    >>> roots = [1630412, 1630413, 1630414, 1630415, 1630416]
    >>> h_ind_ok, halo_ok = tru.check_tree_complete(tt, roots[5:10],
                                            nout_ini = 0, nout_fi =81)
    >>> print(halo_ok.shape)
    >>> (5, 82)

    -------

    '''
    if nout_fi is None:
        nout_fi = max(tree.data['nout'])
    if nout_ini is None:
        nout_ini = min(tree.data['nout'])

    # Make sure all parameters are given
    complete_list=np.zeros(len(halo_list), dtype=bool)
    #for ihal,halo in enumerate(halo_list):
#    idxlist, idlist = treeutils.get_main_prg(tree, halo_list, nout_fi=0, nout_ini=nout_ini)

    # idlist must be in 2-d format even only one halo is given.
    idxlist, idlist = get_main_prg(tree, haloids=halo_list,
                                   nout_ini=nout_ini, nout_fi=nout_fi)

    # If no tree branch is found, false returned.
    if idxlist is False:
        return False, False
    else:
        for i in range(len(halo_list)):
            if len(np.where(idxlist[i] == 0)[0]) > 1:
                complete = False
                #print(i,'Broken')
            elif len(np.where(idxlist[i] == 0)[0]) == 0:
                #print(i,'Complete')
                complete = True
            complete_list[i] = complete

    return complete_list, idlist[complete_list]


def final_halo_list(data,nout=None):
    '''
    Returns list of galaxies at the final snapshot by default,
    or from the desired nout.
    Example:
    >>>
    '''
    if nout not in locals():
        try:
            nout = max(data['nout'])
            ii = data['nout'] == nout
            return data['id'][ii] # RS ... or standard
        except:
            nout = max(data['NOUT'])
            ii = data['NOUT'] == nout
            return data['HALNUM'][ii] # HM



def tm2ct(tm):
    import numpy as np
    dtype=[('aexp', '<f8'), ('id', '<i8'), ('desc_aexp', '<i8'),
           ('desc_id', '<i8'), ('nprog', '<i8'), ('pid', '<i8'),
 ('mvir', '<f8'), ('rvir', '<f8'), ('rs', '<i8'), ('vrms', '<f8'),
 ('mmp', '<f8'), ('aexp_last_MM', '<f8'), ('vmax', '<i8'),
 ('x', '<f8'), ('y', '<f8'), ('z', '<f8'),
 ('vx', '<f8'), ('vy', '<f8'), ('vz', '<i8'),
 ('jx', '<i8'), ('jy', '<i8'), ('jz', '<f8'),
 ('spin', '<i8'), ('b_id', '<i8'), ('d_id', '<i8'), ('tree_root_id', '<i8'),
 ('Orig_halo_id', '<i8'), ('nout', '<i8'), ('next_coprogenitor_d_id', '<i8'),
 ('last_progenitor_d_id', '<f8'), ('ax', '<f8'), ('ay', '<f8'), ('az', '<f8'),
 ('tree', '<i8', 50)]


#    newtree = np.zeros(len(hm), dtype=dtype)
    tm['id'] = tm.pop('HALNUM')
    tm['nout'] = tm.pop('NOUT')

    return newtree
