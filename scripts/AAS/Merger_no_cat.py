
# coding: utf-8

# In[3]:

import tree.ctutils as ctu
from tree import treeutils
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
    
    Mass ratio should refer to the time before any stripping happen.
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
    id_prgs = atree['id'][i_prgs]
    mass_prgs = atree['m'][i_prgs]
    
    #mass_prgs_norm = mass_prgs / sum(mass_prgs)

    return mass_prgs


def load_tree(wdir, is_gal=False, no_dump=False):
    import pickle
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


# In[30]:

import utils.match as mtc
import matplotlib.pyplot as plt
import pandas as pd

is_gal = True
nout_fi = 187

#Last merger
import matplotlib.pyplot as plt
nout_ini = 100 # recent merger = nout = 140 or later.

# Load tree
is_gal = True

# all catalogs
verbose=False

#
most_recent_only = False

#clusters = ['39990', '36415', '10002', '05427', '36413', '01605']
#clusters=['05427']

# final result arrays

gal_list=[]
mr_list=[]
nout_list=[]

#for cluster in clusters:
#wdir = '/home/hoseung/Work/data/' + clusters[0] + '/'
wdir = './'

alltrees = load_tree(wdir, is_gal=is_gal)

ft = alltrees.data[alltrees.data['nout'] == nout_fi]
allgals = ft['id'][ft['m'] > 1e10]

#catalog = pickle.load(open(wdir + '/catalog_GM/' + 'catalog187.pickle', 'rb'))
#print("Cluster ",cluster)
for idx in allgals:
    #gal = cat['id']
    #if verbose: print("analyzing merger events of galaxy ", gal)

    # Convert halo id to tree id
    #idx = id2idx(alltrees.data, gal, 187)
    #idx = cat['idx']

    # full tree of a galaxy
    atree = ctu.extract_a_tree(alltrees.data, idx)

    # main progenitor tree
    main = ctu.extract_main_tree(alltrees.data, idx)

    x_nout = main['nout'].flatten()
    x_nout = x_nout[x_nout > nout_ini]

    mass_ratios_single = np.zeros(len(x_nout))
    for i, nout in enumerate(x_nout):
        # merger ratio
        i_prgs = np.where(atree['desc_id'] == main['id'][i])[0]

        # multiple prgs = merger
        if len(i_prgs) > 1:
            if verbose: print(" {} mergers at nout = {}".format(len(i_prgs), nout))
            id_prgs = atree['id'][i_prgs]
            mass_prgs = atree['m'][i_prgs]
            m_r = mass_prgs / max(mass_prgs)
            if verbose:
                print(" Mass ratios : ", m_r)
            mass_ratios_single[i] = max([mass_prgs[1:] / max(mass_prgs)][0])
        else:
            mass_ratios_single[i] = 0

    ind_ok = np.where(mass_ratios_single > 0.1)[0]
    #print("all ind_ok", ind_ok)
    if len(ind_ok) > 0:
        # if a satellite oscillates around the host, 
        # it could be identified as multiple mergers with short time interval. 
        # leave only the first passage / merger.
        good =[]
        for i in range(len(ind_ok)-1):
            if ind_ok[i+1] > ind_ok[i] + 2:
                good.append(ind_ok[i])
        good.append(ind_ok[-1])
        ind_ok = good
#        if most_recent_only:
#            ind_ok = max(ind_ok) # most recent 

#        print("  galaxy {}, Last nout {}, Merger ratio 1:{:.1f}".format(idx,
#                                                                     x_nout[ind_ok],
#                                                                       1./mass_ratios_single[ind_ok]))
        mr = 1./mass_ratios_single[ind_ok]

        gal_list.append(idx)
        mr_list.append(mr)
        nout_list.append(x_nout[ind_ok])


        
"""
fig, ax = plt.subplots(1)

ax.scatter(nout_list, mr_list)
ax.set_title("last merger Vs final lambda")
ax.set_ylabel(r"$\lambda _R$")
ax.set_xlabel("Last merger")
for i,gal_name in enumerate(gal_list):
    ax.text(nout_list[i]+0.5, mr_list[i]+0.1, str(gal_name))
plt.show()
"""
with open(wdir + 'merger_list.txt', 'w') as f:
#    print("Major mergers in this cluster")
    for gal, nout, mr in zip(gal_list, mr_list, nout_list):
        for ni, mi in zip(nout, mr):
            f.write("{}  {}  {} \n".format(gal, ni, mi))
        

