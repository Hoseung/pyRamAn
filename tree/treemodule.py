# -*- coding: utf-8 -*-
"""
treemodule

Created on Sun Jun 14 06:35:45 2015

@author: hoseung
"""
import numpy as np

def _is_ascii(filename):
    return filename.split(".")[-1] == "dat"
    
class CTree(object):
    """
        compatible with ConsistentTrees 1.01
    """

    def __init__(self, filename=None):
        if filename is not None:
            self.load(filename=filename)
    
    def _add_info(self):
        #self.pboxsize = 199.632011
        self.pboxsize = 200.0

    def _load_ascii(self, filename):
        cnt_header = 0
        
        datatype =[ 'f8','i8','i8','i8','i8','i8','i8','i8','i8','f8'\
                   ,'f8','f8','i8','f8','f8','f8','i8','f8','f8','f8'\
                   ,'f8','f8','f8','f8','f8','f8','f8','i8','i8','i8'\
                   ,'i8','i8','i8','i8','i8','f8','i8']#,\      


        with open(filename, 'rb') as f:   
            for i in range(180):
                line = f.readline()
                line = line.decode('utf-8')
                if line[0] != '#':
                    self.ntrees = int(line) # The first line after the header is the number of trees.
                    cnt_header = f.tell()
                    break
            f.seek(cnt_header)
            self.data = np.genfromtxt(f,dtype=datatype)
    
        self.data.dtype.names=(\
             'aexp','id','desc_aexp','desc_id','nprog','pid','upid','desc_pid','phantom','sam_mvir'\
            ,'mvir','rvir','rs','vrms','mmp','aexp_last_MM','vmax','x','y','z'\
            ,'vx','vy','vz','jx','jy','jz','spin','b_id','d_id','tree_root_id'\
            ,'Orig_halo_id','nout','next_coprogenitor_d_id','last_progenitor_d_id'\
            ,'last_mainleaf_depthfirst_id', 'tidal_force', 'tidal_id')

        print("Loading Consistent Tree data from ASCII is done")

    def _load_pickle(self, filename):
        import pickle
        try:
            with open(filename, "rb") as f:
                self.data = pickle.load(f)
                print("Loading Consistent Tree data from ASCII is done")
        except IOError:
            print("Error, No such file.", filename)
        
    def load(self, filename=None):
        if filename is None:
#            from tkinter import tk
#            tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
            from tkinter.filedialog import askopenfilename
            filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file

        if _is_ascii(filename):
            self._load_ascii(filename)
        else:
            self._load_pickle(filename)

        # The tree output file is written in 'bytes' string rather than string - this is the modern concept of text.
        # So the b' is at the very begining of lines. Python3 now distinguishes between string and byte string.
        # 
        # Python3 strings are unicode by default. 
        # You need to specify the encoding of the text file. 
        # Of course you can there are built-in methods to detect the encoding of a text file.
        
        # The output of Consistent tree is 'utf-8' 
        # Additional complexity is that
        # numpy genfromtxt always want it to be byte strings.
        # So you need 'rb' rather than 'r'
        self._tree_ids()
           
    def _tree_ids(self):
        i = np.where(self.data['nout'] == max(self.data['nout']))[0]
        self.trees_idx = self.data['id'][i]
        self.trees_id = self.data['Orig_halo_id'][i]

    """ get_main_prg in treeutils is a working version. 
        Use that instead.
        def get_main_prg(trees, haloid=None, haloind=None, unique_id=True):


    def get_main_prg(self, ids=None, original=True):
        if ids is None:
            if original:
                ids = self.trees_id
            else:
                ids = self.trees_idx

        if type(ids) is not list: ids = [ ids ]
        for thisid in ids:
            tree = self.data[np.where(self.data['tree_root_id'] == thisid)]
            prgs=[]
            next_h = tree[0]['id']
            while next_h != -1:
                print(next_h)
                i_hal = np.where(tree['id'] == next_h)
                halo = tree[i_hal]
                next_h = halo['last_progenitor_d_id']
                prgs.append(next_h)
    """                
            # Now, what are the main progenitors?

    def show_data(self, data, ind):
        '''
        This prints list of filed names and values of a numpy array in two column format.
        Sometimes there are too many fields that listing the names in a row makes it difficult
        to match the field name and the field value.
        Example)
        >>>
        '''
        for ii,jj in zip(self.data.dtype.names,data[ind]):
            print("%s : %f" % (ii,jj))

def load_tree(wdir, is_gal=False, no_dump=False, load_ascii=False):
    import pickle
    #from tree import treemodule
    import tree.ctutils as ctu
    from general import defaults
    
    df = defaults.Default()
    tree_path = df.tree_path(is_gal=is_gal)

    if load_ascii:
        alltrees = CTree()
        alltrees.load(filename= wdir + tree_path + 'tree_0_0_0.dat')
        # Fix nout -----------------------------------------------------
        nout_max = alltrees.data['nout'].max()
        alltrees.data['nout'] += 187 - nout_max
        print("------ NOUT fixed")
        alltrees.data = ctu.augment_tree(alltrees.data, wdir, is_gal=is_gal)
        print("------ tree data extended")
        if not no_dump:
            pickle.dump(alltrees, open(wdir + df.ext_tree_pickle(is_gal=is_gal), "wb" ))
    else:
        try:
            alltrees = pickle.load(open(wdir + df.ext_tree_pickle(is_gal=is_gal), "rb" ))
            print("Loaded an extended tree")
        except:
            alltrees = CTree()
            alltrees.load(filename= wdir + tree_path + 'tree_0_0_0.dat')
            # Fix nout -----------------------------------------------------
            nout_max = alltrees.data['nout'].max()
            alltrees.data['nout'] += 187 - nout_max
            print("------ NOUT fixed")
            alltrees.data = ctu.augment_tree(alltrees.data, wdir, is_gal=is_gal)
            print("------ tree data extended")
            if not no_dump:
                pickle.dump(alltrees, open(wdir + df.ext_tree_pickle(is_gal=is_gal), "wb" ), alltrees)
                
    return alltrees


def rs2codeunit(rst):
    """ 
        Check and clean up.         
        
        nout in Consistent Tree by default starts from 0 regardless of 
        the original simulation snapshot number. 
        This function assumes the nouts are already fixed. 
        In practice, it should be fixed when reading from ASCII and pickling it.
    
    """
    import numpy as np
    nouts = np.unique(rst['nout'])
#    wdir = '/home/hoseung/Work/data/AGN2/'
    for nout in nouts:
#        inow = np.where(rst['nout'] eq nout)
#        info = load.info.Info(base = wdir, nout=nout)

#%%
        dir_halo = base + "rhalo/rockstar_halos/"
        f_tree = base + "rhalo/tree.pickle"
        with open(f_tree, "rb") as ft:
            rstree = pickle.load(ft)

#inout = np.where(rstree['nout'] eq 81)
#print(rstree['x'][inout])

            rstree['x'] *= 200/199.632011
            rstree['y'] *= 200/199.632011
            rstree['z'] *= 200/199.632011

#print(hall.)    