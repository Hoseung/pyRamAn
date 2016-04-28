# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 11:00:31 2015

halo / galaxy calss including basic data load functionality. 

Parameters
----------
nout : int
    snapshot number
base : str
    work directory
is_gal : logical
    load galaxy if true
return_id : logical
    If True, load and return constituent particle id of each halo. 
return_id_list : int list
    specify halos of which particle id is returned.
    If return_id is True but return_id_list is None, 
    particle id of all halos are returned.


MODIFICATIONS:
2015. 08. 08
    Floats are all stored as doubles even though the original data is float32.
    If center position is stored in float32, then is used to select particles
    inside the region, then max(x) - min(x) > 2 * region['radius'] is possible!

2016.03.26
    if return_id = True, but no return_id_list is set, then return id lists
     of all halos by dafault.

@author: hoseung
"""

class HaloMeta():
    """HaloMeta class.
    all except halofinder output.

 
    Parameters
    ----------
    nout : snapshot number, int
    base : working directory, str
    info : info object
    halofinder : "RS" or "HM". Full names also work. case insensitive.
 
    Examples
    --------
    It is better to specify nout, base, halofinder from the beginning. 
    All three are necessary to load a halo output
    
    >>> h = tree.halomodule.Halo(nout=132, halofinder="RS", base='~/data/AGN2/')    
    
    given nout and base, info is autoloaded if not explicitely given.
    
    """
    def __init__(self, nout=None, base=None, info=None, halofinder='HM',
                 load=False, is_gal=False, return_id=False, outdir=None,
                 vebose=False):
       self.nout = nout
       self.verbose = verbose
       self.set_base(base)
       self.info = None
       self.set_info(info=info)
       self.run_params = {"none":None}
       if isinstance(halofinder, str):
           self.set_halofinder(halofinder)
       self.aexp = 0.0
       self.age = 0.0
       self.sim = {"none":None}
       self.nhalo = 0
       self.nsub = 0
       self.npart_all = 0
       self.massp = 0 # in case of single level DMO run.
       self.unit={"mass":None, "lengh":None, "velocity":None}
       self.is_gal = is_gal
       self.return_id = return_id
       self.return_id_list=None

       if outdir is None:
           if is_gal:
               self.gal_find_dir = 'GalaxyMaker/'
           else:
               self.dm_find_dir= 'halo/'
       else:
           if is_gal:
               self.gal_find_dir = outdir
           else:
               self.dm_find_dir= outdir
       
       try:
           self.set_nout(nout)
       except:
           pass
       if load:
           self.load()
       
    def set_info(self, info):
        if info is None:
            try:
                self.load_info()
            except:
            print("[Halo.set_info] Are these parameters right?")
        else:
            self.info = info

    def load_info(self, nout=None):
        import load.info
        if nout is None:
            nout = self.nout
        if self.verbose : print("[Halo.load_info] loading info")
#        print("[Halo.load_info]", nout, self.base)
        self.info = load.info.Info(nout=nout, base=self.base, load=True)    
        if self.verbose : print("[Halo.load_info] info is loaded")

    def set_halofinder(self, halofinder):
        from utils import util
        options = ['HM', 'halomaker', 'rockstar', 'rs']
        guess =  util.fuzzymatch(halofinder, answer_list=options)
        if guess in ["HM", 'halomaker']:
            answer = "HaloMaker"
        elif guess in ['rockstar', 'rs']:
            answer = 'Rockstar'
        else:
            print("Don't understand what you mean:" + halofinder)
            answer = None
        self.halofinder = answer

    def set_base(self, base):
        self.base = base

    def set_sim(self, sim):
        # nouts and aexps
       pass

    def set_nout(self, nout):
        self.nout = nout
        self._set_aexp(nout)

    def _set_aexp(self, nout):
#        self.aexp =
        pass       
        

class Halo(HaloMeta):
    """
    HM :
        r = distance of the furthest particle from the center of the halo.
        rvir = radius where the density is 200 times the _______?
    Notes
    -----
       Rockstar halo id starts from 0. CT id too.

    """
#    def __init__(self):  __ini__ is also inherited from HaloMeta

    def set_data(self, data):
        if data is None:
            self.load()
        else:
            self.data = data

    def load(self):
        if self.halofinder is 'Rockstar':
            self.load_rs()
        elif self.halofinder is 'HaloMaker':
            if self.return_id:
                self.load_hm_old()
            else:
                self.load_hm()
            if self.info is None:
                import load
                info = load.info.Info(base = self.base, nout = self.nout, load = True)
                self.set_info(info)
        
            self.normalize()

    def load_hm_sav(self, nout=None, base=None, info=None):
        if nout is None:
            nout = self.nout
        if base is None:
            base = self.base
        snout = str(self.nout).zfill(3)
        try:
            self.data = readsav(base + 'halo/halo' + snout + '.sav')['h']
            from scipy.io.idl import readsav
        except:
            print("Cannot specify a file to read")
            print("trying to read {0}".format(base + 'halo/halo' + snout + '.sav'))
            print("+++++++++++++++++++++++")

        self.refactor_hm()
        if info is not None:
            self.set_info(info)
#        print("load done")
        self.normalize()
        # load .sav file
    
    def set_return_id_list(self, ids):
            self.return_id_list = ids

    def load_hm(self, nout=None, base=None, info=None):
        if nout is None:
            nout = self.nout
        if base is None:
            base = self.base
        snout = str(self.nout).zfill(3)
        if self.is_gal:
            fn = base + self.gal_find_dir + 'gal/tree_bricks' + snout
        else:
            fn = base + self.dm_find_dir + 'DM/tree_bricks' + snout
        try:
            dtype_halo = [('np', '<i4'), ('id', '<i4'), ('level', '<i4'),
                          ('host', '<i4'), ('sub', '<i4'), ('nsub', '<i4'),
                        ('nextsub', '<i4'),
                        ('m', '<f4'), ('mvir', '<f4'),
                        ('r', '<f4'), ('rvir', '<f4'), 
                        ('tvir', '<f4'), ('cvel', '<f4'), 
                        ('x', '<f4'), ('y', '<f4'), ('z', '<f4'),
                        ('vx', '<f4'), ('vy', '<f4'), ('vz', '<f4'),
                        ('ax', '<f4'), ('ay', '<f4'), ('az', '<f4'),
                        ('sp', '<f4'), ('idx', '<i4'),
                        ('p_rho', '<f4'),('p_c', '<f4'), 
                        ('energy', '<f8', (3,)), ('radius', '<f8', (4,))]

            if self.is_gal:
                dtype_halo += [('sig', '<f4'), ('sigbulge', '<f4'),
                               ('mbulge', '<f4'), ('hosthalo', '<i4'),
                               ('g_nbin', '<i4'), ('g_rr', '<f4', (100,)),
                               ('g_rho', '<f4', (100,))]

            import tree.load_c.rd_hal as rd_halo
            import numpy as np
            temp = rd_halo.read_file(fn.encode(), int(self.is_gal))# as a byte str.
            
            self.nbodies, self.halnum, self.subnum,\
                self.massp, self.aexp, self.omegat, self.age = temp[0:7]
            ntot = self.halnum + self.subnum
            self.data = np.recarray(ntot, dtype=dtype_halo)
            self.data['np'], self.data['id'],\
            levels, ang, energy, \
            self.data['m'],\
            radius, pos,\
            self.data['sp'], vel = temp[7:17]
            vir, profile = temp[17:19]
            if self.is_gal:
                temp_gal = temp[19]

            self.data['energy'] = energy.reshape((ntot,3))
            levels = levels.reshape((ntot,5))
            self.data['level'], self.data['host'], \
            self.data['sub'], self.data['nsub'], \
            self.data['nextsub'] = levels[:,0], \
                levels[:,1], levels[:,2],levels[:,3], levels[:,4]
            pos = pos.reshape((ntot,3))
            self.data['x'], self.data['y'], self.data['z'] \
                    = pos[:,0],pos[:,1],pos[:,2]
            vel = vel.reshape((ntot,3))
            self.data['vx'], self.data['vy'], self.data['vz'] \
                    = vel[:,0],vel[:,1],vel[:,2]
            ang = ang.reshape((ntot,3))
            self.data['ax'], self.data['ay'], self.data['az'] \
                    = ang[:,0],ang[:,1],ang[:,2]
            self.data['r'] = radius[::4].copy()

#            copy so that memory are continuous. (right?)
            self.data['rvir'],self.data['mvir'], \
                    self.data['tvir'],self.data['cvel'] = vir[::4].copy(),\
                vir[1::4].copy(),vir[2::4].copy(),vir[3::4].copy()
            self.data['p_rho'], self.data['p_c'] =\
                    profile[::2].copy(),profile[1::2].copy() # profile rho and concentration
            if self.is_gal:
                self.data['sig'], self.data['sigbulge'], self.data['mbulge'] =\
                        temp_gal[::3].copy(), temp_gal[1::3].copy(), temp_gal[2::3].copy()
                self.data['g_nbin'] = temp[20]
                self.data['g_rr'] = temp[21].reshape(ntot,100)
                self.data['g_rho']= temp[22].reshape(ntot,100)

        except:
            print("Something wrong")

    def load_hm_old(self, nout=None, base=None, info=None):
        if nout is None:
            nout = self.nout
        if base is None:
            base = self.base
        snout = str(self.nout).zfill(3)
        if self.is_gal:
            fn = base + self.gal_find_dir + 'gal/tree_bricks' + snout
        else:
            fn = base + self.dm_find_dir + 'DM/tree_bricks' + snout
        try:
            f = open(fn, "rb")
            import numpy as np
            from load.utils import read_fortran, skip_fortran
            self.nbodies = read_fortran(f, np.dtype('i4'), 1)
            self.massp = read_fortran(f, np.dtype('f4'), 1)
            self.aexp = read_fortran(f, np.dtype('f4'), 1)
            self.omegat = read_fortran(f, np.dtype('f4'), 1)
            self.age = read_fortran(f, np.dtype('f4'), 1)
            self.halnum, self.subnum = read_fortran(f, np.dtype('i4'), 2)

            dtype_halo = [('np', '<i4'), ('id', '<i4'), ('level', '<i4'),
                          ('host', '<i4'), ('sub', '<i4'), ('nsub', '<i4'),
                        ('nextsub', '<i4'),
                        ('m', '<f4'), ('mvir', '<f4'),
                        ('r', '<f4'), ('rvir', '<f4'), 
                        ('tvir', '<f4'), ('cvel', '<f4'), 
                        ('x', '<f4'), ('y', '<f4'), ('z', '<f4'),
                        ('vx', '<f4'), ('vy', '<f4'), ('vz', '<f4'),
                        ('ax', '<f4'), ('ay', '<f4'), ('az', '<f4'),
                        ('sp', '<f4'), ('idx', '<i4'),
                        ('p_rho', '<f4'),('p_c', '<f4')]
            if self.is_gal:
                dtype_halo += [('sig', '<f4'), ('sigbulge', '<f4'),
                               ('mbulge', '<f4'), ('hosthalo', '<i4')]

            tothal = self.halnum + self.subnum
            self.data = np.recarray(tothal, dtype=dtype_halo)
#            print("retunr id?", self.return_id)
            if self.return_id:
                self.idlists=[]
                self.hal_idlists=[]
          
#            print("idlists", self.idlists)          
            for i in range(tothal):
                nph = read_fortran(f, np.dtype('i4'), 1)
                self.data['np'][i] = nph

                tmp = read_fortran(f, np.dtype('i4'), nph) # id list. 
                hnu = read_fortran(f, np.dtype('i4'), 1)[0]
                self.data['id'][i] = hnu
#                print(self.return_id_list, "return_id_list")
                if self.return_id:
                    if self.return_id_list is not None:
                        if hnu in self.return_id_list:
                            self.idlists.append(tmp)
                            self.hal_idlists.append(hnu)
                    else:
                        #pass
                  # By default save all id lists. 
                        self.idlists.append(tmp)
            
                read_fortran(f, np.dtype('i4'), 1) #timestep
                self.data['level'][i], self.data['host'][i], \
                self.data['sub'][i], self.data['nsub'][i], \
                self.data['nextsub'][i] = read_fortran(f, np.dtype('i4'), 5)
                self.data['m'][i] = read_fortran(f, np.dtype('f4'), 1)
                self.data['x'][i], self.data['y'][i], self.data['z'][i] \
                    = read_fortran(f, np.dtype('f4'), 3)
                self.data['vx'][i], self.data['vy'][i], self.data['vz'][i] \
                    = read_fortran(f, np.dtype('f4'), 3)
                self.data['ax'][i], self.data['ay'][i], self.data['az'][i] \
                    = read_fortran(f, np.dtype('f4'), 3)                
                self.data['r'][i] = read_fortran(f, np.dtype('f4'), 4)[0]
                read_fortran(f, np.dtype('f4'), 3)#energies
                self.data['sp'][i] = read_fortran(f, np.dtype('f4'), 1)
                if self.is_gal:
                    self.data['sig'][i], self.data['sigbulge'][i], self.data['mbulge'][i] \
                        = read_fortran(f, np.dtype('f4'), 3)
                    #skip_fortran(f)
                    
                self.data['rvir'][i], self.data['mvir'][i],self.data['tvir'][i],self.data['cvel'][i] \
                        = read_fortran(f, np.dtype('f4'), 4)
                self.data['p_rho'][i], self.data['p_c'][i] = read_fortran(f, np.dtype('f4'), 2)# profile rho and concentration
                if self.is_gal:
                    # stellar surface profile
                    skip_fortran(f) # nbin
                    skip_fortran(f) # rr
                    skip_fortran(f) # rho
    
            f.close()
    
            if self.return_id_list is None:
                self.hal_idlists = self.data['id'] 
#            self.refactor_hm()
        except IOError:
            print("Couldn't find {}".format(fn))

    def refactor_hm(self):
        """
        refactor HaloMaker halo into Rockstar format. (mainly names are modified.)
        """
        import numpy as np
        data = self.data
        """
        self.aexp = data.aexp
        self.age = data.age
        self.npart_all = data.nbodies
        self.massp = data.massp
        self.nhalo = data.halnum[0]
        self.nsub = data.subnum[0]
        """
        # refactor the recarray.
        names = ["id", "np", "m", "mvir", "r", "rvir",
                 "x", "y", "z", "vx", "vy", "vz",
                 "level", "host", "sub", "nsub", "nextsub",
                 "spinx", "spiny", "spinz", "angx", "angy", "angz"]

        dtypes = "int64, int64, float64, float64, float64, float64, \
                  float64, float64, float64, float64, float64, float64,\
                  int8, int64, int64, int32, int64, \
                  float64, float64, float64, float64, float64, float64"

        self.data = np.rec.fromarrays([data["HNU"][0], data["NP"][0],
                                 data["M"][0], data["MVIR"][0],
                                 data["R"][0], data["RVIR"][0],
                                 data["P"][0][0], data["P"][0][1], data["P"][0][2],
                                 data["V"][0][0], data["V"][0][1], data["V"][0][2],
                                 data["HHOST"][0][0], data["HHOST"][0][1],
                                 data["HHOST"][0][2], data["HHOST"][0][3], data["HHOST"][0][4],
                                 data["SP"][0][0], data["SP"][0][1], data["SP"][0][2],
                                 data["ANG"][0][0], data["ANG"][0][1], data["ANG"][0][2]],
                                 dtype = dtypes)
        self.data.dtype.names = names

    def normalize(self):
        """
        normalize qunatities in comoving scale.
        Works only with internal info. So set Halo.info first.
        """
        if self.halofinder is 'Rockstar':
            self.normalize_rs()
        elif self.halofinder is 'HaloMaker':
            self.normalize_hm()
            # modify units
            # All in comoving scale. (But mass.. hmm..)
            # Length, distance in kpc/h
            # mass in Msun/h

    def normalize_rs(self):
        """
        Normalizes positions and lenghts in code unit ([0,1]),
        and mass in solar mass (physical).
        """
        self.data['x'] = self.data['x'] / self._boxsize
        self.data['y'] = self.data['y'] / self._boxsize
        self.data['z'] = self.data['z'] / self._boxsize
        self.data['rvir'] = self.data['rvir'] / self._boxsize / 1000 # in code unit.
#        self.data['m'] = self.data['m'] * 1e11 # / info.h
#        self.data['mvir'] = self.data['mvir'] * 1e11 # / info.h
        self.unit.update({"Mass":"Msun (physical)",
                          "Length":"code unit (0 - 1), (comoving)",
                            "velocity":"km/s physical"})

    def normalize_hm(self):
        """
        Normalizes positions and lenghts in code unit ([0,1]),
        and mass in solar mass (physical).
        """
        info = self.info
        self.data['x'] = self.data['x'] / info.pboxsize + 0.5 #* info.cboxsize * 1e3
        self.data['y'] = self.data['y'] / info.pboxsize + 0.5 #* info.cboxsize * 1e3
        self.data['z'] = self.data['z'] / info.pboxsize + 0.5 #* info.cboxsize * 1e3
        self.data['r'] = self.data['r'] / info.pboxsize  #* info.cboxsize * 1e3
        self.data['rvir'] = self.data['rvir'] / info.pboxsize #* info.cboxsize * 1e3
        self.data['m'] = self.data['m'] * 1e11 # / info.h
        self.data['mvir'] = self.data['mvir'] * 1e11 # / info.h
        self.unit.update({"Mass":"Msun (physical)",
                          "Length":"code unit (0 - 1), (comoving)",
                            "velocity":"km/s physical"})

    def load_rs_ascii(self, nout=None, base=None, info=None, pickle=True):
        print("Loading Rockstar halos from ascii files")
        if pickle:
            print("Halos will be automatically pickled.")
        from tree import rshalo
        if nout is None:
            nout = self.nout
        if base is None:
            base = self.base
        fname = base + "rockstar_halos/" + "halos_" + str(nout) + "."
        self.data = rshalo.read_halo_all(fname, sort=True)
        if pickle:
            self.pickle_halo(fname=fname)

    def load_rs(self, nout=None, base=None, info=None):
        import pickle
        import numpy as np
        if nout is None:
            nout = self.nout
        if base is None:
            base = self.base
        fname = base + 'rhalo/halos_py/halos_' + str(nout).zfill(3) + '.pickle'
        with open(fname, 'rb') as f:
            data = pickle.load(f)
            self.data = np.rec.fromrecords(np.hstack(data), dtype=data.dtype)
            self._boxsize = 199.632011
            # Temporarly, to wrap the data in recarray.
        self.normalize_rs()

    def pickle_halo(self, fname=None):
        import pickle
        with open(fname, mode='rw'):
            pickle.dump(self, fname)

    def derive_from(self, old_halo, ind=None):
        """
        derive a subset from another halo instance.
        1) all meta data will be overwritten.
        2) ind specifies the part of .data to be inherited.

        Last modified
        """
        if ind is None:
            ind = range(len(old_halo.data))
        # How can I force a integer to be a single element list?
        for key, val in old_halo.__dict__.items():
            if key == "data" or key == "nhalo":
                continue
            self.__dict__.update({key:val})
        self.data = old_halo.data[ind]
        try:        
            self.nhalo = len(ind)
        except:
            self.nhalo = len([ind])


