# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 11:00:31 2015

halo / galaxy calss including basic data load functionality.

@author: hoseung
"""

from tree import rd_hal as rd_halo
import numpy as np
from utils.io_module import read_fortran, skip_fortran
from load.info import Info
import struct
class HaloMeta():
    """
    HaloMeta class.


    all except halofinder output.


    Attributes
    ----------

    Methods
    -------
    set_info(self, info)
    set_nout(self, nout)

    set_halofinder(self, halofinder)



    """
    def __init__(self, nout=None, base='./', info=None, halofinder='HM',
                 load=True, is_gal=False, return_id=False, outdir=None, fn=None,
                 verbose=False, double=False, pure=False):
        """

        Parameters
        ----------
        nout : int
            snapshot number
        base : str
            working directory
        info : info object

        halofinder : {"RS", "HM"}
            Full names also work. case insensitive.

        return_id : logical or sequence
            If True, load (as hcat.hal_idlists) constituent particle id of all halo.
            If halo ids are given, load (as hcat.hal_idlists) particle ids of the given halos.

        Examples
        --------
        It is better to specify nout, base, halofinder from the beginning.
        All three are necessary to load a halo output

        >>> h = tree.halomodule.Halo(nout=132, halofinder="RS", base='~/data/AGN2/')

        given nout and base, info is auto-loaded if not explicitely given.

        """
        self.fn = fn
        self.nout = nout
        self.verbose = verbose
        self.double = double
        self.pure = pure
        self.base = base
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
        self.convert=True
        if return_id is False:
            self.return_id = False
        else:
            if hasattr(return_id, "__len__"):
                # String also has __len__, but let's just ignore such cases.
                self._return_id_list = return_id
            else:
                self._return_id_list = None # None = load all halo's ids.
            self.return_id = True

        if outdir is None:
            if is_gal:
                self.gal_find_dir = 'GalaxyMaker/'
            else:
                self.dm_find_dir= 'halo/'
        else:
            if is_gal:
                self.gal_find_dir = outdir
            else:
                self.dm_find_dir = outdir

        try:
            self.set_nout(nout)
        except:
            pass

        if load:
            self.load()


    def set_info(self, info):
        if info is None:
            try:
                self._load_info()
            except:
                print("[Halo.set_info] Couldn't load info file.")
                raise # automatically raise the most recent error
        else:
            self.info = info


    def _load_info(self):
        #if self.verbose:
        if True:
            print("[Halo.load_info] loading info")
            print("[Halo.load_info] nout = {}, base ={}".format(self.nout, self.base))
        self.info = Info(nout=self.nout, base=self.base, load=True)
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


    def set_nout(self, nout):
        self.nout = nout
        self._set_aexp(nout)


class Halo(HaloMeta):
    """
    Class to hold Halo finder catalog.

    HM :
        r = distance of the furthest particle from the center of the halo.
        rvir = radius where the density is 200 times the _______?

    Notes
    -----
    Rockstar halo id starts from 0. CT id too.
    For the details of Halofinders refer Aubert 2004, Behroozi 2013.


    """
    def __init__(self, convert=True, **kwargs):
        self.convert = convert
        super(Halo, self).__init__(**kwargs)


    def _check_params(self):
        assert (self.base is not None), "No working directory given : {}".format(self.base)
        assert (self.nout is not None), "No nout given : {}".format(self.nout)
        #assert (self.base is not None), "No  : {}".format(self.base)

    def set_data(self, data):
        if data is None:
            self.load()
        else:
            self.data = data

    def load(self, nout=None, base=None, info=None, pure=None, double=None):
        """
        There are nout, base keywords.
        But self.nout and self.base are already available.
        Determine the priority among them.
        """
        if double is None:
            double = self.double
        if pure is None:
            pure = self.pure

        if self.fn is None:
            self._check_params()
        if self.halofinder is 'Rockstar':
            self.load_rs()
        elif self.halofinder is 'HaloMaker':
            if self.fn is None:
                if nout is None:
                    nout = self.nout
                if base is None:
                    base = self.base
                try:
                    self.data = pickle.load(open(base + self.gal_find_dir + "gal_pickle/gcat_{}.pickle".format(nout), "rb"))
                    return
                except:
                    pass

                snout = str(self.nout).zfill(3)
                if self.is_gal:
                    self.fn = base + self.gal_find_dir + 'gal/tree_bricks' + snout
                else:
                    self.fn = base + self.dm_find_dir + 'DM/tree_bricks' + snout
            print(self.fn)
            if self.verbose:
                print("Loading file:", self.fn)
            print("Loading...", double, pure)
            self.load_hm(self.fn, double=double, pure=pure)
            if self.info is None:
                info = Info(base = self.base, nout = self.nout, load=True)
                self.set_info(info)
        if self.convert:
            self.normalize()
        else:
            print("Not converting unit!")

    def load_hm(self, fn, double=None, pure=None):
        if double == None:
            double = self.double
        if pure == None:
            pure = self.pure

        if double:
            dtype_float = "<f8"
        else:
            dtype_float = "<f4"

        dtype_halo = [('np', '<i4'), ('id', '<i4'), ('level', '<i4'),
                      ('host', '<i4'), ('sub', '<i4'), ('nsub', '<i4'),
                      ('nextsub', '<i4'),
                      ('m', dtype_float), ('mvir', dtype_float),
                      ('r', dtype_float), ('rvir', dtype_float),
                      ('tvir', dtype_float), ('cvel', dtype_float),
                      ('x', dtype_float), ('y', dtype_float), ('z', dtype_float),
                      ('vx', dtype_float), ('vy', dtype_float), ('vz', dtype_float),
                      ('ax', dtype_float), ('ay', dtype_float), ('az', dtype_float),
                      ('sp', dtype_float), ('idx', '<i4'),
                      ('p_rho', dtype_float),('p_c', dtype_float),
                      ('energy', '<f8', (3,)), ('abc', '<f8', (3,))]

        if self.is_gal:
            dtype_halo += [('sig', dtype_float), ('sigbulge', dtype_float),
                           ('mbulge', dtype_float), ('hosthalo', '<i4'),
                           ('g_nbin', '<i4'), ('g_rr', dtype_float, (100,)),
                           ('g_rho', dtype_float, (100,))]


        f = open(fn, "rb")
        if pure:
            brick_data = f.read()
            offset, halnum, subnum = load_header(brick_data, double=double)
            self.data = np.zeros(halnum+subnum, dtype=dtype_halo)
            for i in range(halnum+subnum):
                offset = load_a_halo(brick_data, offset, self.data[i], is_gal=self.is_gal, double=double)
            f.close()
        else:
            self.nbodies = read_fortran(f, np.dtype('i4'), 1)[0]
            f.close()
            #self.nbodies = rd_halo.read_nbodies(fn.encode())
            temp = rd_halo.read_file(fn.encode(), self.nbodies, int(self.is_gal))# as a byte str.

            allID, __, self.halnum, self.subnum,\
                self.massp, self.aexp, self.omegat, self.age = temp[0:8]
            ntot = self.halnum + self.subnum
            self.data = np.recarray(ntot, dtype=dtype_halo)
            self.data['np'], self.data['id'],\
            levels, ang, energy, \
            self.data['m'],\
            radius, pos,\
            self.data['sp'], vel = temp[8:18]
            vir, profile = temp[18:20]
            if self.is_gal:
                temp_gal = temp[20]

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

#           copy so that memory is continuous. (Not tested!)
            self.data['rvir'],self.data['mvir'], \
                    self.data['tvir'],self.data['cvel'] = vir[::4].copy(),\
                vir[1::4].copy(),vir[2::4].copy(),vir[3::4].copy()
            self.data['p_rho'], self.data['p_c'] =\
                    profile[::2].copy(),profile[1::2].copy() # profile rho and concentration
            if self.is_gal:
                self.data['sig'], self.data['sigbulge'], self.data['mbulge'] =\
                        temp_gal[::3].copy(), temp_gal[1::3].copy(), temp_gal[2::3].copy()
                self.data['g_nbin'] = temp[21]
                self.data['g_rr'] = temp[22].reshape(ntot,100)
                self.data['g_rho']= temp[23].reshape(ntot,100)

        if self.return_id:
            self.idlists=[]
            self.hal_idlists=[]
            iskip=0
            for hid, hnp in zip(self.data["id"],self.data["np"]):
                if self._return_id_list is not None:
                    if hid in self._return_id_list:
                        self.idlists.append(allID[iskip:iskip+hnp])
                        self.hal_idlists.append(hid)
                else:
                    # for every halo.
                    self.idlists.append(allID[iskip:iskip+hnp])
                    self.hal_idlists.append(hid)
                iskip += hnp


    def refactor_hm(self):
        """
        refactor HaloMaker halo into Rockstar format. (mostly name modifications.)
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
                      data["HHOST"][0][2], data["HHOST"][0][3],
                      data["HHOST"][0][4],
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
        1) all meta data is overwritten.
        2) ind specifies the part of .data to be copied.

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





def load_header(brick_data, double=False):
    offset = 4
    if double:
        nbytes = 8
        dtype_float="d"
    else:
        nbytes = 4
        dtype_float="f"

    nbodies = struct.unpack("i", brick_data[4:8])[0]
    offset += 12
    massp = struct.unpack(dtype_float, brick_data[offset:offset+nbytes])[0]
    offset += 8 + nbytes
    aexp = struct.unpack(dtype_float, brick_data[offset:offset+nbytes])[0]
    offset += 8 + nbytes
    omegat = struct.unpack(dtype_float, brick_data[offset:offset+nbytes])[0]
    offset += 8 + nbytes
    age = struct.unpack(dtype_float, brick_data[offset:offset+nbytes])[0]
    offset += 8 + nbytes
    halnum = struct.unpack("i", brick_data[offset:offset+4])[0]
    subnum = struct.unpack("i", brick_data[offset+4:offset+8])[0]
    return offset+16, halnum, subnum




def load_a_halo(brick_data, offset, dd, is_gal=True, double=False):
    if double:
        nbytes = 8
        dtype_float="d"
    else:
        nbytes = 4
        dtype_float="f"

    npart = struct.unpack("i", brick_data[offset:offset+4])[0]
    dd["np"]=npart
    offset += 12  # 12 = 4 + 8
    ids = struct.unpack_from("<{}i".format(npart), brick_data[offset:offset+4*npart])
    offset += 4*npart + 8
    dd["id"] = struct.unpack("i", brick_data[offset:offset+4])[0]
    #offset += 12
    #dd["nstep"] = struct.unpack("i", brick_data[offset:offset+4])[0]
    offset += 24
    dd["level"],dd["host"],dd["sub"],dd["nsub"],dd["nextsub"]\
    = struct.unpack_from("<5i", brick_data[offset:offset+20])
    offset += 28
    dd["m"] = struct.unpack(dtype_float, brick_data[offset:offset+nbytes])[0]
    offset += 8 + nbytes
    dd["x"],dd["y"],dd["z"] = struct.unpack_from("<3"+dtype_float, brick_data[offset:offset+3*nbytes])
    offset += 8 + 3*nbytes
    dd["vx"],dd["vy"],dd["vz"] = struct.unpack_from("<3"+dtype_float, brick_data[offset:offset+3*nbytes])
    offset += 8 + 3*nbytes
    dd["ax"],dd["ay"],dd["az"] = struct.unpack_from("<3"+dtype_float, brick_data[offset:offset+3*nbytes])
    offset += 8 + 3*nbytes
    radius= struct.unpack_from("<4"+dtype_float, brick_data[offset:offset+4*nbytes])
    dd["r"],dd["abc"] = radius[0], radius[1:]
    offset += 8 + 4*nbytes
    dd["energy"] = struct.unpack_from("<3"+dtype_float, brick_data[offset:offset+3*nbytes])
    offset += 8 + 3*nbytes
    dd["sp"] = struct.unpack(dtype_float, brick_data[offset:offset+nbytes])[0]
    offset += 8 + nbytes
    if is_gal:
        dd["sig"], dd["sigbulge"], dd["mbulge"]\
        = struct.unpack_from("<3"+dtype_float, brick_data[offset:offset+3*nbytes])
        offset += 8+ 3*nbytes
    dd["mvir"],dd["rvir"],dd["tvir"],dd["cvel"]\
    = struct.unpack_from("<4"+dtype_float, brick_data[offset:offset+4*nbytes])
    offset += 8+4*nbytes
    dd["p_rho"],dd["p_c"] = struct.unpack_from("<2"+dtype_float, brick_data[offset:offset+2*nbytes])
    offset += 8+2*nbytes
    if is_gal:
        g_nbin = struct.unpack("i", brick_data[offset:offset+4])[0]
        dd["g_nbin"]=g_nbin
        offset += 12
        dd["g_rr"] = struct.unpack_from("<{}".format(g_nbin)+dtype_float, brick_data[offset:offset+g_nbin*nbytes])
        offset += 8 + g_nbin*nbytes
        dd["g_rho"] = struct.unpack_from("<{}".format(g_nbin)+dtype_float, brick_data[offset:offset+g_nbin*nbytes])
        offset += 8 + g_nbin*nbytes

    return offset
