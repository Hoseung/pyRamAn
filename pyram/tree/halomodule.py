# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 11:00:31 2015

halo / galaxy calss including basic data load functionality.

@author: hoseung
"""

import numpy as np
from numpy.core.records import fromarrays
import struct

from ..load.info import Info
from ..load.dtypes import get_halo_dtype, add_dtypes
from ..utils.io_module import read_fortran
from .readhtm import readhtm as readh

AdaptaHOP_params=dict(
    gdir='GalaxyMaker/gal/',
    gfn = lambda nout : f'tree_bricks{nout:03d}',
    ddir='halo/DM/',
    dfn = lambda nout : f'tree_bricks{nout:03d}',
)

HOP_params=dict(
    gdir='GalaxyMaker/gal_HOP/',
    gfn = lambda nout : f'tree_brick_{nout:03d}',
    ddir='halo/DM_HOP/',
    dfn = lambda nout : f'tree_brick_{nout:03d}',
)



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
                 verbose=False, double=False, pure_python=False, read_mbp=False,
                 HOP=False, add_fields=None):
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

        add_fields : None
            New fields is added to the halo catalog by using the "add_dtypes" function.
            The add_fields must be in the following format:
            (name, type, shape, existing field, offset wrt the existing field)
            For example,
            ("pos", dtype_float, (3,), "x", 0),
            ("vel", dtype_float, (3,), "vx", 0),
            ("lvec", dtype_float, (3,), "ax", 0)

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
        self.pure_python = pure_python
        self.base = base
        self.HOP = HOP
        self.hal_param = HOP_params if self.HOP else AdaptaHOP_params
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
        self.read_mbp = read_mbp
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
                self.gal_find_dir = self.hal_param['gdir']
            else:
                if self.read_mbp:
                    self.dm_find_dir = 'halo/DM_mbp'
                else:
                    self.dm_find_dir= self.hal_param['ddir']
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
            self.load(add_fields=add_fields)


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
        from ..utils import util
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

    def load(self, nout=None, base=None, info=None, pure_python=None, double=None,
             add_fields=None):
        """
        There are nout, base keywords.
        But self.nout and self.base are already available.
        Determine the priority among them.
        """
        if double is None:
            double = self.double
        if pure_python is None:
            pure_python = self.pure_python

        if self.fn is None:
            self._check_params()
        if self.halofinder == 'Rockstar':
            self.load_rs()
        elif self.halofinder == 'HaloMaker':
            # Get the filename
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
                    self.fn = base + self.gal_find_dir# + 'tree_bricks' + snout
                else:
                    self.fn = base + self.dm_find_dir# + self.hal_param['gfn'](self.nout)

            if self.verbose:
                print(self.fn)
                print("Loading file:", self.fn)
                print("double :", double, "Pure Python:", pure_python, "read_mbp",self.read_mbp)
            self.load_hm(self.fn, double=double, pure_python=pure_python, add_fields=add_fields)
            if self.info is None:
                info = Info(base = self.base, nout = self.nout, load=True)
                self.set_info(info)
        if self.convert:
            self.normalize()
        else:
            print("Not converting unit!")

    def load_hm(self, fn, double=None, pure_python=None, add_fields=None):
        if double == None:
            double = self.double
        if pure_python == None:
            pure_python = self.pure_python
        fname = fn+ self.hal_param['gfn'](self.nout)
        print("FNAME:", fname)
        f = open(fname, "rb")
        dtypes_halo = get_halo_dtype(is_gal=self.is_gal,
                                    double=double,
                                    read_mbp=self.read_mbp,
                                    new_fields=add_fields,
                                    auto_add_field=False)
        if pure_python:
            brick_data = f.read()
            offset, header_info = load_header(brick_data, double=double)
            nbodies, aexp, omegat, age, halnum, subnum = header_info

            #self.data = np.zeros(halnum+subnum,
            #                     dtype=dtypes_halo)
            for i in range(halnum+subnum):
                offset = load_a_halo(brick_data, offset, self.data[i],
                                     is_gal=self.is_gal, double=double)
            f.close()
        else:

            self.nbodies = read_fortran(f, np.dtype('i4'), 1)[0]
            f.seek(0)

            ###############################
            # Read header
            #with open(fname, "rb") as f:
            brick_data = f.read()
            
            offset, header_info = load_header(brick_data, double=double)
            self.nbodies, self.aexp, self.omegat, self.age, self.nhalo, self.nsub = header_info
            f.close()
            self.return_id=0
            readh.read_bricks(fname[:-3], self.is_gal, self.nout, self.nout+1, self.return_id, double)

            if(not double):
                self.data = fromarrays([*readh.integer_table.T, *readh.real_table.T], dtype=dtypes_halo)
                dtype_float = "<f4"
            else:
                self.data = fromarrays([*readh.integer_table.T, *readh.real_table_dp.T], dtype=dtypes_halo)
                dtype_float = "<f8"

            #array = HaloMaker.unit_conversion(array, snap)
            add_dtype = [("pos", dtype_float, (3,), "x", 0),
                        ("vel", dtype_float, (3,), "vx", 0),
                        ("lvec", dtype_float, (3,), "ax", 0)]
            newdt = add_dtypes(self.data.dtype, add_dtype)
            self.data.dtype = np.dtype(newdt)

            if(self.data.size==0):
                print("No tree_brick file found, or no halo found in %s" % fn)

        if self.return_id:
            """
            TODO
            To be compatible with AHF,
            self.idlists.append(dict(hid=hid,
                                     iddm = allID[iskip:iskip+hnp]))

            This will break a few things in analysis pipeline, though.
            """
            self.idlists=[]
            self.hal_idlists=[]
            iskip=0
            #######################
            All_ID = np.array(readh.part_ids)
            for hid, hnp in zip(self.data["id"],self.data["np"]):
                if self._return_id_list is not None:
                    if hid in self._return_id_list:
                        self.idlists.append(All_ID[iskip:iskip+hnp])
                        self.hal_idlists.append(hid)
                else:
                    # for every halo.
                    self.idlists.append(All_ID[iskip:iskip+hnp])
                    self.hal_idlists.append(hid)
                iskip += hnp

        #############################
        readh.close()


    def refactor_hm(self):
        """
        refactor HaloMaker halo into Rockstar format. (mostly name modifications.)
        """
        import numpy as np
        data = self.data
        # refactor the recarray.
        names = ["id", "np", "m", "mvir", "r", "rvir",
                 "x", "y", "z", "vx", "vy", "vz",
                 "level", "host", "sub", "nsub", "nextsub",
                 "spinx", "spiny", "spinz", "angx", "angy", "angz"]

        _dtypes = "int64, int64, float64, float64, float64, float64, \
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
                      dtype = _dtypes)
        self.data.dtype.names = names

    def normalize(self):
        """
        normalize qunatities in comoving scale.
        Works only with internal info. So set Halo.info first.
        """
        if self.halofinder == 'Rockstar':
            self.normalize_rs()
        elif self.halofinder == 'HaloMaker':
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
        self.unit.update({"Mass":"Msun (physical)",
                          "Length":"code unit (0 - 1), (comoving)",
                            "velocity":"km/s physical"})

    def normalize_hm(self):
        """
        Normalizes positions and lenghts in code unit ([0,1]),
        and mass in solar mass (physical).
        """
        info = self.info
        self.data['x'] = self.data['x'] / info.pboxsize + 0.5
        self.data['y'] = self.data['y'] / info.pboxsize + 0.5
        self.data['z'] = self.data['z'] / info.pboxsize + 0.5
        self.data['r'] = self.data['r'] / info.pboxsize
        self.data['rvir'] = self.data['rvir'] / info.pboxsize
        self.data['m'] = self.data['m'] * 1e11
        self.data['mvir'] = self.data['mvir'] * 1e11
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
        #if pickle:
        #    self.pickle_halo(fname=fname)

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

    def cal_purity(self):
        """
        Calculate purity by the mean mass of the DM particles.
        """
        mean_m = self.data["m"] / self.data["np"]
        M_hires = mean_m.min()
        try:
            self.data["purity"] = np.clip(1 - (mean_m / M_hires -1) / 8, 0,1)
        except:
            from numpy.lib import recfunctions as rfc
            self.data = rfc.append_fields(self.data, "purity",
                                          np.clip(1 - (mean_m / M_hires -1) / 8, 0,1),
                                           dtypes=float, usemask=False)

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
    return offset+16, (nbodies, aexp, omegat, age, halnum, subnum)




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
