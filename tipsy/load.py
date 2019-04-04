import numpy as np
from collections import namedtuple
from .config import Changa_config
from load.dtypes import add_dtypes

ccfg = Changa_config()

class TipsyHeader():
    def __init__(self, f=None, swap=True):
        self.aexp = 0
        self.zred = 0
        self.nbodies = 0
        self.ndim = 0
        self.nsph = 0
        self.ndark = 0
        self.nstar = 0
        self._setup=False
        if not f is None:
            self.read_header(f, swap)

    def __set_from_dict(self, dic):
        keys=dic.keys()
        try:
            self.aexp = dic["a"]
        except:
            self.aexp = dic["aexp"]

    def read_header(self, f, swap=False):
        """
        To do
        1. Make byteswap optional. - done

        """
        if not self._setup:
            if swap:
                self.aexp = np.fromfile(f, dtype=np.float64, count=1).byteswap()[0]
                self.nbodies = np.fromfile(f, dtype=np.int32, count=1).byteswap()[0]
                self.ndim = np.fromfile(f, dtype=np.int32, count=1).byteswap()[0]
                self.nsph = np.fromfile(f, dtype=np.int32, count=1).byteswap()[0]
                self.ndark = np.fromfile(f, dtype=np.int32, count=1).byteswap()[0]
                self.nstar = np.fromfile(f, dtype=np.int32, count=1).byteswap()[0]
                # 4 Bytes of padding to make the data 8-Byte aligned.
                np.fromfile(f, dtype=np.int32, count=1).byteswap()[0]
            else:
                self.aexp = np.fromfile(f, dtype=np.float64, count=1)[0]
                self.nbodies = np.fromfile(f, dtype=np.int32, count=1)[0]
                self.ndim = np.fromfile(f, dtype=np.int32, count=1)[0]
                self.nsph = np.fromfile(f, dtype=np.int32, count=1)[0]
                self.ndark = np.fromfile(f, dtype=np.int32, count=1)[0]
                self.nstar = np.fromfile(f, dtype=np.int32, count=1)[0]
                # 4 Bytes of padding to make the data 8-Byte aligned.
                np.fromfile(f, dtype=np.int32, count=1)[0]

            self.zred = 1./self.aexp -1
            self._setup = True
            return f.tell()

def convert_num_if_possible(s):
    """
    convert a string of integer or float into int or float value.
    It will be useful to convert  "yes"/"no" into True/False.
    """
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except:
            return s

def cal_units(param, header):
    from collections import namedtuple
    """
    (mass, time, density, velocity, ...) in (g, sec, g/cc and so on) converts
    a code unit into a more familiar unit.

    UnitA_in_UnitB converts quantity A in code unit into unit B.

    NOTE
    ----
    If I want to get rid of entailing units, the target unit system must be
    very general (for example everything in cgs like info.unit_d / unit_l of
    RAMSES analysis)
    What Pynbody currently does is :
    dunit converts the code unit into kpc,
    munit converts the code unit into Msun, and so on...
    It makes sense when we study galaxies. But the set of (kpc, Msun, km/s)
    is somewhat arbitrary.
    For example, both Msun/kpc^3 and g/cc are very often used.

    Think about it.

    """
    aexp = header.aexp
    G = 6.67408e-11 #m3 kg-1 s-2
    G_cgs = 6.67408e-11 * 1e6 * 1e-3

    Msun_in_g = 1.989e33

    kpc_in_cm = 3.08567758e21
    kpc_in_km = kpc_in_cm *1e-5

    yr_in_s = 31556952 # 365.2425 days
    Gyr_in_sec = 1e9*yr_in_s

    cosmo = bool(param['bComove'])

    mass_in_msun = np.float64(param["dMsolUnit"])
    # assume boxsize == 1.
    length_in_kpc = np.float64(param["dKpcUnit"])
    density_in_Msol_kpc2 = mass_in_msun/length_in_kpc**3 # Msol kpc^-3

    Msolkpc2_in_cgs = Msun_in_g/kpc_in_cm**3
    density_in_cgs = density_in_Msol_kpc2 * Msolkpc2_in_cgs

    time_in_sec = 1./np.sqrt(density_in_cgs*G_cgs)
    vel_in_kms = length_in_kpc*kpc_in_km/time_in_sec # in km/s

    time_in_gyr = time_in_sec * Gyr_in_sec**-1 # in Gyr

    if cosmo:
        hub = np.float64(param["dHubble0"])
        hubble = round(hub * 10 * vel_in_kms / length_in_kpc, 4) 
        # 1/100 * kms/Mpc = 10 * kms/kpc
        # 67.74, not 67.7403103305142


    Info = namedtuple("info", ["m_in_msun",
                               "t_in_gyr",
                               "l_in_kpc",
                               "d_in_gcc",
                               "v_in_kms",
                               "h0",
                               "cboxsize",
                               "Om", "Ol", "aexp", "zred"])
    info = Info(mass_in_msun,
    time_in_gyr,
    length_in_kpc * aexp,
    density_in_cgs,
    vel_in_kms,
    hubble,
    length_in_kpc * hubble * 1e-3,
    param["dOmega0"],
    param["dLambda"],
    aexp,
    1/aexp -1
    )
    return info


class TipsySim():
    def __init__(self,
                fn=None,
                base=None,
                nout=None,
                load=False,
                ignore_param=False,
                IC=False,
                BYTE_SWAP=True):

        """
        parameters
        ----------
        base =
        nout =
        fn   =
        load = False
        ignore_param = False
            Do not load .param. No unit conversion constant will be available.
        BYTE_SWAP = True

        """

        self._BYTE_SWAP=BYTE_SWAP
        self.param = {}
        self.header = TipsyHeader()
        """
        Incomplete. Haven't decided what to do when fn == None.
        """
        self._gas_aux =[]
        self._star_aux =[]
        self._dm_aux =[]

        if base is not None:
            self.base = base
        if nout is not None:
            self.nout = nout

        self.IC = IC
        if self.IC:
            ignore_param = True

        self.set_fn(fn)
        self._f = open(self._fn, "rb")
        self._size_header = self.header.read_header(self._f, self._BYTE_SWAP)
        if not ignore_param:
            self.read_param()

        #print("BYTESWAP", self._BYTE_SWAP)
        if load:
            self.load_data()

        if not ignore_param:
            self.info = cal_units(param=self.param, header=self.header)

    def set_fn(self, fn):
        """

        """
        if self.IC:
            self._fn = fn
        else:
            try:
                self._fn = "{}.{:06d}".format(self.base, self.nout)
            except:
                self.nout = int(fn.split(".")[-1])
                self.base = fn[:-6]
                self._fn = fn


    def set_aux_list(self, aux_list):
        """
        check if given aux file is available

        Note
        ----
        dtype_aux_gas.append((aux,dt_aux.get(name)[0], 1, "", 0))
        -> 1, "", 0 is needed as input to load.dtype.add_dtype()
        """
        dtype_aux_star = []
        dtype_aux_gas = []
        dtype_aux_dm = []
        dt_aux = ccfg.dtype_aux.fields
        for aux in aux_list:
            if aux in ccfg.gas_properties:
                self._gas_aux.append(aux)
                dtype_aux_gas.append((aux,dt_aux.get(aux)[0], 1, "", 0))
            if aux in ccfg.star_properties:
                self._star_aux.append(aux)
                dtype_aux_star.append((aux,dt_aux.get(aux)[0], 1, "", 0))
            if aux in ccfg.dm_properties:
                self._dm_aux.append(aux)
                dtype_aux_dm.append((aux,dt_aux.get(aux)[0], 1, "", 0))

        self._dtype_star = add_dtypes(ccfg.dtype_star, dtype_aux_star)
        self._dtype_gas  = add_dtypes(ccfg.dtype_gas, dtype_aux_gas)
        self._dtype_dm   = add_dtypes(ccfg.dtype_dm, dtype_aux_dm)


    def load_data(self, gas=True, star=True, dm=True):
        """
        Todo : skip unwanted populations.

        The resulting array is read-only if a string, which is immutable, is in the buffer.
        """
        f = self._f

        dtype_gas  = ccfg.dtype_gas
        dtype_dm   = ccfg.dtype_dm
        dtype_star = ccfg.dtype_star

        if gas:
            f.seek(self._size_header)
            gas_tmp= np.frombuffer(f.read(ccfg.dtype_gas.itemsize * self.header.nsph),
                        dtype=ccfg.dtype_gas).copy()
            if len(self._gas_aux) > 0:
                self.gas = np.zeros(self.header.nsph, dtype=self._dtype_gas)
                # Cannot use dtype.descr because it's not defined for types with overlapping or out-of-order fields
                for (key, type) in gas_tmp.dtype.fields.items():
                    self.gas[key] = gas_tmp[key]
            else:
                self.gas = gas_tmp

        if dm:
            f.seek(self._size_header + ccfg.dtype_gas.itemsize * self.header.nsph)
            dm_tmp = np.frombuffer(f.read(ccfg.dtype_dm.itemsize * self.header.ndark),
                                dtype=ccfg.dtype_dm).copy()
            if len(self._dm_aux) > 0:
                #dtype_dm_new = append_dtype(ccfg.dtype_star, self.list_aux)
                self.dm = np.zeros(self.header.ndark, dtype=self._dtype_dm)
                for (key, type) in dm_tmp.dtype.fields.items():
                    self.dm[key] = dm_tmp[key]
            else:
                self.dm = dm_tmp

        if star:
            f.seek(self._size_header + \
                   ccfg.dtype_gas.itemsize * self.header.nsph +\
                   ccfg.dtype_dm.itemsize * self.header.ndark)
            star_tmp = np.frombuffer(f.read(ccfg.dtype_star.itemsize * self.header.nstar),
                                dtype=ccfg.dtype_star).copy()
            if len(self._star_aux) > 0:
                #dtype_star_new = append_dtype(ccfg.dtype_star, self.list_aux)
                self.star = np.zeros(self.header.nstar, dtype=self._dtype_star)
                for (key, type) in star_tmp.dtype.fields.items():
                    self.star[key] = star_tmp[key]
            else:
                self.star = star_tmp

        if self._BYTE_SWAP:
            if gas:  self.gas = self.gas.newbyteorder(">")
            if dm:   self.dm = self.dm.newbyteorder(">")
            if star: self.star = self.star.newbyteorder(">")

        f.close()
        # add ID?

    def read_param(self, fn_param=None):
        if fn_param is None:
            fn = self._fn
            try:
                # Default name
                fn_param = fn.replace(fn.split(".")[-1], "param")
                f = open(fn_param, "r")
            except:
                # Look for any .param file in the directory
                from glob import glob
                l = glob(fn.replace(fn.split("/")[-1], "/*.param"))
                if len(l) == 1:
                    f = open(l[0], "r")
                elif len(l) == 0:
                    raise IOError("No .param file is found in the directory")
                elif len(l) > 1:
                    raise IOError("There are multiple .param files in the directory")
        else:
            f = open(fn_param, "r")

        for line in f:
            try:
                if line[0] != "#":
                    s = line.split("#")[0].split()
                    self.param[s[0]] = convert_num_if_possible(" ".join(s[2:]))

            except IndexError as ValueError:
                pass
        f.close()
