import numpy as np

dtype_gas = np.dtype({  'm': (('<f4', 1), 0),
             'pos': (('<f4', (3,)), 4),
               'x': (('<f4', 1), 4),
               'y': (('<f4', 1), 8),
               'z': (('<f4', 1), 12),
             'vel': (('<f4', (3,)), 16),
              'vx': (('<f4', 1), 16),
              'vy': (('<f4', 1), 20),
              'vz': (('<f4', 1), 24),
             'rho': (('<f4', 1), 28),
            'temp': (('<f4', 1), 32),
             'eps': (('<f4', 1), 36),
           'metal': (('<f4', 1), 40),
             'phi': (('<f4', 1), 44)})

dtype_dm = np.dtype({  'm': (('<f4', 1), 0),
             'pos': (('<f4', (3,)), 4),
               'x': (('<f4', 1), 4),
               'y': (('<f4', 1), 8),
               'z': (('<f4', 1), 12),
             'vel': (('<f4', (3,)), 16),
              'vx': (('<f4', 1), 16),
              'vy': (('<f4', 1), 20),
              'vz': (('<f4', 1), 24),
             'eps': (('<f4', 1), 28),
             'phi': (('<f4', 1), 32)})

dtype_star = np.dtype({  'm': (('<f4', 1), 0),
             'pos': (('<f4', (3,)), 4),
               'x': (('<f4', 1), 4),
               'y': (('<f4', 1), 8),
               'z': (('<f4', 1), 12),
             'vel': (('<f4', (3,)), 16),
              'vx': (('<f4', 1), 16),
              'vy': (('<f4', 1), 20),
              'vz': (('<f4', 1), 24),
           'metal': (('<f4', 1), 28),
           'tform': (('<f4', 1), 32),
             'eps': (('<f4', 1), 36),
             'phi': (('<f4', 1), 40)})

class TipsyHeader():
    def __init__(self, f=None):
        self.aexp = 0
        self.zred = 0
        self.nbodies = 0
        self.ndim = 0
        self.nsph = 0
        self.ndark = 0
        self.nstar = 0
        self.pad = 0
        self._setup=False
        if not f is None:
            self.read_header(f)

    def read_header(self, f):
        """
        To do
        1. Make byteswap optional.

        """
        if not self._setup:
            self.aexp = np.fromfile(f, dtype=np.float64, count=1).byteswap()[0]
            self.nbodies = np.fromfile(f, dtype=np.int32, count=1).byteswap()[0]
            self.ndim = np.fromfile(f, dtype=np.int32, count=1).byteswap()[0]
            self.nsph = np.fromfile(f, dtype=np.int32, count=1).byteswap()[0]
            self.ndark = np.fromfile(f, dtype=np.int32, count=1).byteswap()[0]
            self.nstar = np.fromfile(f, dtype=np.int32, count=1).byteswap()[0]
            self.pad = np.fromfile(f, dtype=np.int32, count=1).byteswap()[0]
            self.zred = 1./self.aexp -1
            self._setup = True


def convert_num_if_possible(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except:
            return s

class Info():
    def __init__(self, param):
        """
        (mass, time, density, velocity, ...) in (g, sec, g/cc and so on) converts
        a code unit into a more familiar unit.

        UnitA_in_UnitB converts a quantites in unit A into unit B.

        NOTE
        ----
        If I want to get rid of entailing units, the target unit system must be
        very general (for example everything in cgs like info.unit_d / unit_l of
        RAMSES analysis)
        What Pynbody currently does is :
        dunit converts the code unit into kpc,
        munit converts the code unit into Msun, and so on...
        It makes sens when we study galaxies. But the set of (kpc, Msun, km/s)
        is somewhat arbitrary.

        Think about it.

        """
        G = 6.67408e-11 #m3 kg-1 s-2
        G_cgs = 6.67408e-11 * 1e6 * 1e-3

        Msun_in_g = 1.989e33

        kpc_in_cm = 3.08567758e21
        kpc_in_km = kpc_in_cm *1e-5

        yr_in_s = 31556952 # 365.2425 days
        Gyr_in_sec = 1e9*yr_in_s

        self.cosmo = bool(param['bComove'])

        mass_in_msun = np.float64(param["dMsolUnit"])
        # assume boxsize == 1.
        boxtokpc = np.float64(param["dKpcUnit"])
        density_in_Msol_kpc2 = mass_in_msun/length_in_kpc**3 # Msol kpc^-3

        Msolkpc2_in_cgs = Msun_in_g/kpc_in_cm**3
        density_in_cgs = density_in_Msol_kpc2 * Msolkpc2_in_cgs

        time_in_sec = 1./np.sqrt(unit_d_cgs*G_cgs)
        vel_in_kms = length_in_kpc*kpc_in_km/time_in_sec # in km/s

        time_in_gyr = time_in_sec * Gyr_in_sec**-1 # in Gyr

        if self.cosmo:
            hub = np.float64(param["dHubble0"])
            hubble = hub * 10 * vel_in_kms / length_in_kpc # 1/100 * kms/Mpc = 10 * kms/kpc


class TipsySim():
    def __init__(self, fn=None):
        self.param = {}
        self.header = TipsyHeader()
        """
        Incomplete. Haven't decided what to do when fn == None.
        """
        if fn is not None:
            self.set_fn(fn)
            self._f = open(self._fn, "rb")
            self.header.read_header(self._f)
            self.read_param()

    def set_fn(self, fn):
        """

        """
        self._fn = fn

    def load_data(self):
        f = self._f
        self.gas = np.frombuffer(f.read(np.dtype(dtype_gas).itemsize * self.header.nsph),
                    dtype=dtype_gas).newbyteorder(">")

        self.dm = np.frombuffer(f.read(np.dtype(dtype_dm).itemsize * self.header.ndark),
                            dtype=dtype_dm).newbyteorder(">")

        self.star = np.frombuffer(f.read(np.dtype(dtype_star).itemsize * self.header.nstar),
                            dtype=dtype_star).newbyteorder(">")
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




def cal_units(param):
    from collections import namedtuple
    """
    (mass, time, density, velocity, ...) in (g, sec, g/cc and so on) converts
    a code unit into a more familiar unit.

    UnitA_in_UnitB converts a quantites in unit A into unit B.

    NOTE
    ----
    If I want to get rid of entailing units, the target unit system must be
    very general (for example everything in cgs like info.unit_d / unit_l of
    RAMSES analysis)
    What Pynbody currently does is :
    dunit converts the code unit into kpc,
    munit converts the code unit into Msun, and so on...
    It makes sens when we study galaxies. But the set of (kpc, Msun, km/s)
    is somewhat arbitrary.
    For example, both Msun/kpc^3 and g/cc are very often used.

    Think about it.

    """
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
        hubble = hub * 10 * vel_in_kms / length_in_kpc # 1/100 * kms/Mpc = 10 * kms/kpc


    Info = namedtuple("info", ["m_in_msun",
                               "t_in_gyr",
                               "l_in_kpc",
                               "d_in_gcc",
                               "v_in_kms",
                               "hubble"])
    info = Info(mass_in_msun, time_in_gyr, length_in_kpc, density_in_cgs, vel_in_kms, hubble)
    return info
