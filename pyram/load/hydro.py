# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 21:32:55 2015

@author: hoseung
"""
import numpy as np
from .sim import Simbase
from ..utils.io_module import read_header
import struct
import gc
from .info import Info
from .pure_python import amr2cell_py
from .a2c import a2c

def generate_fname(nout,path="",ftype="",cpuid=1,ext=""):

    if len(path) > 0:
        if path[-1] != "/":
            path=path+"/"

    if nout == -1:
        from glob import glob
        filelist = sorted(glob(path+"output*"))
        number = filelist[-1].split("_")[-1]
    else:
        number = str(nout).zfill(5)

    infile = path+"output_"+number
    if len(ftype) > 0:
        infile += "/"+ftype+"_"+number
        if cpuid >= 0:
            infile += ".out"+str(cpuid).zfill(5)

    if len(ext) > 0:
        infile += ext

    return infile


class Dummy():
    def __init__(self):
        pass


class Hydro(Simbase):
    def __init__(self,
                 config,
                 nout=None,
                 info=None,
                 region=None,
                 ranges=None,
                 load=False,
                 #cosmo=True,
                 cpus=None,
                 cpu_fixed=None,
                 fields=None,
                 nvarh=6,
                 amr2cell_params={}):
        """
        Parameters
        ----------
        region : dict-like
            d
        ranges : array-like (3 by 2)
            region preceeds(?) ranges.
        nvarh : int
            Desired number of hydro variable to load.
            nvarh=6 loads rho, vel, temp, and metal.
            nvarh=12 load all chemical components (HAGN)

        Note
        ----
            It is important to distingush the number of hydro-variables to read and
            the number of hydrovaraibles in the simulation.

            Recent version of RAMSES generates "hydro_file_descriptor.txt" by default.
            This file contains nvarh and the type of each hydro-variable.
            With this file, nvarh input is not needed.

            self.amr is set inside the Simbase instance, which is evaluated before
            the Hydro.init(). Thus, if Hydro.init() "initializes" as self.amr = None,
            the already-set amr instance is gone. So don't do that.

        """
        super(Hydro, self).__init__()
        self.config=config
        #self.cosmo = cosmo
        if info is None:
            assert nout is not None, "either info or nout is required"
            print("[Hydro.__init__] Loading info")
            #info = Info(nout=nout, cosmo=cosmo)
            info = Info(base=base,
                        nout=nout, 
                        cosmo=self.config.cosmo, 
                        data_dir=data_dir)
        self.info = info
        self.nout = info.nout
        self.cpus = cpus
        self.cpu_fixed=cpu_fixed
        self.fields=fields
        try:
            self.ncpu = len(self.cpus)
        except:
            self.ncpu = 0

        try:
            self.base = self.info.base
        except:
            self.base = base

        snout = str(info.nout).zfill(5)
        # file name
        #self.out_dir = 'snapshots/'
        self.set_fnbase(self.base, self.info.data_dir, 'hydro')
        #self.data_dir = info.base + '/' + self.out_dir \
        #                         + 'output_' + snout \
        #                         + '/hydro_' + snout \
        #                         + '.out'
        self._get_basic_info()
        self.set_info(info)
        self.header.nvarh=nvarh
        if region is not None:
            ranges = region.ranges
        if ranges is not None:
            self.set_ranges(ranges=ranges)
        elif hasattr(info, "ranges"):
            if info.ranges is not None:
                self.ranges=info.ranges
                self.cpus = info.cpus
        else:
            self.set_ranges([[-9e9,9e9]] * 3)

        if load:
            print("amr2cell_params", amr2cell_params)
            self.amr2cell(**amr2cell_params)

    def _get_basic_info(self):
        try:
            f = open(self._fnbase + '00001', "rb")
        except:
            from glob import glob
            hydros = glob(self._fnbase + "*")
            f = open(hydros[0], "rb")

        self.header = Dummy()
        self._read_hydro_header(f)

    def read_hydro_dscr(self):
        def shorten(name):
            if 'density' in name:
                return name.replace('density', 'rho')
            elif 'velocity_' in name:
                return name.replace('velocity_', 'v')
            elif 'pressure' in name:
                return name.replace('pressure', 'temp')
            elif 'metallicity' in name:
                return name.replace('metallicity', 'metal')
            elif 'refinement_scalar' in name:
                return name.replace('refinement_scalar', 'ref')
            else:
                return name
            
        def to_dtype(dt):
            """
            I have seen 'd' only so far.
            Others are just a guess ;)
            """
            dt=dt.replace(" ","")
            if dt == 'd':
                return 'f8'
            if dt == 'f':
                return 'f4'
            if dt == 'i':
                return 'i4'
            if dt == 'l':
                return 'i8'

        dt = []
        with open(self.base + '/output_00001/hydro_file_descriptor.txt', 'r') as f:
            for l in f.readlines():
                if not l.startswith("#"):
                    ivar, varname, vartype = l.split(",")
                    # No whitespace allowed in a field name
                    varname = shorten(varname.replace(" ",""))
                    # detach '\n' and indicate 1-element
                    dt.append((varname, to_dtype(vartype.split()[0]), 1))
        
        from . import dtypes
        self.config.dtype_c = dtypes.add_dtypes(dtypes.hydro_dtypes_default, dt)

    def set_info(self, info):
        self.info = info

    def set_cell_dtype(self, dtype):
        new_ver = self.config.sim_type in ['fornax']
        if new_ver:
            self.read_hydro_dscr()
        else:
            pass
            # refer to original dtype in amr2cell() 
        
    def _read_hydro_header(self, f, verbose=False):
        # Global
        h1 = read_header(f,
                         np.dtype(
                             [('ncpu', 'i'),
                              ('nvarh', 'i4'),
                              ('ndim', 'i4'),
                              ('nlevelmax', 'i4'),
                              ('nboundary', 'i4'),
                              ('gamma', 'f8')]))

        if h1['nboundary'] == 0 and verbose:
            print(' Periodic boundary condition')

        # if assign
        self.header.ncpu = h1['ncpu']
        self.header.ndim = h1['ndim']
        self.header.nlevelmax = h1['nlevelmax']
        self.header.nboundary = h1['nboundary']
        self.header.gamma = h1['gamma']
        self.header.nvarh_org = h1['nvarh']

    def load(self, pure=False, **kwargs):
        if pure:
            amr2cell_py(self, **kwargs)
        else:
            self.amr2cell(**kwargs)

    def amr2cell(self, lmax=None, icpu=0, cpu=True, ref=False,
                 verbose=False, return_meta=False, fields=None,
                 ranges=None, nvarh=None, additional_field=None):
        """
        Loads AMR and HYDRO and output hydro data into particle-like format(cell).

        Parameters
        ----------
        cpu : bool, optional
            If True, cpu number of each cell is stored.
        icpu : int, array-like, optional
            list of cpus to load, has no effect...
        lmax : int, optional
            Limit the maximum level of hydro data returned.
        return_meta : bool, optional
            If True, returns meta data instead. (Why would I want that??)
        verbose : bool, optional

        """
        if verbose: print("[hydro.amr2cell], self.cpus = ", self.cpus)
        if nvarh is None:
            if self.header.nvarh is None:
                nvarh = self.header.nvarh_org
                self.header.nvarh = nvarh
            else:
                nvarh = self.header.nvarh

        nlevelmax = self.header.nlevelmax
        if lmax is None:
            lmax = nlevelmax

        if verbose: print('[hydro.amr2cell] >>> working resolution (lmax) =', lmax)

        if ranges is not None: self.set_ranges(ranges=ranges)

        xmi, xma = self.ranges[0]
        ymi, yma = self.ranges[1]
        zmi, zma = self.ranges[2]

        work_dir = self.info.base + '/' + self.info.data_dir + 'output_' + str(self.info.nout).zfill(5)
        if verbose:
            print("[hydro.amr2cell] Ranges", xmi, xma, ymi, yma)
            print("[hydro.amr2cell] Ranges", xmi, xma, zmi,zma)
            print("[hydro.amr2cell] cpus", self.cpus)

        if verbose: print("[hydro.amr2cell] before a2c_count..  lmax =", lmax)

        ngridtot, nvarh = a2c.a2c_count(work_dir, xmi, xma, ymi, yma, zmi, zma, lmax, self.cpus)
        if verbose: print("[hydro.amr2ell] a2c_count done")
        if verbose: print("[hydro.amr2cell]ranges", xmi, xma, ymi, yma, zmi, zma)
        #return
        if return_meta:
            return (ngridtot, nvarh, work_dir, xmi, xma, ymi, yma, zmi, zma, lmax)
        else:
            xarr = np.zeros((ngridtot,3), order="F", dtype=np.float64)
            varr = np.zeros((ngridtot,nvarh), order="F", dtype=np.float64)
            dxarr = np.zeros(ngridtot, order="F", dtype=np.float64)
            cpuarr = np.zeros(ngridtot, order="F", dtype=np.int32)
            refarr = np.zeros(ngridtot, order="F", dtype=np.int32)
            a2c.a2c_load(xarr, dxarr, varr, cpuarr, refarr, work_dir, xmi, xma, ymi, yma, zmi, zma,\
                                lmax, ngridtot, nvarh, self.cpus)
            # nvarh + 2 because fortran counts from 1, and nvarh=5 means 0,1,2,3,4,5.
            # nvarh_org of NH = 7 : rho, u,v,w, T, metal, zoomin_tag

            dtype_cell = {'pos': (('<f8', (3,)), 0),
                            'x': (('<f8', 1), 0),
                            'y': (('<f8', 1), 8),
                            'z': (('<f8', 1), 16),
                           'dx': (('<f8', 1), 24),
                         'var0': (('<f8', 1), 32),
                          'rho': (('<f8', 1), 32),
                          'vel': (('<f8', (3,)), 40),
                           'vx': (('<f8', 1), 40),
                           'vy': (('<f8', 1), 48),
                           'vz': (('<f8', 1), 56),
                         'var1': (('<f8', 1), 40),
                         'var2': (('<f8', 1), 48),
                         'var3': (('<f8', 1), 56),
                         'var4': (('<f8', 1), 64),
                         'temp': (('<f8', 1), 64),
                         'var5': (('<f8', 1), 72),
                        'metal': (('<f8', 1), 72)}
            dt_off = 72
            if cpu:
                dtype_cell.update({'cpu': (('<f8',1),dt_off+8)})
                dt_off += 8 
            if ref:
                dtype_cell.update({'ref': (('bool',1),dt_off+4)})
                dtype_cell.update({'var6': (('bool',1),dt_off+4)})
                dt_off += 4
            if additional_field is not None:
                """
                Todo: Make last_var_now a managed attribute
                """
                last_var_now = [s for s in dtype_cell.keys() if 'var' in s][-1]
                last_var_now = int(last_var_now.split('var')[1])
                if isinstance(additional_field, dict):
                    additional_field = [additional_field]
                try:
                    for af in additional_field:
                        itemsize = np.dtype(af["dtype"]).itemsize*af["dim"]
                        last_var_now += 1
                        dtype_cell.update({f'var{last_var_now}':((af["dtype"],af["dim"]),dt_off+itemsize)})
                        dtype_cell.update({af["name"]:((af["dtype"],af["dim"]),dt_off+itemsize)})
                        dt_off += itemsize
                except:
                    print("additional field is given, but can't update the dtype")
                    print("Excepted: {'name':'asdf', 'dtype':'<f8', 'dim':3}")
                    print("Or, array/list of such dicts")
                    print("Currently given : ", additional_field)

        self.cell = np.zeros(ngridtot, dtype=dtype_cell)
        self.cell['pos'] = xarr
        self.cell['dx'] = dxarr
        # Ignore the last hydro variable, which is zoom-in flag

        for i in range(nvarh-1):
            self.cell['var' + str(i)] = varr[:,i]
        if cpu:
            self.cell['cpu'] = cpuarr#cell[3]
        if ref:
            self.cell["ref"] = refarr#cell[4]
        #xarr=0; dxarr=0; cpuarr=0; refarr=0; varr=0
        del xarr, varr, dxarr, cpuarr, refarr
        gc.collect()


    #def amr2cell_py(self, lmax=None,
