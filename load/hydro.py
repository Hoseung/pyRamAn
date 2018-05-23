# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 21:32:55 2015

@author: hoseung
"""
import numpy as np
from load.sim import Simbase
from utils.io_module import read_header
import struct


def generate_fname(nout,path="",ftype="",cpuid=1,ext=""):

    if len(path) > 0:
        if path[-1] != "/":
            path=path+"/"

    if nout == -1:
        filelist = sorted(glob.glob(path+"output*"))
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

def get_binary_data(fmt="",
                    ninteg=0,
                    nlines=0,
                    nfloat=0,
                    nstrin=0,
                    nquadr=0,
                    nlongi=0,
                    content=None,
                    correction=0):

    offset = 4*ninteg + 8*(nlines+nfloat+nlongi) + nstrin + nquadr*16 + 4 + correction
    byte_size = {"i":4,"d":8,"q":8}
    if len(fmt) == 1:
        mult = 1
    else:
        mult = eval(fmt[0:len(fmt)-1])
    pack_size = mult*byte_size[fmt[-1]]

    return struct.unpack(fmt, content[offset:offset+pack_size])


class Dummy():
    def __init__(self):
        pass


class Hydro(Simbase):
    def __init__(self,
                 nout=None,
                 info=None,
                 amr=None,
                 region=None,
                 ranges=None,
                 load=False,
                 cosmo=True,
                 cpus=None,
                 cpu_fixed=None,
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
            It is important to distingush the number of hydrovaraibles to read and
            the number of hydrovaraibles in the simulation.

            Recent version of RAMSES generates "hydro_file_descriptor.txt" by default.
            This file contains nvarh and the type of each hydro-variable.
            With this file, nvarh input is not needed.


        """
        super(Hydro, self).__init__()
        self.cosmo = cosmo
        if info is None:
            assert nout is not None, "either info or onut is required"
            from load.info import Info
            print("[Hydro.__init__] Loading info")
            info = Info(nout=nout, cosmo=cosmo)
        self.info = info
        self.nout = info.nout
        self.cpus = cpus
        self.cpu_fixed=cpu_fixed
        try:
            self.ncpu = len(self.cpus)
        except:
            self.ncpu = 0

        snout = str(info.nout).zfill(5)
        # file name
        self.out_dir = 'snapshots/'
        self._fnbase = info.base + '/' + self.out_dir \
                                 + 'output_' + snout \
                                 + '/hydro_' + snout \
                                 + '.out'
        self._get_basic_info()
        self.set_info(info)
        self.header.nvarh=nvarh
        if region is not None:
            ranges = region['ranges']
        if ranges is not None:
            self.set_ranges(ranges=ranges)
        elif hasattr(info, "ranges"):
            if info.ranges is not None:
                self.set_ranges(ranges=info.ranges)
        else:
            self.set_ranges([[-9e9,9e9]] * 3)

        try:
            self.amr = amr
        except NameError:
            print("Loading amr first! \n")

        if load:
            print("amr2cell_params", amr2cell_params)
            self.amr2cell(**amr2cell_params)

    def _get_basic_info(self):
        try:
            f = open(self._fbase + '00001', "rb")
        except:
            from glob import glob
            hydros = glob(self._fnbase + "*")
            f = open(hydros[0], "rb")

        self.header = Dummy()
        self._read_hydro_header(f)

    def set_info(self, info):
        self.info = info

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
            self.amr2cell_py(**kwargs)
        else:
            self.amr2cell(**kwargs)

    def amr2cell(self, lmax=None, icpu=0, cpu=True, ref=False,
                 verbose=False, return_meta=False,
                 ranges=None, nvarh=None):
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
        print("[hydro.amr2cell], self.cpus   - 0", self.cpus)
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
        #print("[hydro.amr2cell], self.cpus   - 1", self.cpus)
        # Set ranges
        #print("[hydro.amr2cell], self.ranges - 1", self.ranges)
        xmi, xma = self.ranges[0]
        ymi, yma = self.ranges[1]
        zmi, zma = self.ranges[2]

        work_dir = self.info.base + '/' + self.out_dir + 'output_' + str(self.info.nout).zfill(5)
        if verbose:
            print("[hydro.amr2cell] Ranges", xmi, xma, ymi, yma)
            print("[hydro.amr2cell] Ranges", xmi, xma, zmi,zma)
            print("[hydro.amr2cell] cpus", self.cpus)

        from load import a2c
        if verbose: print("[hydro.amr2cell] before a2c_count..  lmax =", lmax)
        out = a2c.a2c_count(work_dir, xmi, xma, ymi, yma, zmi, zma, lmax, self.cpus)
        if verbose: print("[hydro.amr2ell] a2c_count done")
        if verbose: print("[hydro.amr2cell]ranges", xmi, xma, ymi, yma, zmi, zma)
        #return
        if return_meta:
            return (out[0], work_dir, xmi, xma, ymi, yma, zmi, zma, lmax)
        else:
            cell = a2c.a2c_load(work_dir, xmi, xma, ymi, yma, zmi, zma,\
                                lmax, out[0], nvarh+2, self.cpus)
            # nvarh + 2 because fortran counts from 1, and nvarh=5 means 0,1,2,3,4,5.
            #dtype_cell = [('x', '<f8'), ('y', '<f8'), ('z', '<f8'), ('dx', '<f8')]
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
                dtype_cell.update({'ref': (('bool',1),dt_off+8)})
            #for i in range(nvarh):
            #    dtype_cell.append( ('var' + str(i), '<f8'))

        self.cell = np.zeros(len(cell[1]), dtype=dtype_cell)
        self.cell['x'] = cell[0][:,0]
        self.cell['y'] = cell[0][:,1]
        self.cell['z'] = cell[0][:,2]
        self.cell['dx'] = cell[1]
        for i in range(nvarh):
            self.cell['var' + str(i)] = cell[2][:,i]
        if cpu:
            self.cell['cpu'] = cell[3]
        if ref:
            self.cell["ref"] = cell[4]


    def amr2cell_py(self, lmax=None,
                    icpu=0, cpu=False,
                    read_level = False,
                    read_ref = False,
                    verbose=False, debug=False, nvar_read=6):
        import load
        from utils.io_module import read_fortran, skip_fortran

        read_cpu = cpu
        # if both amr and hydro does not exist, create.
        #################
        cpus = self.cpus
        ndim = self.info.ndim
        if not verbose: print("CPU list:", cpus)
        twotondim = 2**ndim
        nvarh = self.header.nvarh
        ncpu = self.header.ncpu
        nboundary = self.header.nboundary
        nlevelmax = self.header.nlevelmax
        if lmax is None:
            lmax = nlevelmax

        print(' >>> working resolution (lmax) =', lmax)

        # Set ranges
        xmi = self.ranges[0][0]
        xma = self.ranges[0][1]
        ymi = self.ranges[1][0]
        yma = self.ranges[1][1]
        zmi = self.ranges[2][0]
        zma = self.ranges[2][1]

        # Size of arrays
        if icpu == 0:
            listmax = ncpu
        elif icpu > 0:
            listmax = ncpu + nboundary
        ngridarr = np.zeros((nlevelmax, listmax), dtype=np.int32)

        nx = self.amr.header.ng[0]
        ny = self.amr.header.ng[1]
        nz = self.amr.header.ng[2]
        xbound = [int(nx/2), int(ny/2), int(nz/2)]

        # Allocate work arrays
        twotondim = 2**self.info.ndim
        xcent = np.zeros([8,3],dtype=np.float64)
        xg    = np.zeros([self.info.ngridmax,3],dtype=np.float64)
        son   = np.zeros([self.info.ngridmax,twotondim],dtype=np.int32)
        var   = np.zeros([self.info.ngridmax,twotondim,nvar_read+6],dtype=np.float64)
        xyz   = np.zeros([self.info.ngridmax,twotondim,self.info.ndim],dtype=np.float64)
        ref   = np.zeros([self.info.ngridmax,twotondim],dtype=np.bool)

        iprog = 1
        istep = 10
        ncells_tot = 0

        # Control which variable to read here
        var_read = [True]*nvar_read

        work_dir = self.info.base + '/' + self.out_dir + 'output_{:05d}/'.format(self.info.nout)

        data_pieces = dict()
        npieces = 0

        read_header = True

        ncoarse = np.product(self.amr.header.ng)
        nboundary = self.amr.header.nboundary
        key_size = self.amr.header.key_size
        noutput = self.amr.header.nnouts

        ngridlevel = np.zeros([self.info.ncpu_tot+nboundary,self.info.lmax],dtype=np.int32)

        for kcpu in range(1, self.info.ncpu_tot+1):
            # Read binary AMR file
            fn_amr = work_dir +"amr_{:05d}.out{:05d}".format(self.info.nout, kcpu)
            if kcpu in cpus or kcpu == 1:
                with open(fn_amr, mode='rb') as famr:
                    amrContent = famr.read()

                    # Read binary HYDRO file
                    fn_hydro = fn_amr.replace("amr", "hydro")
                    with open(fn_hydro, mode='rb') as fhydro: # b is important -> binary
                        hydroContent = fhydro.read()
                        # All hydro files have the same header, right?
                        if read_header:
                            fhydro.seek(0)
                            self._read_hydro_header(fhydro)
                            read_header = False

            #print("Reading....... ", kcpu)

            nstrin = 0
            #ngridlevel = header_icpu.numbl

            # Read the number of grids
            ninteg = 14+(2*self.info.ncpu_tot*self.info.lmax)
            nfloat = 18+(2*noutput)+(2*self.info.lmax)
            nlines = 21
            ngridlevel[:self.info.ncpu_tot,:] = np.asarray(get_binary_data(fmt="%ii"%(self.info.ncpu_tot*self.info.lmax),\
             content=amrContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)).reshape(self.info.lmax,self.info.ncpu_tot).T


            # Offset for AMR
            ninteg1 = 14+(3*self.info.ncpu_tot*self.info.lmax)+(10*self.info.lmax)+(3*nboundary*self.info.lmax)+5+3*ncoarse
            nfloat1 = 18+(2*noutput)+(2*self.info.lmax)
            nlines1 = 21+2+3*min(1,nboundary)+1+1+1+3
            nstrin1 = 128 + key_size

            # Offset for HYDRO
            ninteg2 = 5
            nfloat2 = 1
            nlines2 = 6
            nstrin2 = 0

            # Loop over levels
            for ilevel in range(lmax):

                # Geometry
                dxcell=0.5**(ilevel+1)
                dx2=0.5*dxcell
                for ind in range(twotondim):
                    iz=int((ind)/4)
                    iy=int((ind-4*iz)/2)
                    ix=int((ind-2*iy-4*iz))
                    xcent[ind,0]=(float(ix)-0.5)*dxcell
                    xcent[ind,1]=(float(iy)-0.5)*dxcell
                    xcent[ind,2]=(float(iz)-0.5)*dxcell

                # Cumulative offsets in AMR file
                ninteg_amr = ninteg1
                nfloat_amr = nfloat1
                nlines_amr = nlines1
                nstrin_amr = nstrin1

                # Cumulative offsets in HYDRO file
                ninteg_hydro = ninteg2
                nfloat_hydro = nfloat2
                nlines_hydro = nlines2
                nstrin_hydro = nstrin2

                # Loop over domains
                for jcpu in range(nboundary+self.info.ncpu_tot):

                    ncache = ngridlevel[jcpu, ilevel]#header_icpu.numbl[ilevel,jcpu]

                    # Skip two lines of integers
                    nlines_hydro += 2
                    ninteg_hydro += 2

                    if ncache > 0:

                        if jcpu == kcpu-1:
                            # xg: grid coordinates
                            ninteg = ninteg_amr + ncache*3
                            nfloat = nfloat_amr
                            nlines = nlines_amr + 3
                            nstrin = nstrin_amr
                            for n in range(self.info.ndim):
                                offset = 4*ninteg + 8*(nlines+nfloat+n*(ncache+1)) + nstrin + 4
                                xg[:ncache,n] = struct.unpack("%id"%(ncache),           amrContent[offset:offset+8*ncache])

                            # son indices
                            ninteg = ninteg_amr + ncache*(4 + 2 * self.info.ndim)
                            nfloat = nfloat_amr + ncache*self.info.ndim
                            nlines = nlines_amr + 4 + 3*self.info.ndim
                            nstrin = nstrin_amr
                            for ind in range(twotondim):
                                offset = 4*(ninteg+ind*ncache) + 8*(nlines+nfloat+ind) + nstrin + 4
                                son[:ncache,ind] = struct.unpack("%ii"%(ncache), amrContent[offset:offset+4*ncache])
                                # var: hydro variables
                                jvar = 0
                                for ivar in range(nvarh):
                                    if var_read[ivar]:
                                        offset = 4*ninteg_hydro + 8*(nlines_hydro+nfloat_hydro+(ind*nvarh+ivar)*(ncache+1)) + nstrin_hydro + 4
                                        var[:ncache,ind,jvar] = struct.unpack("%id"%(ncache), hydroContent[offset:offset+8*ncache])
                                        jvar += 1

                                # var: coordinates and cell sizes
                                var[:ncache,ind,-6] = float(ilevel+1) # level
                                for n in range(self.info.ndim):
                                    xyz[:ncache,ind,n] = xg[:ncache,n] + xcent[ind,n]-xbound[n]
                                    var[:ncache,ind,-5+n] = xyz[:ncache,ind,n]*self.info.boxlen
                                var[:ncache,ind,-2] = dxcell*self.info.boxlen
                                var[:ncache,ind,-1] = kcpu # CPU number
                                # ref: True if the cell is unrefined
                                ref[:ncache,ind] = np.logical_not(np.logical_and(son[:ncache,ind] > 0,ilevel < lmax-1))

                            # Select only the unrefined cells that are in the region of interest
                            if self.info.ndim == 1:
                                cube = np.where(np.logical_and(ref[:ncache,:], \
                                                np.logical_and((xyz[:ncache,:,0]+dx2)>=xmi, \
                                                               (xyz[:ncache,:,0]-dx2)<=xma)))
                            elif self.info.ndim == 2:
                                cube = np.where(np.logical_and(ref[:ncache,:], \
                                                np.logical_and((xyz[:ncache,:,0]+dx2)>=xmi, \
                                                np.logical_and((xyz[:ncache,:,1]+dx2)>=ymi, \
                                                np.logical_and((xyz[:ncache,:,0]-dx2)<=xma, \
                                                               (xyz[:ncache,:,1]-dx2)<=yma)))))
                            elif self.info.ndim == 3:
                                cube = np.where(np.logical_and(ref[:ncache,:], \
                                                np.logical_and((xyz[:ncache,:,0]+dx2)>=xmi, \
                                                np.logical_and((xyz[:ncache,:,1]+dx2)>=ymi, \
                                                np.logical_and((xyz[:ncache,:,2]+dx2)>=zmi, \
                                                np.logical_and((xyz[:ncache,:,0]-dx2)<=xma, \
                                                np.logical_and((xyz[:ncache,:,1]-dx2)<=yma, \
                                                               (xyz[:ncache,:,2]-dx2)<=zma)))))))
                            else:
                                print("Bad number of dimensions")
                                return 0

                            cells = var[cube]
                            ncells = np.shape(cells)[0]
                            if ncells > 0:
                                ncells_tot += ncells
                                npieces += 1
                                # Add the cells in the master dictionary
                                data_pieces["piece"+str(npieces)] = cells

                        # Now increment the offsets while looping through the domains
                        ninteg_amr += ncache*(4+3*twotondim+2*self.info.ndim)
                        nfloat_amr += ncache*self.info.ndim
                        nlines_amr += 4 + 3*twotondim + 3*self.info.ndim


                        nfloat_hydro += ncache*twotondim*nvarh
                        nlines_hydro += twotondim*nvarh

                # Now increment the offsets while looping through the levels
                ninteg1 = ninteg_amr
                nfloat1 = nfloat_amr
                nlines1 = nlines_amr
                nstrin1 = nstrin_amr

                ninteg2 = ninteg_hydro
                nfloat2 = nfloat_hydro
                nlines2 = nlines_hydro
                nstrin2 = nstrin_hydro

        all_array = np.concatenate(list(data_pieces.values()), axis=0)

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
        if read_cpu:
            dtype_cell.update({'cpu': (('<f8',1),dt_off+8)})
            dt_off += 8
        if read_level:
            dtype_cell.update({'level': (('<f8',1),dt_off+8)})
            dt_off +=8
        if read_ref:
            dtype_cell.update({'ref': (('bool',1),dt_off+8)})


        # This part
        self.cell = np.zeros(len(all_array), dtype=dtype_cell)

        for i in range(nvarh):
            self.cell['var' + str(i)] = all_array[:,i]

        if read_level:
            self.cell['level'] = all_array[:,nvar_read+0]

        if read_cpu:
            self.cell['cpu'] = all_array[:,nvar_read+5]

        self.cell['x'] = all_array[:,nvar_read+1]
        self.cell['y'] = all_array[:,nvar_read+2]
        self.cell['z'] = all_array[:,nvar_read+3]
        self.cell['dx'] = all_array[:,nvar_read+4]

        if read_ref:
            self.cell["ref"] = all_array[:,7]

        #return master_data_array
