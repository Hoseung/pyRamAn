# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 04:50:37 2016

@author: hoseung
"""
from load.utils import read_header, read_fortran
import numpy as np


class Header():
    def __init__(self):
        self.nout = -1
        self.gid = -1
        self.hid = -1
        

class Units():
    def __init__(self):
        self.time = ""
        self.length =""
        self.position ="" # code, center, galaxy
        self.velocity = ""
        self.density = ""
#        self.unit_system = "cgs"
  
    def set_system(self, name):
        self.time = name
        self.length =name
        self.position =name
        self.velocity = name
        self.density = name
        

class Dummy():
    def __init__(self):
        pass

class Gal():
    def __init__(self, nout, idgal, wdir='./', idhal = -1, load=True, info=None, rscale=1.5):
        """
        
        Parameters
        ----------
        load : logical (Default = True)
            load avaialable (star, dm, cell) data on creating a instance.
        
        """
        self.star = None # data -> star
        self.header = None
        self.dm = None
        self.cell = None
        self.nout = nout
        self.gid = idgal
        self.hid = idhal
        self.units = Dummy()
        self.units.star = Units()
        self.units.dm = Units()
        self.units.cell = Units()
        self.units.header = Units()
        self.wdir = wdir
        #self.info = Dummy()
        self.set_info(info)
        self.rscale=rscale
        if load:
            self.load()

    def _check_info(self):
        return hasattr(self.info, "unit_l")

    def _get_minimal_info(self, info):
        """
        exclude method.
        """
        keep_list=["msun", "unit_Z", "unit_l", "unit_d", "unit_t", "unit_v",\
                   "unit_T2", "pboxsize","kms","unit_nH",\
                   "base", "nout", "aexp","zred","h","H0","time",\
                   "ob","ol","om","tGyr","unit_flux","boxtokpc"]
#        if self.info is None:
#            self.info = Dummy()
        for name, val in info.__dict__.items():
            if name in keep_list:
                setattr(self.info, name, val)
#                print(self.info.unit_d)

    def set_info(self, info=None):
        if info is None:
            from load.info import Info
            self.info = Info(self.nout, base=self.wdir)
        #self._get_minimal_info(info) 
        # Most of the information are needed. 
        # No point filtering few attributes.


    def get_rgal(self):
        """
        get rgal as ptp() of star position. 
        If necessary, only the resulting value is converted. 
        """
        if not self._check_info():
            print("Aborting...")
        if self.units.star.length == "code":
            self.rgal = 0.5 * max([self.star['x'].ptp(), self.star['y'].ptp(), self.star['z'].ptp()])
        elif self.units.star.length == "gm":
            self.rgal = 0.5 * max([self.star['x'].ptp(), self.star['y'].ptp(), self.star['z'].ptp()]) / self.info.pboxsize
            

    def load(self, star="gm", dm="gm", cell="gm",
             info=None, rscale=None, radius=None):
        """

        Notes
        -----
        when loading raw data the choice of data 1.5 * rgal,
        where rgal is the maximum ptp() of stellar particles, is arbitrary. 
        I don't know how to GUESS where the galaxy gas component will truncate.

        Info must be available by this time.
        If not, load one.


        
        >>> GM_gal.header.dtype
        >>> dtype([('my_number', '<i4'), ('level', '<i4'), ('mgal', '<f8'), ('xg', '<f8', (3,)), ('vg', '<f8', (3,)), ('lg', '<f8', (3,)), ('npart', '<i4')])

        >>> hmo.Halo().data.dtype
        >>> dtype((numpy.record, [('np', '<i4'), ('id', '<i4'), ('level', '<i4'), ('host', '<i4'), ('sub', '<i4'), ('nsub', '<i4'), ('nextsub', '<i4'), ('m', '<f4'), ('mvir', '<f4'), ('r', '<f4'), ('rvir', '<f4'), ('tvir', '<f4'), ('cvel', '<f4'), ('x', '<f4'), ('y', '<f4'), ('z', '<f4'), ('vx', '<f4'), ('vy', '<f4'), ('vz', '<f4'), ('ax', '<f4'), ('ay', '<f4'), ('az', '<f4'), ('sp', '<f4'), ('idx', '<i4'), ('p_rho', '<f4'), ('p_c', '<f4'), ('energy', '<f8', (3,)), ('radius', '<f8', (4,))]))
        """
        from utils import sampling
        #if info is None and self.info is None:
        #    self.
        if info is not None and not hasattr(self.info, "unit_l"):
            self._get_minimal_info(info) 

        if rscale is not None:
            self.rscale = rscale
        if star == "gm":
            self.header, self.star = _rd_gal(self.nout, self.gid, wdir=self.wdir, metal=True)
            self.units.star.set_system("gm")
            if self.info is not None:
                self.center_code = self.header['xg'] / self.info.pboxsize + 0.5
                self.get_rgal()
                self.region = sampling.set_region(centers=self.center_code,
                                                  radius=self.rscale * self.rgal)
        elif star == "raw":
            # load catalog
            import tree.halomodule as hmo
            gcat = hmo.Halo(nout=self.nout, base=self.wdir, is_gal=True)
            thisgal = gcat.data[gcat.data["id"] == self.gid]
            # In which format should the header be?
            # Raw

            # Header is originally a numpy array, but a dict should be enough.
            self.header = dict(my_number = thisgal["np"],
                               mgal = thisgal["id"],
                               xg = (thisgal["x"], thisgal["y"], thisgal["z"]),
                               vg = (thisgal["vx"], thisgal["vy"], thisgal["vz"]),
                               npart = thisgal["np"])
            self.rgal = thisgal["r"] # in code unit.
            if radius is None:
                radius = self.rgal

            self.region = sampling.set_region(centers=self.header["xg"],
                                              radius=self.rscale * radius)
            from load.part import Part
            pp = Part(info=self.info, ptypes=['star id pos vel time metal'],
                      region=self.region, load=True)
            self.star = pp.star
            self.units.star.set_system("code")

        if dm == "gm":
            self.dm = rd_dm(self.nout, self.gid, wdir=self.wdir)
            self.units.dm.set_system("gm")
        elif dm == 'raw':
            from load.part import Part
            pp = Part(info=self.info, ptypes=['dm id pos vel'],
                  region=self.region, load=True)
            self.dm = pp.dm
            self.units.dm.set_system("code")
        if cell == "gm":
            self.cell = rd_cell(self.nout, self.gid, wdir=self.wdir, metal=True)
            self.units.cell.set_system("gm")
        elif cell == "raw":
            from load.hydro import Hydro
            hh = Hydro(info=self.info, region=self.region)
            hh.amr2cell()
            self.cell = hh.cell
            self.units.cell.set_system("code")


    def time2gyr(self, info=None):
        """
        Only stars have time.
        """
        if self.units.star.time == "Gyr":
            print("stellar age already in Gyr unit")
        else:
            if info is None:
                try:
                    import load
                    self.info = load.info.Info(nout=self.nout, base=self.wdir)
                    info = self.info
                except:
                    print("Failed to load info")
                    return
 
            import utils.cosmology
            self.star['time'] = utils.cosmology.time2gyr(self.star['time'],
                                 z_now = info.zred,
                                 info=info)
            self.unit.time = "Gyr"
    

    def gm2code(self, header=True, star=True, dm=True, cell=True):
        if header == True and self.header is not None:
            if self.unit.header.position == "gm":
                try:
                    self.header['xg'] = self.header['xg'] / self.info.pboxsize + 0.5
                    self.unit.header.position = "code"
                    self.unit.header.length = "code"
                except AttributeError:
                    print("No .header attribute")
            else:
                print("Header position not in GM unit")
        if star == True and self.star is not None:
            if self.unit.star.position == "gm":
                try:
                    for field in ["x", "y", "z"]:
                        self.star[field] = self.star[field] / self.info.pboxsize + 0.5
                    self.unit.star.position = "code"
                    self.unit.star.length = "code"
                except AttributeError:
                    print("No .star attribute")
            else:
                print("star position not in GM unit")
        if dm == True and self.dm is not None:
            if self.unit.dm.position == "gm":
                try:
                    for field in ["x", "y", "z"]:
                        self.dm[field] = self.dm[field] / self.info.pboxsize + 0.5
                    self.unit.dm.position = "code"
                    self.unit.dm.length = "code"
                except AttributeError:
                    print("No .dm attribute")
            else:
                print("dm position not in GM unit")
        if cell == True and self.cell is not None:
            if self.unit.cell.position == "gm":
                try:
                    for field in ["x", "y", "z"]:
                        self.cell[field] = self.cel[field] / self.info.pboxsize + 0.5
                    self.unit.cell.position = "code"
                except AttributeError:
                    print("No .cell attribute")
            else:
                print("cell position not in GM unit")
            if self.unit.cell.length == "gm":
                self.cell['dx'] *= self.info.pboxsize
                self.unit.cell.length = "code"


def _rd_gal(nout, idgal, wdir="./", metal=True,
          nchem=0, long=True, fname=None):
    """
    Parameters
    ----------
    nout : int
        Snashot number
    idgal : int
        GalaxyMaker id. (not idx)

    Examples
    --------
        >>> gal = rd_gal(187, 123)

    Notes
    -----
    header xg in Mpc (physical, centered at 0.5, 0.5, 0.5 of the simualtion volume)
    """
    if fname is None:
        idgal = str(idgal).zfill(7)
        dir_nout = "GAL_" + str(nout).zfill(5)
        fname = wdir + 'GalaxyMaker/' +  dir_nout + '/gal_stars_' + idgal
   
    #header, data =
    return rd_gm_star_file(fname)
 

def rd_gal(nout, idgal, wdir="./", metal=True,
          nchem=0, long=True, fname=None):
    """
    Parameters
    ----------
    nout : int
        Snashot number
    idgal : int
        GalaxyMaker id. (not idx)

    Examples
    --------
        >>> gal = rd_gal(187, 123)

    Notes
    -----
    header xg in Mpc (physical, centered at 0.5, 0.5, 0.5 of the simualtion volume)
    """
    if fname is None:
        idgal = str(idgal).zfill(7)
        dir_nout = "GAL_" + str(nout).zfill(5)
        fname = wdir + 'GalaxyMaker/' +  dir_nout + '/gal_stars_' + idgal
   
    print("[rd_GM.rd_gal] fname=", fname)
    header, data = rd_gm_star_file(fname)
        
    gal = Gal(nout, idgal, wdir=wdir, load=False)
    gal.star = data
    gal.header = header
    gal.gid = header['my_number']
    gal.units.star.set_system("gm")
    gal.wdir = wdir
    return gal


def rd_dm(nout, idgal, wdir="./", long=True, fname=None):
    """

    header xg in Mpc (physical, centered at 0.5, 0.5, 0.5 of the simualtion volume)
    """
    if fname is None:
        idgal = str(idgal).zfill(7)
        dir_nout = "HAL_" + str(nout).zfill(5)      
        fname = wdir + 'halo/' +  dir_nout + '/hal_dms_' + idgal
    return rd_gm_dm_file(fname)   


def rd_gm_dm_file(fname, long=True):
    # Header structure
    dtype_header = np.dtype([('my_number', 'i4'),
                             ('level', 'i4'),
                             ('mgal', 'f8'),
                             ('xg', 'f8', (3,)),
                             ('vg', 'f8', (3,)),
                             ('lg', 'f8', (3,)),
                             ('npart', 'i4')])
    
    # variable type
    dtype_data=[('x', '<f8'),
                ('y', '<f8'),
                ('z', '<f8'),
                ('vx', '<f8'),
                ('vy', '<f8'),
                ('vz', '<f8'),
               ('id', '<i4'),
               ('m', '<f8')]
    
    with open(fname, "rb") as f:
        header = read_header(f, dtype=dtype_header)
        header['mgal'] *= 1e11 # mass fof galaxy
        
        # data array
        data = np.recarray(header["npart"], dtype=dtype_data)
        data['x'] = read_fortran(f, np.dtype('f8'), header["npart"]) # row-major
        data['y'] = read_fortran(f, np.dtype('f8'), header["npart"])
        data['z'] = read_fortran(f, np.dtype('f8'), header["npart"])
    
        data['vx'] = read_fortran(f, np.dtype('f8'), header["npart"]) # row-major
        data['vy'] = read_fortran(f, np.dtype('f8'), header["npart"])
        data['vz'] = read_fortran(f, np.dtype('f8'), header["npart"])
        
        
        data['m'] = read_fortran(f, np.dtype('f8'), header["npart"]) # row-major
        data['id'] = read_fortran(f, np.dtype('i4'), header["npart"])
        
    return data


def rd_gm_star_file(fname, metal=True, nchem=0, long=True):
    # Header structure
    dtype_header = np.dtype([('my_number', 'i4'),
                             ('level', 'i4'),
                             ('mgal', 'f8'),
                             ('xg', 'f8', (3,)),
                             ('vg', 'f8', (3,)),
                             ('lg', 'f8', (3,)),
                             ('npart', 'i4')])
    
    # variable type
    dtype_data=[('x', '<f8'),
                ('y', '<f8'),
                ('z', '<f8'),
                ('vx', '<f8'),
                ('vy', '<f8'),
                ('vz', '<f8'),
                ('id', '<i4'),
                ('m', '<f8'),
                ('time', '<f8')]
    if metal:
        dtype_data.append(('metal', '<f8'))
        if nchem > 0:
            dtype_data.append(('cp', '<f8', (nchem,)))

    with open(fname, "rb") as f:
        header = read_header(f, dtype=dtype_header)
        header['mgal'] *= 1e11 # mass fof galaxy
        
        # data array
        data = np.recarray(header["npart"], dtype=dtype_data)
        data['x'] = read_fortran(f, np.dtype('f8'), header["npart"]) # row-major
        data['y'] = read_fortran(f, np.dtype('f8'), header["npart"])
        data['z'] = read_fortran(f, np.dtype('f8'), header["npart"])
    
        data['vx'] = read_fortran(f, np.dtype('f8'), header["npart"]) # row-major
        data['vy'] = read_fortran(f, np.dtype('f8'), header["npart"])
        data['vz'] = read_fortran(f, np.dtype('f8'), header["npart"])
        
        
        data['m'] = read_fortran(f, np.dtype('f8'), header["npart"]) # row-major
        data['id'] = read_fortran(f, np.dtype('i4'), header["npart"])
        data['time'] = read_fortran(f, np.dtype('f8'), header["npart"]) # row-major
        if 'metal' in data.dtype.names:
            data['metal'] = read_fortran(f, np.dtype('f8'), header["npart"])
            if nchem > 0:
                for i in range(nchem):
                    data['cp'][:,i] = read_fortran(f, np.dtype('f8'), header["npart"])
    return header, data



def rd_cell(nout, idgal, wdir="./", metal=True, nchem=0,
            fname=None):
    """
    A warpper of rd_gm_cell_file that (only) provides the file name. 
    header xg in Mpc (physical, centered at 0.5, 0.5, 0.5 of the simualtion volume)
    """
    if fname is None:
        idgal = str(idgal).zfill(7)
        dir_nout = "CELL_" + str(nout).zfill(5) 
        fname = wdir + 'GalaxyMaker/' +  dir_nout + '/gal_cells_'+ str(idgal).zfill(7)
    return rd_gm_cell_file(nout, idgal, fname, metal=metal, nchem=nchem)


def rd_gm_cell_file(nout, idgal, fname, metal=True, nchem=0):
    """
    Read GalaxyMaker format cell dump data into Recarray.
    No unit conversion performed.

    Notes
    -----
    Cell may have 'cpu' field. take care of this. 
    """
    import  utils.io_module as io
    with open(fname, 'rb') as f:
        nout0 = io.read_fortran(f, dtype=np.int32, check=False)[0]
#        assert nout == nout0, "given nout ({}) and loaded nout ({}) do not match".format(nout, nout0)
        gid = io.read_fortran(f, dtype=np.int32, check=False)[0] 
#        assert idgal == gid, "given idgal ({}) and loaded idgal ({}) do not match".format(idgal, gid)
        
        ncell = io.read_fortran(f, dtype=np.int32, check=False)[0]
        cell = np.zeros(ncell, dtype=[('x', '<f8'),('y', '<f8'),('z', '<f8'),
                                     ('dx', '<f8'),('var0', '<f8'),('var1', '<f8'),
                                     ('var2', '<f8'),('var3', '<f8'),('var4', '<f8'),
                                     ('var5', '<f8')])
        cell['x']  = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['y']  = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['z']  = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['dx'] = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['var0'] = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['var1'] = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['var2'] = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['var3'] = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['var4'] = io.read_fortran(f, dtype=np.float64, n=ncell)
        cell['var5'] = io.read_fortran(f, dtype=np.float64, n=ncell)
    return cell
