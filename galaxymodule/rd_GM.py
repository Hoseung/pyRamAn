# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 04:50:37 2016

Inherits galaxymodule.galaxy.Galaxy class

A galaxy is built with stars, DMs, and cell data bu default. The three components
are assuemd to present, but an additional functionality is desired in case one
or more components are missing.

Unlike cells, there are two possibilities in determining member particles.
Eighter use GalaxyMaker/HaloMaker idlist or take everythin within a certain radius.
As a galaxy is resolved into mulitple substructures at the NH resolution, the membership of
a "host" galaxy maybe not compatble with observations. By default, I will take
all member particles of all member substructures as the member of the "host" galaxy.

@author: hoseung
"""
from galaxymodule.galaxy import Galaxy
from utils.io_module import read_header, read_fortran
import numpy as np
import pickle

class Header():
    def __init__(self):
        self.nout = -1
        self.gid = -1
        self.hid = -1


class Units():
    """
    A class designed to track the units of gal data.

    raw, gm, relative_p

    Parameters
    ----------
    name = {"gm", "raw", "relative_p", "cgs"}
    length={"cMpc", "pMpc", "ckpc", "pkpc", "code"}
    position={"box corner", "box center", "galaxy"}
    velocity={"code", "kms"}
    mass={"Msun", "1e11Msun", "code"}
    density={"code", "cgs"}
    time = {"conformal", "Gyr"}
    metal = {"solar", "abs"}

    vel_org = {"box", "galaxy"}
    """
    def __init__(self, name=None,
                       length=None,
                       position=None,
                       velocity=None,
                       mass=None,
                       density=None,
                       time=None,
                       metal=None,
                       vel_org=None):
        """Only string or None are acceptable arguments."""
        self.name = name
        self.length = length
        self.position = position # code, center, galaxy
        self.velocity = velocity
        self.mass = mass
        self.density = density
        self.time = time
        self.metal = metal
        self.vel_org = vel_org

unit_gm = Units(name="gm",
                length="pMpc",
                position="box center",
                velocity="kms",
                mass="1e11Msun",
                density="code",
                time="conformal",
                metal="code",
                vel_org="box")

unit_rel = Units(name="relative_p",
                length="pkpc",
                position="galaxy",
                velocity="kms",
                mass="Msun",
                density="code",
                time="Gyr",
                metal="code",
                vel_org="galaxy")

unit_raw = Units(name="code",
                length="code",
                position="box corner",
                velocity="kms",
                mass="Msun",
                density="code",
                time="conformal",
                metal="code",
                vel_org="box")

unit_systems=[unit_gm, unit_rel, unit_raw]


class Dummy():
    def __init__(self):
        pass

class Gal(Galaxy):
    def __init__(self, nout, info, idgal=None,
                 gcat=None, hcat=None,
                 gcat_subs=[], hcat_subs=[],
                 type_dm="gm", type_cell="gm",
                 base='./', idhal = -1, load=True, rscale=1.5):
        """

        Parameters
        ----------
        load : logical (Default = True)
            load avaialable (star, dm, cell) data on creating the instance.
        idgal : int
            galaxy ID
        gcat : gcat.data
            gcat if avaialble.
            Well, actually it does not work properly without a gcat.
        halo : hcat.data
            If there is no CELL, DM, or GAL files avaialble, automatically fallback
            to read dm and cell from the raw data. halo information is required in that case.
            And... again, this does not work properly, yet...??
        type_dm : ["gm", "raw"], "gm" by default.
            indicate the type of intended DM.
            "gm" tries to read DM dump file, while "raw" tries to read from raw DM data.
        type_cell : ["gm", "raw"], "gm" by default.
            Same as type_dm.
            Note that there is no type_star, as GAL dump are always assuemd to be avaialble.
        idhal : defaults to -1
            ID of the mathing halo must be provided explicitly.

        wdir : "./" by default.

        load : True
            Load actual data. If false, only meta data are loaded.
        rscale : 1.5
            scale factor in grabbing galaxy data from raw DM/CELL data.

        Notes
        -----
        By default, loaded data (star, dm, cell) are in "gm" units.


        ToDo
        ----
        GAL, CELL files' naming convention is confusing. (17.08.20)

        """
        assert(not(idgal == None and gcat == None)), ("either idgal or a gcat"
        " is needed.")
        #assert(not(info == None)), "Need info, use info=gcat.info"
        if idgal is None:
            idgal = gcat["id"]

        super(Gal, self).__init__(gcat=gcat, info=info, hcat=hcat)
        try:
            self.hid = self.hcat["id"]
        except:
            self.hid = idhal
            pass
        self.star = None
        self.cell = None
        self.dm = None
        self.header = None
        self.nout = nout
        self.gid = idgal
        self.gcat_subs = gcat_subs
        self.hcat_subs = hcat_subs
        self.units = Dummy()
        self.units.star = Units()
        self.units.dm = Units()
        self.units.cell = Units()
        self.units.header = Units()
        self.wdir = wdir
        self.info = info
        self.debug=False
        #self.set_info(info)
        self.rscale=rscale
        if load:
            if idhal < 0:
                type_dm = None
            self.load(type_dm=type_dm, type_cell=type_cell)
        # try loading cell:


    def _check_info(self):
        return hasattr(self.info, "unit_l")

    def _get_minimal_info(self, info):
        """
        from an info instance, exclude methods and leave variables.
        """
        keep_list=["msun", "unit_Z", "unit_l", "unit_d", "unit_t", "unit_v",\
                   "unit_T2", "pboxsize","kms","unit_nH",\
                   "base", "nout", "aexp","zred","h","H0","time",\
                   "ob","ol","om","tGyr","unit_flux","boxtokpc"]
#
        for name, val in info.__dict__.items():
            if name in keep_list:
                setattr(self.info, name, val)

    def get_rgal(self):
        """
        Set Gal.rgal as ptp() of star position in kpc.
        If necessary, only the resulting value is converted.

        When we try to extract components from the whole data,
        the
        Rgal is required to guess the size of the region that encompasses all galactic components.

        """
        if not self._check_info():
            print("Aborting...")

        if self.units.star.name == "code":
            self.rgal = 0.5 * max([self.star['x'].ptp(), self.star['y'].ptp(), self.star['z'].ptp()])
            self.rgal_code = self.rgal
        elif self.units.star.name == "gm":
            self.rgal = 0.5 * max([self.star['x'].ptp(), self.star['y'].ptp(), self.star['z'].ptp()]) / self.info.pboxsize
            self.rgal_code = self.rgal * 1e-3
        elif self.units.star.name == "relative_p":
            self.rgal = 0.5 * max([self.star['x'].ptp(), self.star['y'].ptp(), self.star['z'].ptp()])
            self.rgal_code = self.rgal / self.info.boxtokpc


    def load(self, type_star="gm", type_dm="gm", type_cell="gm",
             info=None, rscale=None, radius=None,
             rd_star_params=None,
             rd_cell_params=None,
             rd_dm_params=None):
        """
        Load per-galaxy data (if exists).
        Automatically skips missing components.

        Parameters
        ----------
        type_star : string {"gm", "raw"}
            unit of stellar particle data, "gm" by default.
        type_dm : string {"gm", "raw"}
            unit of DM particle data, "gm" by default.
        type_cell : string {"gm", "raw"}
            unit of gas data, "gm" by default.

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
        from utils.sampling import Region

        if rscale is not None:
            self.rscale = rscale
        try:
            if type_star == "gm":
                self.header, self.star = _rd_gal(self.nout, self.gid,
                             base=self.wdir, metal=True, **rd_star_params)
                if len(self.gcat_subs) > 0:
                    ss_all = [self.star]
                    for sub in self.gcat_subs:
                        h, ss = _rd_gal(self.nout, sub["id"],
                                 base=self.wdir, metal=True, **rd_star_params)
                        ss_all.append(ss)
                    self.star = np.concatenate(ss_all)

                # convert units to kpc from the center
                self.center_code = self.header['xg'] / self.info.pboxsize + 0.5

                self.star["pos"]-=(self.center_code - 0.5)*self.info.pboxsize
                self.star["pos"]*=1e3
                self.star["m"] *= 1e11

                self.units.header = unit_gm
                self.units.star = unit_rel

                self.get_rgal()
                self.region = Region(centers=self.center_code,
                                              radius=self.rscale * self.rgal_code)
            elif type_star == "raw":
            # load catalog
                import tree.halomodule as hmo
                gcat = hmo.Halo(nout=self.nout, base=self.wdir, is_gal=True)
                thisgal = gcat.data[gcat.data["id"] == self.gid]
                # Raw

                # Header is originally a numpy array, but a dict should be enough.
                self.header = dict(my_number = thisgal["np"],
                                   mgal = thisgal["id"],
                                   xg = (thisgal["x"], thisgal["y"], thisgal["z"]),
                                   vg = (thisgal["vx"], thisgal["vy"], thisgal["vz"]),
                                   npart = thisgal["np"])
                self.units.header = unit_gm
                self.rgal = thisgal["r"] # in code unit.
                if radius is None:
                    radius = self.rgal

                self.region = Region(centers=self.header["xg"],
                                                  radius=self.rscale * radius)

                from load.part import Part
                pp = Part(info=self.info, ptypes=['star id pos vel time metal'],
                          region=self.region, load=True)
                self.star = pp.star
                self.units.star.name="code"
            self._has_star=True
        except FileNotFoundError as e:
            print("File Not Found:", e.filename)
            print("No stellar data loaded")
            self.star = None
            pass

        try:
            if type_dm == "gm":
                self.dm = rd_dm(self.nout, self.hid, base=self.wdir,
                                **rd_dm_params)
                self.units.dm.name="gm"
            elif type_dm == 'raw':
                from load.part import Part
                pp = Part(info=self.info, ptypes=['dm id pos vel'],
                      region=self.region, load=True)
                self.dm = pp.dm
                self.units.dm.name="code"
            self._has_dm=True
        except FileNotFoundError as e:
            print("File Not Found:", e.filename)
            print("No DM data loaded")
            self.dm = None
            pass

        try:
            if type_cell == "gm":
                self.cell = rd_cell(self.nout, self.gid, base=self.wdir,
                                    metal=True, **rd_cell_params)
                self.units.cell= unit_gm
                self._has_cell=True
            elif type_cell == "raw":
                from load.hydro import Hydro
                hh = Hydro(info=self.info, region=self.region)
                hh.amr2cell()
                self.cell = hh.cell
                self.units.cell= unit_raw
                self._has_cell=True
            elif type_cell == False:
                self.cell = None
                pass

        except FileNotFoundError as e:
            print("File Not Found:", e.filename)
            print("No CELL data loaded")
            self.cell = None
            pass

    def gm2code(self, header=True, star=True, dm=True, cell=True):
        if header == True and self.header is not None:
            if self.units.header.position == "gm":
                try:
                    self.header['xg'] = self.header['xg'] / self.info.pboxsize + 0.5
                    self.units.header.position = "code"
                    self.units.header.length = "code"
                except AttributeError:
                    print("No .header attribute")
            else:
                print("Header position not in GM unit")
        if star == True and self.star is not None:
            if self.units.star.position == "gm":
                try:
                    for field in ["x", "y", "z"]:
                        self.star[field] = self.star[field] / self.info.pboxsize + 0.5
                    self.units.star.position = "code"
                    self.units.star.length = "code"
                except AttributeError:
                    print("No .star attribute")
            else:
                print("star position not in GM unit")
        if dm == True and self.dm is not None:
            if self.units.dm.position == "gm":
                try:
                    for field in ["x", "y", "z"]:
                        self.dm[field] = self.dm[field] / self.info.pboxsize + 0.5
                    self.units.dm.position = "code"
                    self.units.dm.length = "code"
                except AttributeError:
                    print("No .dm attribute")
            else:
                print("dm position not in GM unit")
        if cell == True and self.cell is not None:
            if self.units.cell.position == "gm":
                try:
                    for field in ["x", "y", "z"]:
                        self.cell[field] = self.cel[field] / self.info.pboxsize + 0.5
                    self.units.cell.position = "code"
                except AttributeError:
                    print("No .cell attribute")
            else:
                print("cell position not in GM unit")
            if self.units.cell.length == "gm":
                self.cell['dx'] *= self.info.pboxsize
                self.units.cell.length = "code"



def time2gyr(self, info=None):
    """
    Only stars have time.
    """
    if self.units.star.time == "Gyr":
        print("stellar age already in Gyr unit")
    else:
        if info is None:
            if hasattr(self, "info"):
                info = self.info
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
        self.units.time = "Gyr"




def _rd_gal(nout, idgal, base="./", metal=True,
          nchem=0, long=True, fname=None, additional_fields=None):
    """
    Just read a GM file.
    All this function does is to generate a file name
    from nout and idgal.


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
    return rd_gm_star_file(fname, additional_fields=additional_fields)


def rd_gal(nout, idgal, info=None, base="./", metal=True,
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

    ...

    """
    if fname is None:
        idgal = str(idgal).zfill(7)
        dir_nout = "GAL_" + str(nout).zfill(5)
        fname = wdir + 'GalaxyMaker/' +  dir_nout + '/gal_stars_' + idgal

    print("[rd_GM.rd_gal] fname=", fname)
    header, data = rd_gm_star_file(fname)

    gal = Gal(nout, idgal, info=info, base=wdir, load=False)
    gal.star = data
    gal.header = header
    gal.gid = header['my_number']
    gal.units.star.name = "gm"
    gal.wdir = wdir
    return gal


def rd_dm(nout, idgal, base="./", long=True, fname=None):
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
    dtype_data = {'pos': (('<f8', (3,)), 0),
                  'x': (('<f8', 1), 0),
                  'y': (('<f8', 1), 8),
                  'z': (('<f8', 1), 16),
                 'id': (('<i8', 1), 24),
                  'm': (('<f8', 1), 32),
                'vel': (('<f8', (3,)), 40),
                 'vx': (('<f8', 1), 40),
                 'vy': (('<f8', 1), 48),
                 'vz': (('<f8', 1), 56)}

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


def rd_gm_star_file(fname, metal=True, nchem=0,
                    additional_fields=None):
    """
    Load one GAL file, return header and the data as is.

    Parameters
    ----------
    metal : boolean {True}
        if True, read metallicity from the raw file
    nchem : int {0}
        Number of chemical elements in addition to the total metallicity.

    additional_fields : dtype in a specific format
        appends fields to the default stellar particle fields.



    Examples
    --------
    new_dtype = {"ellip": (('<f8', 1), 0),
           "flux_r": (('<f8', 1), 8),
           "flux_g": (('<f8', 1), 16),
           "flux_b": (('<f8', 1), 24)}
    rd_star_params = dict(additional_fields=new_dtype)
    gg.load(type_cell="none", type_dm="none", rd_star_params=rd_star_params)
    print(gg.star.dtype)
        dtype((numpy.record,
        {'names':
        ['pos','x','y','z','id','m','vel','vx','vy','vz','time','metal','ellip',
        'flux_r','flux_g','flux_b'],
        'formats':
        [('<f8', (3,)),'<f8','<f8','<f8','<i8','<f8',('<f8', (3,)),'<f8','<f8',
        '<f8','<f8','<f8','<f8','<f8','<f8','<f8'],
        'offsets':
        [0,0,8,16,24,32,40,40,48,56,72,80,80,88,96,104],
        'itemsize':
        112}))



    NOTE: There is no point automatically calculating the offset as sometimes
    I want duplicate offset for multiple fields. Unless the dtype grow huge,
    it is easier to manually manage the dtype.

    I could do something like this if REALLY needed.
    dtype_data = {'pos': ('<f8', (3,)),
                'x': ('<f8', 1),
                'y': ('<f8', 1),
                'z': ('<f8', 1),
               'id': ('<i8', 1),
                'm': ('<f8', 1),
              'vel': ('<f8', (3,)),
               'vx': ('<f8', 1),
               'vy': ('<f8', 1),
               'vz': ('<f8', 1),
             'time': ('<f8', 1)}

    for field in dtype_data.items():
        dtype_data.update({field[0]: (field[1][0], d_off)})
        d_off += np.dtype(field[1][0]).itemsize * np.squeeze(field[1][1])
    """
    # Header structure
    dtype_header = np.dtype([('my_number', 'i4'),
                             ('level', 'i4'),
                             ('mgal', 'f8'),
                             ('xg', 'f8', (3,)),
                             ('vg', 'f8', (3,)),
                             ('lg', 'f8', (3,)),
                             ('npart', 'i4')])

    # variable type
    dtype_data = {'pos': (('<f8', (3,)), 0),
                    'x': (('<f8', 1), 0),
                    'y': (('<f8', 1), 8),
                    'z': (('<f8', 1), 16),
                   'id': (('<i8', 1), 24),
                    'm': (('<f8', 1), 32),
                  'vel': (('<f8', (3,)), 40),
                   'vx': (('<f8', 1), 40),
                   'vy': (('<f8', 1), 48),
                   'vz': (('<f8', 1), 56),
                 'time': (('<f8', 1), 72)}
    d_off = 72
    if metal:
        dtype_data.update({'metal': (('<f8', 1), d_off+8)})
        d_off +=8
        if nchem > 0:
            dtype_data.update({'cp': (('<f8', (nchem,)), d_off+8)})

    if additional_fields is not None:
        for field in additional_fields.items():
            additional_fields.update({field[0]: (field[1][0], field[1][1]+d_off)})
        dtype_data.update(additional_fields)

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


def rd_cell(nout, idgal, base="./", metal=True, nchem=0,
            fname=None):
    """
    A warpper of rd_gm_cell_file that (only) provides the file name.
    header xg in Mpc (physical, centered at 0.5, 0.5, 0.5 of the simualtion volume)
    """
    if fname is None:
        snout = str(nout).zfill(5)
        dir_nout = "CELL_" + snout
        try:
            fname = wdir + 'GalaxyMaker/' +  dir_nout + '/CELL_'+str(nout)+"_"+str(idgal)+".pickle"
            return pickle.load(open(fname, "rb"))
        except:
            fname = wdir + 'GalaxyMaker/' +  dir_nout + '/CELL_'+str(nout)+"_"+str(idgal).zfill(7)
            return rd_gm_cell_file(nout, idgal, fname, metal=metal, nchem=nchem)


def rd_gm_cell_file(nout, idgal, fname, metal=True, cpu=False, ref=False, nchem=0):
    """
    Read GalaxyMaker format cell dump data into Recarray.
    No unit conversion performed.

    Notes
    -----
    Cell may have 'cpu' field. take care of this.
    """
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

    import  utils.io_module as io
    with open(fname, 'rb') as f:
        nout0 = io.read_fortran(f, dtype=np.int32, check=False)[0]
#        assert nout == nout0, "given nout ({}) and loaded nout ({}) do not match".format(nout, nout0)
        gid = io.read_fortran(f, dtype=np.int32, check=False)[0]
#        assert idgal == gid, "given idgal ({}) and loaded idgal ({}) do not match".format(idgal, gid)

        ncell = io.read_fortran(f, dtype=np.int32, check=False)[0]

        cell = np.zeros(ncell, dtype=dtype_cell)
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
