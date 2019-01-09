import numpy as np
from utils.fancy_print import printw

dtype_halo_ahf = np.dtype([('id','<i4'),
            ('hostHalo','<i4'),
                ('nSub','<i4'),
                ('mvir','<f4'),
                  ('np','<i4'),
                   ('x','<f4'),
                   ('y','<f4'),
                   ('z','<f4'),
                  ('vx','<f4'),
                  ('vy','<f4'),
                  ('vz','<f4'),
                ('rvir','<f4'),
                ('Rmax','<f4'),
                  ('R2','<f4'),
             ('mbp_off','<f4'),
             ('com_off','<f4'),
                ('Vmax','<f4'),
               ('V_esc','<f4'),
                ('sigV','<f4'),
              ('lambda','<f4'),
             ('lambdaE','<f4'),
                ('Lx','<f4'),
                ('Ly','<f4'),
                ('Lz','<f4'),
                ('b','<f4'),
                ('c','<f4'),
                ('Eax','<f4'),
                ('Eay','<f4'),
                ('Eaz','<f4'),
                ('Ebx','<f4'),
                ('Eby','<f4'),
                ('Ebz','<f4'),
                ('Ecx','<f4'),
                ('Ecy','<f4'),
                ('Ecz','<f4'),
           ('overdens','<f4'),
              ('nbins','<i4'),
            ('fMhires','<f4'),
               ('Ekin','<f4'),
                ('Epot','<f4'),
                ('SurfP','<f4'),
                ('Phi0','<f4'),
                ('cNFW','<f4'),
                ('n_gas','<i4'),
                ('M_gas','<f4'),
                ('lambda_gas','<f4'),
                ('lambdaE_gas','<f4'),
                ('Lx_gas','<f4'),
                ('Ly_gas','<f4'),
                ('Lz_gas','<f4'),
                ('b_gas','<f4'),
                ('c_gas','<f4'),
                ('Eax_gas','<f4'),
                ('Eay_gas','<f4'),
                ('Eaz_gas','<f4'),
                ('Ebx_gas','<f4'),
                ('Eby_gas','<f4'),
                ('Ebz_gas','<f4'),
                ('Ecx_gas','<f4'),
                ('Ecy_gas','<f4'),
                ('Ecz_gas','<f4'),
                ('Ekin_gas','<f4'),
                ('Epot_gas','<f4'),
                ('n_star','<i4'),
                ('M_star','<f4'),
                ('lambda_star','<f4'),
                ('lambdaE_star','<f4'),
                ('Lx_star','<f4'),
                ('Ly_star','<f4'),
                ('Lz_star','<f4'),
                ('b_star','<f4'),
                ('c_star','<f4'),
                ('Eax_star','<f4'),
                ('Eay_star','<f4'),
                ('Eaz_star','<f4'),
                ('Ebx_star','<f4'),
                ('Eby_star','<f4'),
                ('Ebz_star','<f4'),
                ('Ecx_star','<f4'),
                ('Ecy_star','<f4'),
                ('Ecz_star','<f4'),
                ('Ekin_star','<f4'),
                ('Epot_star','<f4')])

class AHF_halo():
    def __init__(self, fn=None, load=True, return_id=False):
        """
        Assumes the files follow the default AHF naming scheme.

        To be merged with Tipsy class because a snapshot and a halo catalogue
        are always closely related!
        """

        self.data = None
        self._has_data = False

        self.set_fn(fn)# = fn

        if load: self.load_ahf_data()

        if return_id is False:
            self.return_id = False
            self._return_id_list = None
        else:
            if hasattr(return_id, "__len__"):
                # String also has __len__, but let's just ignore such cases.
                self._return_id_list = return_id
            else:
                self._return_id_list = None # None = load all halo's ids.
            self.return_id = True

        if self.return_id:
            self.load_halo_member()

    def set_fn(self, fn):
        """
        TODO
        AHF filename has redshift in it, which adds an additional layer of
        complexity to guess the file name based on the snapshot base name.

        Note that the file name consists of three parts.
        snapshot_base + redshift + AHF_output_type

        """
        self._fn = fn # Whatever that is...
        last = fn.split(".")[-1]
        if "_halos" in last:
            self._fn_base = fn.split(last)[0]
        else:
            try:
                int(last)
                self._fn_base = fn
            except IOError as e:
                print("file name not understood:" + str(e))

    def load_ahf_data(self, fn=None):
        if fn is None:
            self._fn = self._fn_base + "AHF_halos"

        self.data = np.genfromtxt(self._fn, dtype=dtype_halo_ahf, skip_header=1)
        self._has_data = True

    def physical_unit(self):
        """

        """

    def load_halo_member(self, fn=None):
        """
        Load paticle membership from ".AHF_particles" file.

        Needs the halo catalogue, self.data.

        To match is with the simulation data,
        I need a complete particle list or the particle ID.
        """
        if not self._has_data:
            printw("[AHF.load_halo_member] Load catalog data first.")
            return

        if fn is None:
            # always assum hcat._fn is available.
            # Otherwise asserting npart will fail anyway.
            fn = self._fn.replace("_halos", "_particles")

        self.ind_list=[]

        inds = []
        types = []

        with open(fn, "r")  as f:
            f.readline()
            for hid, hnp in zip(self.data["id"], self.data["np"]):
                l = f.readline()
                nparts = int(l.split()[0])
                assert nparts == hnp

                if self._return_id_list is not None:
                    if hid in self._return_id_list:
                        chunk = np.fromfile(f, dtype=int, sep=" ", count=nparts * 2).reshape(nparts, 2)
                        self.ind_list.append(dict(hid=hid,
                        idgas = chunk[chunk[:,1]==0,0],
                        iddm = chunk[chunk[:,1]==1,0],
                        idstar = chunk[chunk[:,1]==4,0]))
                    else:
                        # Skip n lines. This is a fast way.
                        for nn in npart: next(f)
                else:
                    chunk = np.fromfile(f, dtype=int, sep=" ", count=nparts * 2).reshape(nparts, 2)
                    self.ind_list.append(dict(hid=hid,
                                              idgas = chunk[chunk[:,1]==0,0],
                                              iddm = chunk[chunk[:,1]==1,0],
                                              idstar = chunk[chunk[:,1]==4,0]))

    def load_subs(self, fn=None):
        if not self._has_data:
            printw("[AHF.load_sub] Load catalog data first.")
            return

        if fn is None:
            # always assum hcat._fn is available.
            # Otherwise asserting npart will fail anyway.
            fn = self._fn.replace("_halos", "_substructure")

        sub_ids = []
        hosts = []
        with open(fn, "r") as f:
            while True:
                sss = f.readline().split()
                try:
                    hid, nsub = int(sss[0]), int(sss[1])
                except:
                    break
                subs = np.array(f.readline().split()).astype(np.int32)
                assert len(subs) == nsub
                sub_ids.append(subs)
                hosts.append(hid)

        self.subs = dict(zip(hosts, sub_ids))

    def to_physical():
        """
        convert into physical units and save the state
        """
