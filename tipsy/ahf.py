import numpy as np

dtype_halo_ahf = np.dtype([('ID','<i4'),
            ('hostHalo','<i4'),
                ('nSub','<i4'),
                ('Mvir','<f4'),
                  ('np','<i4'),
                   ('x','<f4'),
                   ('y','<f4'),
                   ('z','<f4'),
                  ('vx','<f4'),
                  ('vy','<f4'),
                  ('vz','<f4'),
                ('Rvir','<f4'),
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

def load_ahf_data(fn):
    return np.genfromtxt(fn,
                  dtype=dtype_halo_ahf,
                  skip_header=1)


def load_halo_member(fn):
    """
    Load paticle membership from ".AHF_particles" file.
    """
    inds = []
    types = []

    with open(fn, "r")  as f:
        f.readline()

        for i in range(len(hal)):
            l = f.readline()
            nparts = int(l.split()[0])
            mem = np.fromfile(f, dtype=int, sep=" ", count=nparts * 2).reshape(nparts, 2)
            inds.append(mem[:, 0])
            types.append(mem[:, 1])

    return inds, types
