# coding: utf-8

# Reading halo bricks
# read binary data (halo properties + particle IDs?)
# def read_binary
# with open(halo_dir + fn_binary, 'rb') as f:
# f.readinto(a) save directly into an array (of the same shape, of course)


def iscomment(s):
   """
   Tests if the argument is a string.
    
   Parameters
   ----------
   s : string        
   """
   return s.startswith('#')
# read ascii data
'''
   #id num_p mvir mbound_vir rvir vmax rvmax vrms x y z vx vy vz Jx Jy Jz E
   #Spin PosUncertainty VelUncertainty bulk_vx bulk_vy bulk_vz BulkVelUnc n_core
   #m200b m200c m500c m2500c Xoff Voff spin_bullock b_to_a c_to_a A[x] A[y] A[z]
   #b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) Rs Rs_Klypin T/|U|
   #M_pe_Behroozi M_pe_Diemer idx i_so i_ph num_cp mmetric
   #a = 0.910000
   #Bounds: (0.000000, 0.000000, 0.000000) - (199.632004, 139.822800, 95.927795)
   #Om = 0.272000; Ol = 0.728000; h = 0.704000
   #FOF linking length: 0.280000
   #Unbound Threshold: 0.500000; FOF Refinement Threshold: 0.700000
   #Particle mass: 5.82735e+07 Msun/h
   #Box size: 199.632011 Mpc/h
   #Total particles processed: 4623859
   #Force resolution assumed: 0.0005 Mpc/h
   #Units: Masses in Msun / h
   #Units: Positions in Mpc / h (comoving)
   #Units: Velocities in km / s (physical, peculiar)
   #Units: Halo Distances, Lengths, and Radii in kpc / h (comoving)
   #Units: Angular Momenta in (Msun/h) * (Mpc/h) * km/s (physical)
   #Units: Spins are dimensionless
   #Units: Total energy in (Msun/h)*(km/s)^2 (physical)
   #Note: idx, i_so, and i_ph are internal debugging quantities
   #Np is an internal debugging quantity.
   #Rockstar Version: 0.99.9-RC3
'''

def read_halo_ascii(filename, return_aexp=False):
    import warnings
    """ Loads a single Rockstar halo output in ASCII format """
    import numpy as np
    datatype = [ 'i8', 'i8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',\
                 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',\
                 'f8', 'f8', 'f8', 'f8', 'f8', 'i8', 'f8', 'f8', 'f8', 'f8',\
                 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',\
                 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'i8', 'i8',\
                 'i8', 'i8', 'f8' ]

    names=('id', 'np', 'mvir', 'mvirb', 'rvir', 'vmax', 'rvmax', 'vrms',
           'x', 'y', 'z', 'vx', 'vy', 'vz', 'jx', 'jy', 'jz',
           'E', 'spin', 'err_x', 'err_v',
           'bulk_vx', 'bulk_vy', 'bulk_vz', 'err_bulkV',
           'n_core', 'm200b', 'm200c', 'm500c', 'm2500c',
           'xoff', 'voff', 'spin_Bullock', 'btoc', 'ctoa',
           'ax', 'ay', 'az', 'btoa500', 'ctoa500',
           'ax500', 'ay500', 'az500', 'rs', 'rs_Klypin', 'toveru',
           'mpe_b', 'mpe_d', 'idx', 'i_so', 'i_ph', 'num_cp', 'mmetric' )

#from itertools import dropwhile, filterfalse
    with open(filename, 'r') as f:
        f.readline()        
        line = f.readline()
        aexp = float(line.split()[-1])        
        #for line in f.dropwhile(iscomment,f):
            #len_file = sum(1 for _ in ifilterfalse(iscomment,f))
            # underscore serves the role of a 'throwaway' variable name
            # to indicate that part of a function result is being
            # deliberately ignored.
            # In this case, _ tells that the loop variable is not going to be used.
            # count number of lines that are not iscomment.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
#      data = np.loadtxt(myfile, unpack=True)
            data = np.genfromtxt(filename, dtype = datatype, names = names, comments='#')
#    data = np.rec.fromrecords(data, names=', '.join(data.dtype.names), dtype=data.dtype)
        data = np.rec.fromrecords(np.hstack(data), dtype=data.dtype)
# Because, unlike tree**.dat, header is clearly separated from the data,
# we can just skip header using the 'comment' option.
    if return_aexp:
        return data, aexp
    else:
        return data


def read_halo_all(filename, sort=False, return_aexp=False, verbose=False, rs_dir='rockstar_halos/'):
    """
    Loads a set of ascii outputs by calling read_halo_ascii.
    Rockstar output for a single snapshot is distributed over multiple files.
    """
## glob does not return error when there is no matching file name.
## it just returns empty list.
## So try - except part makes no sense.
    from glob import glob
    import numpy as np

    # check if the filename is a string
    print(filename.split(sep=rs_dir + 'halos_'))
    former, latter = filename.split(sep=rs_dir + 'halos_')
    nout = latter.split('.')[0]
    try:
        fn_list = sorted(glob(former + rs_dir + 'halos_' + str(nout) + '.*.ascii'))
        if verbose: print(fn_list)
    except:
        print('Invalid file name is given:')
        print(fn_list)
        print('at least, up to DOT is required. -> ..../halos_32.')
    # you don't know in what order the result will be given.
    # It's system specific.
    if return_aexp:
        halo, aexp = read_halo_ascii(fn_list[0], return_aexp=True)
    else:
        halo = read_halo_ascii(fn_list[0], return_aexp=False)
    
    for fn in fn_list[1:]:
        halo = np.append(halo, read_halo_ascii(fn, return_aexp=False))

    # Note that multi-core Rockstar outputs are not strictly ordered.
    # Order of file name does not imply order in IDs.
    if (sort):
        halo.sort(order = 'id')

    if return_aexp:
        return halo, aexp
    else:
        return halo


def get_main_prg(tree, halos, nout_ini=None, nout_fi=None):
    import numpy as np
    inout = np.where(tree[:]["NOUT"] == 0)[0]
    prg_list=np.zeros([len(halos), abs(nout_ini - nout_fi)+1], dtype=int)

    for ihalo, halo in enumerate(halos):
        i_prg = np.where(tree[inout]["HALNUM"] == halo)[0]
        prg_idx = tree[inout[i_prg]]['IDX'][0]
        print(prg_idx)
        prg_list[ihalo][0] = prg_idx
        for i in range(nout_fi, nout_ini):
        # index = idx
            prg_idx = tree[prg_idx]["TREE"][0]
            prg_list[ihalo][i + 1]=prg_idx
            print(prg_idx, tree[prg_idx]["NOUT"])
        # First element of "tree" is the main progenitor.
        # I assume....

    return prg_list


def pickleRS(base='./', nouts=[], nout_ini=1, nout_fi=None):
    """
    Loads Rockstar halo outputs and pickle it.
    output directory = base_dir/rhalo/halos_py    
    """
    from os.path import isdir
    from os import path
    from os import mkdir
    import pickle

    if nout_fi is not None:
        nouts=range(nout_ini,nout_fi)

    if len(nouts) == 0:
        import glob
        flist = glob.glob(base + 'rhalo/rockstar_halos/halos_*.0.ascii')
        nout_list=[]        
        for f in flist:
            t = f.split('halos_')[1]
            nout_list.append(int(t.split('.')[0]))  

    dir_out = path.normpath(base + 'rhalo/halos_py')
    if not isdir(dir_out):
        mkdir(dir_out)
    
    from tree import rshalo
    for i in nouts:
        fn = base + 'rhalo/rockstar_halos/halos_' + str(i) + '.0'
        halo = rshalo.read_halo_all(fn, sort=True)        
        with open(dir_out + '/halos_' + str(i).zfill(3) + '.pickle' , 'wb') as f:
            pickle.dump(halo, f, protocol = pickle.HIGHEST_PROTOCOL)
            
