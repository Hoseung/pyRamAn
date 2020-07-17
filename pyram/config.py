import numpy as np

# avaiable modes: none, ng, nh
output_format = 'output_{snap.iout:05d}'

info_format = {
    'ng': 'info.txt',
}
info_format.update(dict.fromkeys(['nh', 'nh_dm_only', 'none', 'yzics', 'yzics_dm_only', 'iap', 'gem', 'fornax'], 'info_{snap.iout:05d}.txt'))

data_format = {
    'ng': '{{type}}.out{{icpu:05d}}',
}

sinkprop_format = 'sink_{icoarse:05d}.dat'

data_format.update(dict.fromkeys(['nh', 'nh_dm_only', 'none', 'yzics', 'yzics_dm_only', 'iap', 'gem', 'fornax'], '{{type}}_{snap.iout:05d}.out{{icpu:05d}}'))

default = [('x', 'f8'), ('y', 'f8'), ('z', 'f8'), ('vx', 'f8'), ('vy', 'f8'), ('vz', 'f8'), ('m', 'f8')]

# columns for particle table, see readr.f90
part_dtype = {
    'yzics': default + [('time', 'f8'), ('metal', 'f8'), ('id', 'i4'), ('level', 'u1'), ('cpu', 'i4')],
    'yzics_dm_only': default + [('id', 'i4'), ('level', 'u1'), ('cpu', 'i4')],

    'nh': default + [('time', 'f8'), ('metal', 'f8'), ('id', 'i4'), ('level', 'u1'), ('cpu', 'i4')],
    'nh_dm_only' : default + [('id', 'i4'), ('level', 'u1'), ('cpu', 'i4')],

    'none': default + [('time', 'f8'), ('id', 'i4'), ('level', 'u1'), ('cpu', 'i4'), ('family', 'i1'), ('tag', 'i1')],
    'iap': default + [('time', 'f8'), ('metal', 'f8'), ('id', 'i4'), ('level', 'u1'), ('cpu', 'i4'), ('family', 'i1'), ('tag', 'i1')],
    'gem': default + [('time', 'f8'), ('metal', 'f8'), ('id', 'i4'), ('level', 'u1'), ('cpu', 'i4'), ('family', 'i1'), ('tag', 'i1')],
    'fornax': default + [('time', 'f8'), ('metal', 'f8'), ('id', 'i4'), ('level', 'u1'), ('cpu', 'i4'), ('family', 'i1'), ('tag', 'i1')],
    'gem_longint': default + [('time', 'f8'), ('metal', 'f8'), ('id', 'i8'), ('level', 'u1'), ('cpu', 'i4'), ('family', 'i1'),('tag', 'i1')],

    'ng': default + [('id', 'i4'), ('level', 'u1'), ('cpu', 'i4')],
}

sink_prop_dtype_drag = [
    ('id', 'i4'), ('n_star', 'i4'), ('n_dm', 'i4'), ('m', 'f8'), ('x', 'f8'), ('y', 'f8'), ('z', 'f8'), ('vx', 'f8'), ('vy', 'f8'), ('vz', 'f8'),
    ('gas_jx', 'f8'), ('gas_jy', 'f8'), ('gas_jz', 'f8'), ('Mdot', 'f8'), ('Medd', 'f8'), ('dM', 'f8'),
    ('d_avgptr', 'f8'), ('c_avgptr', 'f8'), ('v_avgptr', 'f8'), ('Esave', 'f8'),
    ('jx', 'f8'), ('jy', 'f8'), ('jz', 'f8'), ('spinmag', 'f8'), ('eps_sink', 'f8'),
    ('rho_star', 'f8'), ('rho_dm', 'f8'), ('star_vx', 'f8'), ('star_vy', 'f8'), ('star_vz', 'f8'), ('dm_vx', 'f8'), ('dm_vy', 'f8'), ('dm_vz', 'f8'),
    ('low_star', 'f8'), ('low_dm', 'f8'), ('fast_star', 'f8'), ('fast_dm', 'f8')
]

sink_prop_dtype_drag_fornax = [
    ('id', 'i4'), ('n_star', 'i4'), ('n_dm', 'i4'), ('m', 'f8'), ('x', 'f8'), ('y', 'f8'), ('z', 'f8'), ('vx', 'f8'), ('vy', 'f8'), ('vz', 'f8'),
    ('gas_jx', 'f8'), ('gas_jy', 'f8'), ('gas_jz', 'f8'), ('time', 'f8'), ('Mdot', 'f8'), ('Medd', 'f8'), ('dM', 'f8'),
    ('d_avgptr', 'f8'), ('c_avgptr', 'f8'), ('v_avgptr', 'f8'), ('Esave', 'f8'),
    ('jx', 'f8'), ('jy', 'f8'), ('jz', 'f8'), ('spinmag', 'f8'), ('eps_sink', 'f8'),
    ('rho_star', 'f8'), ('rho_dm', 'f8'), ('star_vx', 'f8'), ('star_vy', 'f8'), ('star_vz', 'f8'), ('dm_vx', 'f8'), ('dm_vy', 'f8'), ('dm_vz', 'f8'),
    ('low_star', 'f8'), ('low_dm', 'f8'), ('fast_star', 'f8'), ('fast_dm', 'f8')
]

sink_prop_dtype = [
    ('id', 'i4'), ('m', 'f8'), ('x', 'f8'), ('y', 'f8'), ('z', 'f8'), ('vx', 'f8'), ('vy', 'f8'), ('vz', 'f8'),
    ('gas_jx', 'f8'), ('gas_jy', 'f8'), ('gas_jz', 'f8'), ('Mdot', 'f8'), ('Medd', 'f8'), ('dM', 'f8'),
    ('d_avgptr', 'f8'), ('c_avgptr', 'f8'), ('v_avgptr', 'f8'), ('Esave', 'f8'),
    ('jx', 'f8'), ('jy', 'f8'), ('jz', 'f8'), ('spinmag', 'f8'), ('eps_sink', 'f8')
]

sink_table_dtype = [('id', 'i8'), ('m', 'f8'), ('x', 'f8'), ('y', 'f8'), ('z', 'f8'), ('vx', 'f8'), ('vy', 'f8'), ('vz', 'f8')]

# columns for hydro quantity table, all float64, see readr.f90
hydro_names = {
    'nh': ['rho', 'vx', 'vy', 'vz', 'P', 'metal', 'refmask'],
    'nh_dm_only': ['rho', 'vx', 'vy', 'vz', 'P', 'metal', 'refmask'],
    'yzics': ['rho', 'vx', 'vy', 'vz', 'P', 'metal'],
    'yzics_dm_only': ['rho', 'vx', 'vy', 'vz', 'P', 'metal'],
    'none': ['rho', 'vx', 'vy', 'vz', 'P'],
    'iap': ['rho', 'vx', 'vy', 'vz', 'P', 'metal', 'refmask'],
    'gem': ['rho', 'vx', 'vy', 'vz', 'P', 'metal', 'refmask'],
    'fornax': ['rho', 'vx', 'vy', 'vz', 'P', 'metal', 'refmask'],
    'ng': ['rho', 'vx', 'vy', 'vz', 'P'],
}

part_family = {
    'sink_tracer': -3,
    'star_tracer': -2,
    'gas_tracer': 0,
    'dm': 1,
    'star': 2,
    'cloud': 3,
    'sink': 3,
    'tracer': [-3, -2, -1, 0],
}

dim_keys = ['x', 'y', 'z']
vel_keys = ['vx', 'vy', 'vz']

# some preferences
progress_bar_limit = 100

verbose_default = 1
#timer = Timer(verbose=verbose_default)
default_box = np.array([[0, 1], [0, 1], [0, 1]])


from .constants import Msol, km, pc, kpc, Mpc, ly, yr, Myr, Gyr, k_B, m_H, gamma
# custom units used for RAMSES snapshot, conversion from code unit to conventional unit.
def custom_units(snap):
    params = snap.params
    l = params['unit_l']
    m = params['unit_m']
    t = params['unit_t']
    d = params['unit_d']

    snap.unit = {
        # Length
        'cm'  : 1E0 / l,
        'm'   : 1E2 / l,
        'km'  : 1E5 / l,
        'pc'  : pc / l,
        'kpc' : kpc / l,
        'Mpc' : Mpc / l,
        'ly'  : ly / l,

        # Mass
        'g'  : 1 / m,
        'kg'  : 1E3 / m,
        'Msol': Msol / m,
        'Msun': Msol / m,

        # Time
        'yr'  : yr / t,
        'Myr' : Myr / t,
        'Gyr' : Gyr / t,

        # Density
        'g/cc': 1E0 / d,
        'H/cc': m_H / d,
        'Msol/Mpc3': Msol / Mpc ** 3 / d,
        'Msol/kpc3': Msol / kpc ** 3 / d,

        # Velocity
        'km/s': 1E5 * t / l,
        'cm/s': t / l,

        # Temperature
        'K'   : t ** 2 / l ** 2 * k_B / m_H,

        # Column density
        'Msol/Mpc2': Msol / Mpc ** 2 / m * l ** 2,
        'Msol/kpc2': Msol / kpc ** 2 / m * l ** 2,
        'Msol/pc2': Msol / pc ** 2 / m * l ** 2,
        'H/cm2': m_H / m * l ** 2,

        # Pressure
        'Ba'  : t**2 * l / m,

        # Flux
        'Msol/yr': Msol / yr / m * t,

        None  : 1
    }
