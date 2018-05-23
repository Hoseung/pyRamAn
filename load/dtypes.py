import numpy as np

def add_dtypes(old_dtypes, add_dtypes):
    """
    append new field, or alias to the current dtypes.
    old_dtypes may be np.dtype or an input to make a np.dtype.

    add_dtypes must be composed of 5 elements:
        name of the field,
        data type,
        shape,
        reference field name,
        offset from the reference field name.

    example
    -------
    old_dtypes = [('nstep', '<i4'),
                  ('id', '<i4'),
                  ('m', dtype_float)]

    add_dtypes = [("x", dtype_float, 1, "xp", 0),
                  ("y", dtype_float, 1, "xp", 8),
                  ("z", dtype_float, 1, "xp", 16),
                  ("vx", dtype_float, 1, "vp", 0),
                  ("vy", dtype_float, 1, "vp", 8),
                  ("vz", dtype_float, 1, "vp", 16)]

    new_dt = add_dtypes(old_dtyes, add_dtypes)
    """
    if not isinstance(old_dtypes, np.dtype):
        old_dtypes = np.dtype(old_dtypes)

    dtype_new =dict(old_dtypes.fields)

    for ad in add_dtypes:
        try:
            offset = old_dtypes.fields.get(ad[3])[1]
        except :
            print("[load/dtype error] Can't find the field ({}) in the old dtypes".format(ad[3]))
            raise
        dtype_new.update({ad[0]: ((ad[1], ad[2]), offset+ad[4])})

    return dtype_new


def get_halo_dtype(is_gal=False, double=False):
    if double:
        dtype_float = "<f8"
        df=8
    else:
        dtype_float = "<f4"
        df=4

    dtype_halo = [('np', '<i4'), ('id', '<i4'), ('level', '<i4'),
                  ('host', '<i4'), ('sub', '<i4'), ('nsub', '<i4'),
                  ('nextsub', '<i4'),
                  ('m', dtype_float), ('mvir', dtype_float),
                  ('r', dtype_float), ('rvir', dtype_float),
                  ('tvir', dtype_float), ('cvel', dtype_float),
                  ('x', dtype_float), ('y', dtype_float), ('z', dtype_float),
                  ('vx', dtype_float), ('vy', dtype_float), ('vz', dtype_float),
                  ('ax', dtype_float), ('ay', dtype_float), ('az', dtype_float),
                  ('sp', dtype_float), ('idx', '<i4'),
                  ('p_rho', dtype_float),('p_c', dtype_float),
                  ('energy', '<f8', (3,)), ('abc', '<f8', (3,))]

    if is_gal:
        dtype_halo += [('sig', dtype_float), ('sigbulge', dtype_float),
                       ('mbulge', dtype_float), ('hosthalo', '<i4'),
                       ('g_nbin', '<i4'), ('g_rr', dtype_float, (100,)),
                       ('g_rho', dtype_float, (100,))]

    add_dtype = [("pos", dtype_float, (3,), "x", 0),
                 ("vel", dtype_float, (3,), "vx", 0),
                 ("lvec", dtype_float, (3,), "ax", 0)]

    return add_dtypes(dtype_halo, add_dtype)


def get_tree_dtypes(BIG_RUN=False, double=True):
    if double:
        dtype_float = "<f8"
        df=8
    else:
        dtype_float = "<f4"
        df=4

    dtype_tree = [('nstep', '<i4'),
                  ('id', '<i4'),
                  ('m', dtype_float),
                  ('macc', dtype_float),
                  ('nsub', '<i4'),
                  ('xp', dtype_float, (3,)),
                  ('vp', dtype_float, (3,)),
                  ('lp', dtype_float, (3,)),
                  ('abc', dtype_float, (4,)),
                  ("ek", dtype_float),
                  ("ep", dtype_float),
                  ("et", dtype_float),
                  ("spin", dtype_float),
                  ('mvir', dtype_float),
                  ('rvir', dtype_float),
                  ('tvir', dtype_float),
                  ('cvel', dtype_float),
                  ('rho_0', dtype_float),
                  ('rs', dtype_float),
                  ('level', '<i4'),
                  ('hosthalo', '<i4'), ('hostsub', '<i4'),
                  ('nextsub', '<i4'), ('idx', '<i4'),
                  ('nprgs', '<i4'),
                  ('f_ind', '<i4'),
                  ('nsons', '<i4'),
                  ('s_ind', '<i4')]
    if not BIG_RUN:
        dtype_tree.append(("np", '<i4'))

    add_dtype = [("x", dtype_float, 1, "xp", 0),
                 ("y", dtype_float, 1, "xp", df),
                 ("z", dtype_float, 1, "xp", 2*df),
                 ("vx", dtype_float, 1, "vp", 0),
                 ("vy", dtype_float, 1, "vp", df),
                 ("vz", dtype_float, 1, "vp", 2*df)]

    return add_dtypes(dtype_tree, add_dtype)
