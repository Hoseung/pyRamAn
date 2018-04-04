

def get_halo_dtype(is_gal=False, double=False):
    if double:
        dtype_float = "<f8"
        df=8
    else:
        dtype_float = "<f4"
        df=4

    dtype_halo = {'np': (('<i4', 1), 0),
                  'id': (('<i4', 1), 4),
                  'level': (('<i4', 1), 8),
                  'host': (('<i4', 1), 12),
                  'sub': (('<i4', 1), 16),
                  'nsub': (('<i4', 1), 20),
                  'nextsub': (('<i4', 1), 24),
                  'm': ((dtype_float, 1), 24+df),
                  'mvir': ((dtype_float, 1), 24+2*df),
                  'r': ((dtype_float, 1), 24+3*df),
                  'rvir': ((dtype_float, 1), 24+4*df),
                  'tvir': ((dtype_float, 1), 24+5*df),
                  'cvel': ((dtype_float, 1), 24+6*df),
                  'x': ((dtype_float, 1), 24+7*df),
                  'y': ((dtype_float, 1), 24+8*df),
                  'z': ((dtype_float, 1), 24+9*df),
                  'pos': ((dtype_float, 3), 24+7*df),
                  'vx': ((dtype_float, 1), 24+10*df),
                  'vy': ((dtype_float, 1), 24+11*df),
                  'vz': ((dtype_float, 1), 24+12*df),
                  'vel': ((dtype_float, 3), 24+10*df),
                  'm': ((dtype_float, 1), 24+13*df),
                  'ax': ((dtype_float, 1), 24+14*df),
                  'ay': ((dtype_float, 1), 24+15*df),
                  'az': ((dtype_float, 1), 24+16*df),
                  'lvec': ((dtype_float, 3), 24+14*df),
                  'sp': ((dtype_float, 1), 24+17*df),
                  'idx': (('<i4', 1), 24+18*df),
                  'p_rho': ((dtype_float, 1), 28+18*df),
                  'p_c': ((dtype_float, 1), 28+19*df),
                  'energy': ((dtype_float, 3), 28+20*df),
                  'abc': ((dtype_float, 3), 28+23*df)}

    if is_gal:
        dtype_halo.update({'sig': ((dtype_float, 1), 28+26*df),
                           'sigbulge': ((dtype_float, 1), 28+27*df),
                           'mbulge': ((dtype_float, 1), 28+28*df),
                           'hosthalo': (("<i4", 1), 28+29*df),
                           'g_nbin': (("<i4", 1), 32+29*df),
                           'g_rr': ((dtype_float, 100), 36+29*df),
                           'g_rho': ((dtype_float, 100), 32+129*df)})
    return dtype_halo


def get_tree_dtypes(BIG_RUN=False):
    dtype_tree = [('nstep', '<i4'),
                  ('id', '<i4'),
                  ('m', '<f8'),
                  ('macc', '<f8'),
                  ('nsub', '<i4'),
                  ('xp', '<f8', (3,)),
                  ('vp', '<f8', (3,)),
                  ('lp', '<f8', (3,)),
                  ('abc', '<f8', (4,)),
                  ("ek", '<f8'),
                  ("ep", '<f8'),
                  ("et", '<f8'),
                  ("spin", '<f8'),
                  ('mvir', '<f8'),
                  ('rvir', '<f8'),
                  ('tvir', '<f8'),
                  ('cvel', '<f8'),
                  ('rho_0', '<f8'),
                  ('rs', '<f8'),
                  ('level', '<i4'),
                  ('hosthalo', '<i4'), ('hostsub', '<i4'),
                  ('nextsub', '<i4'), ('idx', '<i4'),
                  ('nprgs', '<i4'),
                  ('f_ind', '<i4'),
                  ('nsons', '<i4'),
                  ('s_ind', '<i4')]
    if not BIG_RUN:
        dtype_tree.append(("np", '<i4'))

        if True:
            """
            This allows multiple ways of accessing fields.

            tree["x"] is equivalent with tree["pos"][:,0]
            tree["vx"] is equivalent with tree["vel"][:,0]

            Doing this requires an optional parameter "offset" to be given,
            which is calculated as offset+ad[3].
            """
            add_dtype = [("x", "f8", 1, 0, "xp"),
                         ("y", "f8", 1, 8, "xp"),
                         ("z", "f8", 1, 16, "xp"),
                         ("vx", "f8", 1, 0, "vp"),
                         ("vy", "f8", 1, 8, "vp"),
                         ("vz", "f8", 1, 16, "vp")]
            dt = np.dtype(dtype_tree)
            dtype_tree =dict(dt.fields)
            for ad in add_dtype:
                offset = dt.fields.get(ad[4])[1]
                dtype_tree.update({ad[0]: ((ad[1], ad[2]), offset+ad[3])})

    return dtype_tree
