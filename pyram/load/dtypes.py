import numpy as np

"""
    functions to manipulate dtypes of numpy arrays.

    NOTE
    ----
    Field ordering:
        Because numpy arrays are mainly accessed by their field name, there's no
        need to worry about the ordering. Thus, it is enough to put a new field
        to the end of the existing field list.
"""

def add_dtypes(old_dtypes, new_dtypes):
    """
    append new field, or alias to current dtypes.
    old_dtypes may be np.dtype or an input to make a np.dtype.

    add_dtypes must be composed of 5 elements:
        name of the field,
        data type,
        shape,
        reference field name,
        offset from the reference field name.

    example
    -------
    >>> old_dtypes = [('nstep', '<i4'),
                  ('id', '<i4'),
                  ('m', dtype_float)]

    >>> add_dtypes = [("x", dtype_float, 1, "xp", 0),
                  ("y", dtype_float, 1, "xp", 8),
                  ("z", dtype_float, 1, "xp", 16),
                  ("vx", dtype_float, 1, "vp", 0),
                  ("vy", dtype_float, 1, "vp", 8),
                  ("vz", dtype_float, 1, "vp", 16)]

    >>> new_dt = add_dtypes(old_dtypes, add_dtypes)

    If the name of reference field is empty ("" or None),
    the entry is added to the end of the current field list and the offset is ignored.

    To do
    -----
    adapt to new_dtype with minimal information so that
    add_dtypes =[("a1", "<f4"), ("a2", "<i4")]
               == [("a1", "<f4", 1, "", 0), 
                   ("a2", "<i4", 1, "", 0)]

    """
    if not isinstance(old_dtypes, np.dtype):
        old_dtypes = np.dtype(old_dtypes)

    dtype_new =dict(old_dtypes.fields)

    # reorder so that the un-referenced fields appear at last.
    i=0
    simply_add=[]
    for nf in new_dtypes:
        if nf[2] == "" or nf[2] == None:
            simply_add.append(new_dtypes.pop(i))
        else:
            i +=1
    new_dtypes += simply_add

    for nf in new_dtypes:
        if nf[3] in old_dtypes.fields.keys():
            # If the reference field is found.
            offset = old_dtypes.fields.get(nf[3])[1]
            dtype_new.update({nf[0]: ((nf[1], nf[2]), offset+nf[4])})
        else:
            # If not, it could be a default field input
            try:
                # The "last" entry == largest offet
                key_max = max(dtype_new, key=(lambda key: dtype_new[key][-1]))
                off_last = dtype_new[key_max][-1]
                size_last = np.dtype(dtype_new[key_max][0]).itemsize
                dtype_new.update({nf[0]: ((nf[1], nf[2]), off_last+size_last)})
            except:
                print("[load/dtype error] Can't find the field ({}) in the old dtypes".format(nf[3]))
                raise

    return dtype_new


def get_halo_dtype(is_gal=False, double=False, read_mbp=False, new_fields=None, auto_add_field=True):
    """
    Returns dtypes needed for halo catalog.

    Parameters
    ----------
    is_gal : False
        If reading a galaxy catalog, additional fields are used.
    double : False
        NH Galaxy catalogs are in double precision due to the very deep zoom-in.
    read_mbp : False
        Modified version of HaloMaker provides the ID of the most bound particle
        to be used in building trees. (In fact, it's not exactly the most bound
        particle, but something representative of the halo dynaics)
    new_fields : None
        New fields can be added to the halo catalog by using the "add_dtypes" function.
        The new field must be in the following format:
        (name, type, shape, existing field, offset from the existing field)
        For example,
        ("pos", dtype_float, (3,), "x", 0),
        ("vel", dtype_float, (3,), "vx", 0),
        ("lvec", dtype_float, (3,), "ax", 0)
    """
    if double:
        dtype_float = "<f8"
        dtype_float2 = dtype_float # Did sp version have a few dp value?
        df=8
    else:
        dtype_float = "<f4"
        dtype_float2 = dtype_float # Did sp version have a few dp value?
        df=4

    dtype_halo = [('np', '<i4'), ('id', '<i4'), ('level', '<i4'),
                  ('host', '<i4'), ('sub', '<i4'), ('nsub', '<i4'),('nextsub', '<i4'), 
                  #('aexp', dtype_float),
                  ('m', dtype_float), 
                  ('x', dtype_float), ('y', dtype_float), ('z', dtype_float),
                  ('vx', dtype_float), ('vy', dtype_float), ('vz', dtype_float),
                  ('ax', dtype_float), ('ay', dtype_float), ('az', dtype_float),
                  ('r', dtype_float), ('a', dtype_float2),('b', dtype_float2),('c', dtype_float2),
                  ('ek', dtype_float2),('ep', dtype_float2),('et', dtype_float2),
                  ('sp', dtype_float), 
                  ('rvir', dtype_float),('mvir', dtype_float),
                  ('tvir', dtype_float),('cvel', dtype_float),
                  ('p_rho', dtype_float),('p_c', dtype_float),
                  ]
                  # always 8bytes??

    if is_gal:
        dtype_halo.insert(25, ('mbulge', dtype_float))
        dtype_halo.insert(25, ('sigbulge', dtype_float))
        dtype_halo.insert(25, ('sig', dtype_float))
        dtype_halo[3] = ('hosthalo', '<i4') # Substitute 
                       #('g_nbin', '<i4'), ('g_rr', dtype_float, (100,)),
                       #('g_rho', dtype_float, (100,))]

    if read_mbp:
        dtype_halo += [('mbp', '<i8')]

    # Doesn't work with readthm.
    if auto_add_field:
        add_dtype = [("pos", dtype_float, (3,), "x", 0),
                 ("vel", dtype_float, (3,), "vx", 0),
                 ("lvec", dtype_float, (3,), "ax", 0)]

        if new_fields is not None:
            add_dtype += new_fields

        return add_dtypes(dtype_halo, add_dtype)
    else:
        return dtype_halo


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


"""
family: 
0 : gas tracer
1 : DM
-1: DM tracer
2 : star
-2: star tracer
3 : sink
-3: sink tarcer

"""

default = [('x', 'f8'), ('y', 'f8'), ('z', 'f8'), ('vx', 'f8'), ('vy', 'f8'), ('vz', 'f8'), ('m', 'f8')]

part_dtype = {
    'yzics': default + [('time', 'f8'), ('metal', 'f8'), ('id', 'i4'), ('level', 'u1'), ('cpu', 'i4')],
    'yzics_dm_only': default + [('id', 'i4'), ('level', 'u1'), ('cpu', 'i4')],

    'nh': default + [('time', 'f8'), ('metal', 'f8'), ('id', 'i4'), ('level', 'u1'), ('cpu', 'i4')],
    'nh_dm_only' : default + [('id', 'i4'), ('level', 'u1'), ('cpu', 'i4')],

    'none': default + [('time', 'f8'), ('id', 'i4'), ('level', 'u1'), ('cpu', 'i4'), ('family', 'i1'), ('tag', 'i1')],
    'iap': default + [('time', 'f8'), ('metal', 'f8'), ('id', 'i4'), ('level', 'u1'), ('cpu', 'i4'), ('family', 'i1'), ('tag', 'i1')],
    'gem': default + [('time', 'f8'), ('metal', 'f8'), ('id', 'i4'), ('level', 'u1'), ('cpu', 'i4'), ('family', 'i1'), ('tag', 'i1')],
    'fornax': default + [('time', 'f8'), ('metal', 'f8'), ('id', 'i4'), ('level', 'u1'), ('cpu', 'i4'), ('family', 'i1'), ('tag', 'i1')],
    'gem_longint': default + [('time', 'f8'), ('metal', 'f8'), ('id', 'i8'), ('level', 'u1'), ('cpu', 'i4'), ('family', 'i1'), ('tag', 'i1')],

    'ng': default + [('id', 'i4'), ('level', 'u1'), ('cpu', 'i4')],
}




fornax_dtypes={"tracer":{'id': (('<i8', 1), 0),
               'pos': (('<f8', (3,)), 8),
                 'x': (('<f8', 1), 8),
                 'y': (('<f8', 1), 16),
                 'z': (('<f8', 1), 24),
               'vel': (('<f8', (3,)), 32),
                'vx': (('<f8', 1), 32),
                'vy': (('<f8', 1), 40),
                'vz': (('<f8', 1), 48),
               'cpu': (('<i4', 1), 56),
             'level': (('<u1', 1), 60),
               'tag': (('<i1', 1), 61)},
        "dm" : {  'id': (('<i8', 1), 0),
             'pos': (('<f8', (3,)), 8),
               'x': (('<f8', 1), 8),
               'y': (('<f8', 1), 16),
               'z': (('<f8', 1), 24),
             'vel': (('<f8', (3,)), 32),
              'vx': (('<f8', 1), 32),
              'vy': (('<f8', 1), 40),
              'vz': (('<f8', 1), 48),
               'm': (('<f8', 1), 56),
             'cpu': (('<i4', 1), 80),
            'level':(('<u1', 1), 84)},
   "star" : { 'id': (('<i8', 1), 0),
              'pos': (('<f8', (3,)), 8),
                'x': (('<f8', 1), 8),
                'y': (('<f8', 1), 16),
                'z': (('<f8', 1), 24),
              'vel': (('<f8', (3,)), 32),
               'vx': (('<f8', 1), 32),
               'vy': (('<f8', 1), 40),
               'vz': (('<f8', 1), 48),
                'm': (('<f8', 1), 56),
             'time': (('<f8', 1), 64),
            'metal': (('<f8', 1), 72),
              'cpu': (('<i4', 1), 80),
            'level': (('<u1', 1), 84),
              'tag': (('<i1', 1), 85)}, # tag 1 == young star
    "star_tracer":{'id': (('<i8', 1), 0),
               'pos': (('<f8', (3,)), 8),
                 'x': (('<f8', 1), 8),
                 'y': (('<f8', 1), 16),
                 'z': (('<f8', 1), 24),
               'vel': (('<f8', (3,)), 32),
                'vx': (('<f8', 1), 32),
                'vy': (('<f8', 1), 40),
                'vz': (('<f8', 1), 48),
               'cpu': (('<i4', 1), 56),
             'level': (('<u1', 1), 60),
               'tag': (('<i1', 1), 61)},
    "sink" : { 'id': (('<i8', 1), 0),
              'pos': (('<f8', (3,)), 8),
                'x': (('<f8', 1), 8),
                'y': (('<f8', 1), 16),
                'z': (('<f8', 1), 24),
              'vel': (('<f8', (3,)), 32),
               'vx': (('<f8', 1), 32),
               'vy': (('<f8', 1), 40),
               'vz': (('<f8', 1), 48),
                'm': (('<f8', 1), 56),
             'time': (('<f8', 1), 64),
              'cpu': (('<i4', 1), 72),
              'tag': (('<i1', 1), 76)},
       "raw_dtype":[('x', 'f8'),
                   ('y', 'f8'),
                   ('z', 'f8'),
                   ('vx', 'f8'),
                   ('vy', 'f8'),
                   ('vz', 'f8'),
                   ('m', 'f8'),
                   ('time', 'f8'),
                   ('metal', 'f8'),
                   ('id', 'i4'),
                   ('level', 'u1'),
                   ('cpu', 'i4'),
                   ('family', 'i1'),
                   ('tag', 'i1')]}

nh_dtypes={'dtype_star' : { 'id': (('<i8', 1), 0),
                          'pos': (('<f8', (3,)), 8),
                            'x': (('<f8', 1), 8),
                            'y': (('<f8', 1), 16),
                            'z': (('<f8', 1), 24),
                          'vel': (('<f8', (3,)), 32),
                           'vx': (('<f8', 1), 32),
                           'vy': (('<f8', 1), 40),
                           'vz': (('<f8', 1), 48),
                            'm': (('<f8', 1), 56),
                         'time': (('<f8', 1), 64)},
            'dtype_dm' : { 'id': (('<i8', 1), 0),
                          'pos': (('<f8', (3,)), 8),
                            'x': (('<f8', 1), 8),
                            'y': (('<f8', 1), 16),
                            'z': (('<f8', 1), 24),
                          'vel': (('<f8', (3,)), 32),
                           'vx': (('<f8', 1), 32),
                           'vy': (('<f8', 1), 40),
                           'vz': (('<f8', 1), 48),
                            'm': (('<f8', 1), 56)},
           'dtype_sink' :{ 'id': (('<i8', 1), 0),
                          'pos': (('<f8', (3,)), 8),
                            'x': (('<f8', 1), 8),
                            'y': (('<f8', 1), 16),
                            'z': (('<f8', 1), 24),
                          'vel': (('<f8', (3,)), 32),
                           'vx': (('<f8', 1), 32),
                           'vy': (('<f8', 1), 40),
                           'vz': (('<f8', 1), 48)}
                         }
                    
