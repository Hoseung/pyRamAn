import numpy as np
import struct

def load_header(brick_data, double=False):
    offset = 4
    if double:
        nbytes = 8
        dtype_float="d"
    else:
        nbytes = 4
        dtype_float="f"

    nbodies = struct.unpack("i", brick_data[4:8])[0]
    offset += 12
    massp = struct.unpack(dtype_float, brick_data[offset:offset+nbytes])[0]
    offset += 8 + nbytes
    aexp = struct.unpack(dtype_float, brick_data[offset:offset+nbytes])[0]
    offset += 8 + nbytes
    omegat = struct.unpack(dtype_float, brick_data[offset:offset+nbytes])[0]
    offset += 8 + nbytes
    age = struct.unpack(dtype_float, brick_data[offset:offset+nbytes])[0]
    offset += 8 + nbytes
    halnum = struct.unpack("i", brick_data[offset:offset+4])[0]
    subnum = struct.unpack("i", brick_data[offset+4:offset+8])[0]
    return offset+16, halnum, subnum

def load_a_halo(brick_data, offset, dd, is_gal=True, double=False):
    if double:
        nbytes = 8 
        dtype_float="d"
    else:
        nbytes = 4 
        dtype_float="f"

    npart = struct.unpack("i", brick_data[offset:offset+4])[0]
    dd["np"]=npart
    offset += 12  # 12 = 4 + 8
    ids = struct.unpack_from("<{}i".format(npart), brick_data[offset:offset+4*npart])
    offset += 4*npart + 8 
    dd["id"] = struct.unpack("i", brick_data[offset:offset+4])[0]
    offset += 24
    dd["level"],dd["host"],dd["sub"],dd["nsub"],dd["nextsub"]\
    = struct.unpack_from("<5i", brick_data[offset:offset+20])
    offset += 28
    dd["m"] = struct.unpack(dtype_float, brick_data[offset:offset+nbytes])[0]
    offset += 8 + nbytes
    dd["x"],dd["y"],dd["z"] = struct.unpack_from("<3"+dtype_float, brick_data[offset:offset+3*nbytes])
    offset += 8 + 3*nbytes
    dd["vx"],dd["vy"],dd["vz"] = struct.unpack_from("<3"+dtype_float, brick_data[offset:offset+3*nbytes])
    offset += 8 + 3*nbytes
    dd["ax"],dd["ay"],dd["az"] = struct.unpack_from("<3"+dtype_float, brick_data[offset:offset+3*nbytes])
    offset += 8 + 3*nbytes
    radius= struct.unpack_from("<4"+dtype_float, brick_data[offset:offset+4*nbytes])
    dd["r"],dd["abc"] = radius[0], radius[1:]
    offset += 8 + 4*nbytes
    dd["energy"] = struct.unpack_from("<3"+dtype_float, brick_data[offset:offset+3*nbytes])
    offset += 8 + 3*nbytes
    dd["sp"] = struct.unpack(dtype_float, brick_data[offset:offset+nbytes])[0]
    offset += 8 + nbytes
    if is_gal:
        dd["sig"], dd["sigbulge"], dd["mbulge"]\
        = struct.unpack_from("<3"+dtype_float, brick_data[offset:offset+3*nbytes])
        offset += 8+ 3*nbytes
    dd["mvir"],dd["rvir"],dd["tvir"],dd["cvel"]\
    = struct.unpack_from("<4"+dtype_float, brick_data[offset:offset+4*nbytes])
    offset += 8+4*nbytes
    dd["p_rho"],dd["p_c"] = struct.unpack_from("<2"+dtype_float, brick_data[offset:offset+2*nbytes])
    offset += 8+2*nbytes
    if is_gal:
        g_nbin = struct.unpack("i", brick_data[offset:offset+4])[0]
        dd["g_nbin"]=g_nbin
        offset += 12
        dd["g_rr"] = struct.unpack_from("<{}".format(g_nbin)+dtype_float, brick_data[offset:offset+g_nbin*nbytes])
        offset += 8 + g_nbin*nbytes
        dd["g_rho"] = struct.unpack_from("<{}".format(g_nbin)+dtype_float, brick_data[offset:offset+g_nbin*nbytes])
        offset += 8 + g_nbin*nbytes

    return offset, ids


def load_hm(fn, double=True, is_gal=True, return_idlists=[]):
    """
    Return catalog in numpy array, and list of member particles in a list.
    
    >>> catalog, member_ids = load_hm("TREE_DM/tree_bricks500", is_gal=False, return_idlist=[1,3,5,7])

    
    Paramters
    ---------
    double : logical
        if True, assume real are in double precision
    is_gal : logical
        If True, read GalaxyMaker output. If False, read HaloMaker output.
    return_idlists: sequence(list, array, range, tuple)
        Give halo/galaxy ids in a list(sequence) to retrieve member particle ID of the halos.
    
    NOTE
    ----
    Reading tree_bricks in Fortranis 10x faster. 
    But, maybe it's OK to be a bit slow. NH catalogues are small, anyways.
    """
    
    if double:
        dtype_float = "<f8"
    else:
        dtype_float = "<f4"

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

    idlists=[]
    f = open(fn, "rb")
    brick_data = f.read()
    offset, halnum, subnum = load_header(brick_data, double=double)
    gcat = np.zeros(halnum+subnum, dtype=dtype_halo)
    for i in range(halnum+subnum):
        offset,_ = load_a_halo(brick_data, offset, gcat[i], is_gal=is_gal, double=double)
        if gcat[i]["id"] in return_idlists:
            idlists.append(_)
    f.close()
    
    return gcat, idlists
