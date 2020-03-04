import numpy as np

def cal_norm_vec(pop, dest=[0., 0., 1], bound_percentile=50):

    return np.sum(np.cross(pop["pos"],pop["vel"]), axis=0)

def _get_rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    http://stackoverflow.com/a/6802723
    https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula
    """

    #import math
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2)
    b, c, d = -axis*np.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def cal_rotation_matrix(nvec=None, dest=[0., 0., 1]):
    """
    Determine Rotation Matrix by axis-angle rotation.

    parameters
    -------------
    nvec : normal vector.
        If not given, automatically calculated.
    dest : destination vector ( not axis of rotation!).
           +z direction ([0, 0, 1]) by default.
    """
    import numpy.linalg as lag
    import math
    # rotation axis
    dest = np.asarray(dest) # towards +z direction.

    if lag.norm(nvec) != 1.0:
        nvec = nvec / lag.norm(nvec)
    if lag.norm(dest) != 1.0:
        dest = dest / lag.norm(dest)

    print(nvec, dest)
    r_axis = np.cross(nvec, dest)
    angle = math.acos(np.dot(nvec, dest))

    return _get_rotation_matrix(r_axis, angle)
