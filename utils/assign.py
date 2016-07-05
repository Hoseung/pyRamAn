
# coding: utf-8

import numpy as np

def cic(value, x, nx, y=None, ny=1, z=None, nz=1,
        wraparound=False, average=False, debug=False,
        norm_integer=False):
    """ Interpolate an irregularly sampled field using Cloud in Cell
    method.

    This function interpolates an irregularly sampled field to a
    regular grid using Cloud In Cell (nearest grid point gets weight
    1-dngp, point on other side gets weight dngp, where dngp is the
    distance to the nearest grid point in units of the cell size).

    Inputs
    ------
    value: array, shape (N,)
        Sample weights (field values). For a temperature field this
        would be the temperature and the keyword average should be
        True. For a density field this could be either the particle
        mass (average should be False) or the density (average should
        be True).
    x: array, shape (N,)
        X coordinates of field samples, unit indices: [0,NX>.
    nx: int
        Number of grid points in X-direction.
    y: array, shape (N,), optional
        Y coordinates of field samples, unit indices: [0,NY>.
    ny: int, optional
        Number of grid points in Y-direction.
    z: array, shape (N,), optional
        Z coordinates of field samples, unit indices: [0,NZ>.
    nz: int, optional
        Number of grid points in Z-direction.
    wraparound: bool (False)
        If True, then values past the first or last grid point can
        wrap around and contribute to the grid point on the opposite
        side (see the Notes section below).
    average: bool (False)
        If True, average the contributions of each value to a grid
        point instead of summing them.

    Returns
    -------
    dens: ndarray, shape (nx, ny, nz)
        The grid point values.

    Notes
    -----
    Example of default allocation of nearest grid points: nx = 4, * = gridpoint.

      0   1   2   3     Index of gridpoints
      *   *   *   *     Grid points
    |---|---|---|---|   Range allocated to gridpoints ([0.0,1.0> -> 0, etc.)
    0   1   2   3   4   posx

    Example of ngp allocation for wraparound=True: nx = 4, * = gridpoint.

      0   1   2   3        Index of gridpoints
      *   *   *   *        Grid points
    |---|---|---|---|--    Range allocated to gridpoints ([0.5,1.5> -> 1, etc.)
      0   1   2   3   4=0  posx


    References
    ----------
    R.W. Hockney and J.W. Eastwood, Computer Simulations Using Particles
        (New York: McGraw-Hill, 1981).

    Modification History
    --------------------
    IDL code written by Joop Schaye, Feb 1999.
    Avoid integer overflow for large dimensions P.Riley/W.Landsman Dec. 1999
    Translated to Python by Neil Crighton, July 2009.

    Examples
    --------
    >>> nx = 20
    >>> ny = 10
    >>> posx = np.random.rand(size=1000)
    >>> posy = np.random.rand(size=1000)
    >>> value = posx**2 + posy**2
    >>> field = cic(value, posx*nx, nx, posy*ny, ny)
    # plot surface
    """
    from utils.cic_module import update_field_vals

    def findweights(pos, ngrid):
        """ Calculate CIC weights.
        Coordinates of nearest grid point (ngp) to each value. """

        if wraparound:
            # grid points at integer values
            ngp = np.fix(pos + 0.5)
        else:
            # grid points are at half-integer values, starting at 0.5,
            # ending at len(grid) - 0.5
            ngp = np.fix(pos) + 0.5

        # Distance from sample to ngp.
        distngp = ngp - pos

        # weight for higher (right, w2) and lower (left, w1) ngp
        weight2 = np.abs(distngp)
        weight1 = 1.0 - weight2

        # indices of the nearest grid points
        if wraparound:
            ind1 = ngp
        else:
            ind1 = ngp - 0.5
        ind1 = ind1.astype(int)

        ind2 = ind1 - 1
        # Correct points where ngp < pos (ngp to the left).
        ind2[distngp < 0] += 2

        # Note that ind2 can be both -1 and ngrid at this point,
        # regardless of wraparound. This is because distngp can be
        # exactly zero.
        bad = (ind2 == -1)
        ind2[bad] = ngrid - 1
        if not wraparound:
            weight2[bad] = 0.
        bad = (ind2 == ngrid)
        ind2[bad] = 0
        if not wraparound:
            weight2[bad] = 0.

        if wraparound:
            ind1[ind1 == ngrid] = 0

        return dict(weight=weight1, ind=ind1), dict(weight=weight2, ind=ind2)




    nx, ny, nz = (int(i) for i in (nx, ny, nz))
    nxny = nx * ny
    value = np.asarray(value)

    print( 'Resampling %i values to a %i by %i by %i grid' \
		% (len(value), nx, ny, nz))

    # normalise data such that grid points are at integer positions.
    # WHY??
    if norm_integer:
        x = (x - min(x)) / x.ptp() * nx
        if y is not None:
            y = (y - min(y)) / y.ptp() * ny
        if z is not None:
            z = (z - min(z)) / z.ptp() * nz
    
    x1, x2 = findweights(np.asarray(x), nx)
    y1 = z1 = dict(weight=1., ind=0)
    if y is not None:
        y1, y2 = findweights(np.asarray(y), ny)
        if z is not None:
            z1, z2 = findweights(np.asarray(z), nz)

    # float32 to save memory for big arrays (e.g. 256**3)
    field = np.zeros(nx * ny * nz, np.float32)

    if average:
        totalweight = np.zeros(nx * ny * nz, np.float32)
    else:
        totalweight = None
        

    update_field_vals(field, totalweight, x1, y1, z1,
            nx, ny, value, average)
    update_field_vals(field, totalweight, x2, y1, z1,
            nx, ny, value, average)
    if y is not None:
        update_field_vals(field, totalweight, x1, y2, z1,
                nx, ny, value, average)
        update_field_vals(field, totalweight, x2, y2, z1,
                nx, ny, value, average)
        if z is not None:
            update_field_vals(field, totalweight, x1, y1, z2,
                    nx, ny, value, average)
            update_field_vals(field, totalweight, x2, y1, z2,
                    nx, ny, value, average)
            update_field_vals(field, totalweight, x1, y2, z2,
                    nx, ny, value, average)
            update_field_vals(field, totalweight, x2, y2, z2,
                    nx, ny, value, average)

    if average:
        good = totalweight > 0
        field[good] /= totalweight[good]

    return field.reshape((nx, ny, nz)).squeeze()

    def update_field_vals_old(field, totalweight, a, b, c, value, debug=True, norm_integer=False):
        """ This updates the field array (and the totweight array if
        average is True).

        The elements to update and their values are inferred from
        a,b,c and value.
        """
        #print('Updating field vals')
        # indices for field - doesn't include all combinations
        indices = a['ind'] + b['ind'] * nx + c['ind'] * nxny
        # weight per coordinate
        weights = a['weight'] * b['weight'] * c['weight']
        # Don't modify the input value array, just rebind the name.
        value = weights * value
        if average:
            for i,ind in enumerate(indices):
                field[ind] += value[i]
                totalweight[ind] += weights[i]
        else:
            for i, ind in enumerate(indices):
                field[ind] += value[i]
                if debug: print(ind, weights[i], value[i], field[ind])
