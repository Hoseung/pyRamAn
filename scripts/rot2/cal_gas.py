import numpy as np

def init_list_of_objects(size):
    list_of_objects = list()
    for i in range(0,size):
        list_of_objects.append( list() ) #different object reference each time
    return list_of_objects

def get_flat_data(prgs, fields):
    #
    #
    # NOTE
    # varout = [[]] * len(fields)
    # these sub lists point to the same object.
    # So, varout[0] = 0 
    # >>> [[0],[0],[0]]
    # not 
    # >>> [[0], [], []]
    # Thus, init_list_of_objects.
    #
    varout = init_list_of_objects(len(fields))
    for inout, this_step in enumerate(prgs):
        for this_sat in this_step:
            if this_sat is None:
                continue
            for istep, sat in enumerate(this_sat):
                for i,field in enumerate(fields):
                    varout[i].append(sat[field])


    return [np.array(aslit) for aslit in varout]
    #return varout
    
def dividie_cat(catalog, n_sub = 1, ind=0):
    """
    
    """
    ix = n // n_sub**2
    iy = (n - (ii*n_sub**2)) // n_sub
    iz = (n - (ii*n_sub**2)) % n_sub
    xmin = min(gcat.data["x"])
    xmax = min(gcat.data["x"])
    ymin = min(gcat.data["y"])
    ymax = min(gcat.data["y"])
    zmin = min(gcat.data["z"])
    zmax = min(gcat.data["z"])
    dx = 1./n_sub * (xmax-xmin)
    dy = 1./n_sub * (ymax-ymin)
    dz = 1./n_sub * (zmax-zmin)
    
    return gcat.data[np.where(
              (gcat.data["x"] > xmin + dx*ix) * (gcat.data["x"] < xmin + dx*(ix+1))
            * (gcat.data["y"] > ymin + dy*iy) * (gcat.data["y"] < ymin + dy*(iy+1))
            * (gcat.data["z"] > zmin + dz*iz) * (gcat.data["z"] < zmin + dz*(iz+1)))[0]]


def range_from_gcat(data, buffer=0, rscale=1.0):
    return [[min(gal['x'] - rscale * gal['r'] - buffer),
             max(gal['x'] + rscale * gal['r'] + buffer)],
            [min(gal['y'] - rscale * gal['r'] - buffer),
             max(gal['y'] + rscale * gal['r'] + buffer)],
            [min(gal['z'] - rscale * gal['r'] - buffer),
             max( gal['z'] + rscale * gal['r'] + buffer)]]

def range_from_pos_r(pos, rvir, buffer=0, rscale=1.0):
    xmin = min(pos[:,0] - rvir * rscale - buffer)
    ymin = min(pos[:,1] - rvir * rscale - buffer)
    zmin = min(pos[:,2] - rvir * rscale - buffer)

    xmax = max(pos[:,0] + rvir * rscale + buffer)
    ymax = max(pos[:,1] + rvir * rscale + buffer)
    zmax = max(pos[:,2] + rvir * rscale + buffer)
    return ([xmin, xmax], [ymin, ymax], [zmin,zmax])


def range_2code_unit(ranges, pb):
    return np.array([rr/pb + 0.5 for rrs in ranges for rr in rrs]).reshape(3,2)
