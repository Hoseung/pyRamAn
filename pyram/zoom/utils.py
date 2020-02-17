# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 13:37:22 2015

@author: hoseung
"""
def get_pt_id_targets2(part, h, halo_inds, r_frac=1):
    import numpy as np
    assert(np.size(halo_inds) > 0)
    buffer_frac = 1.2 * np.sqrt(r_frac)
    # search radii larger than r_vir will result in more particles.
    npart_tot_guess = np.math.ceil(sum(h.data.np[halo_inds]) * buffer_frac)
    idlist = np.zeros(npart_tot_guess, dtype=np.int32)
    x = part.x
    y = part.y
    z = part.z
    print("Expected # of particles:", npart_tot_guess)
    nhalo = len(halo_inds)
    halo_range_list = np.zeros((nhalo) + 1, dtype=np.int32)
    len_idlist = len(idlist)
    for i, ihalo in enumerate(halo_inds):
        print("# part in this halo", h.data.np[ihalo])
        xr = [h.data.x[ihalo] - h.data.rvir[ihalo] * r_frac,
              h.data.x[ihalo] + h.data.rvir[ihalo] * r_frac]
        yr = [h.data.y[ihalo] - h.data.rvir[ihalo] * r_frac,
              h.data.y[ihalo] + h.data.rvir[ihalo] * r_frac]
        zr = [h.data.z[ihalo] - h.data.rvir[ihalo] * r_frac,
              h.data.z[ihalo] + h.data.rvir[ihalo] * r_frac]
              
        ind_z = np.where((x > xr[0]) & (x < xr[1]) &
                       (y > yr[0]) & (y < yr[1]) &
                       (z > zr[0]) & (z < zr[1]))[0]
        if len(ind_z) > 0:
            print("good")
            # If array is not large enough, append it.
            halo_range_list[i + 1] = halo_range_list[i] + len(ind_z)
            if halo_range_list[i + 1] > len_idlist:
                if i + 1 == nhalo:
                    # If it's the last halo, append by exact difference.
                    npart_more = halo_range_list[i + 1] - len_idlist
                else:
                    # Otherwise, guess from previous halos.
                    # 1.5 * mean npart so far * number of remaining halos
                    npart_more = int(1.5 * (halo_range_list[i] / (i + 1)) * (nhalo - i))
                print("increase the array size by {:d} from {:d}".format(npart_more, len_idlist))
                idlist = np.append(idlist, np.zeros(npart_more, dtype=np.int32))
                len_idlist= len(idlist)

            idlist[halo_range_list[i]:halo_range_list[i+1]] = part.id[ind_z]
            print(halo_range_list[i+1], len(idlist))              

    return halo_range_list, idlist

def get_pt_id_targets(part, h, halo_inds, r_frac=1):
    import numpy as np
    """
    Returns id list of particles in each halo in a SINGLE array.
    Also returned is the list of range of id list for each halo.

    So,
    idlist1 = idslist[ind[0] : ind[1]]
    idlist2 = idslist[ind[1] : ind[2]]

    .. note::
        You don't know how many particles will be found within halos you are interested.
        A good guess would be the sum of halo.np which is number of member particles
        defined by the halo finder.

        So, I allocate the idlist array with length of 1.2 * sum(halo.np[ind]).
    """
    assert(np.size(halo_inds) > 0)
    buffer_frac = 1.3 * np.sqrt(r_frac)
    # search radii larger than r_vir will result in more particles.
    npart_tot_guess = np.math.ceil(sum(h.data.np[halo_inds]) * buffer_frac)
    idlist = np.zeros(npart_tot_guess, dtype=np.int32)
    x = part.x
    y = part.y
    z = part.z
    print("Expected # of particles:", npart_tot_guess)
    nhalo = len(halo_inds)
    halo_range_list = np.zeros((nhalo) + 1, dtype=np.int32)
    len_idlist = len(idlist)
    for i, ihalo in enumerate(halo_inds):
        print("# part in this halo", h.data.np[ihalo])
        xr = [h.data.x[ihalo] - h.data.rvir[ihalo] * r_frac,
              h.data.x[ihalo] + h.data.rvir[ihalo] * r_frac]
        yr = [h.data.y[ihalo] - h.data.rvir[ihalo] * r_frac,
              h.data.y[ihalo] + h.data.rvir[ihalo] * r_frac]
        zr = [h.data.z[ihalo] - h.data.rvir[ihalo] * r_frac,
              h.data.z[ihalo] + h.data.rvir[ihalo] * r_frac]
        ind_x = np.where( (x > xr[0]) & (x < xr[1]))[0]
        # out of 8.1GB particle information()
        if len(ind_x) > 0:
            ind_y = np.where((y[ind_x] > yr[0]) & (y[ind_x] < yr[1]))[0]
            if len(ind_y) > 0:
                ind_z = np.where((z[ind_x[ind_y]] > zr[0]) & (z[ind_x[ind_y]] < zr[1]))[0]
                # If array is not large enough, append it.
                halo_range_list[i + 1] = halo_range_list[i] + len(ind_z)
                if halo_range_list[i + 1] > len_idlist:
                    if i + 1 == nhalo:
                        # If it's the last halo, append by exact difference.
                        npart_more = halo_range_list[i + 1] - len_idlist
                    else:
                        # Otherwise, guess from previous halos.
                        # 1.5 * mean npart so far * number of remaining halos
                        npart_more = int(1.5 * (halo_range_list[i] / (i + 1)) * (nhalo - i))
                    print("increase the array size by {:d} from {:d}".format(npart_more, len_idlist))
                    idlist = np.append(idlist, np.zeros(npart_more, dtype=np.int32))
                    len_idlist= len(idlist)

                idlist[halo_range_list[i]:halo_range_list[i+1]] = part.id[ind_x[ind_y[ind_z]]]
                print(halo_range_list[i+1], len(idlist))
    
    return halo_range_list, idlist


def particle_distribution(part, idlist, buffer = 0):
    from utils import match
    """
    returns center of mass of given particles, xr,yr,zr,
    and a,b,c (in the future)
    """
    ind = match.match_list_ind(part.id, idlist)
    # To increase speed, search ind at one go and slice them.
    print(" Number of matching particles: {} out of {}".format(len(ind), len(idlist)))
    xr = [part["x"][ind].min(), part["x"][ind].max()]
    yr = [part["y"][ind].min(), part["y"][ind].max()]
    zr = [part["z"][ind].min(), part["z"][ind].max()]
    radius = max([xr[1] - xr[0], yr[1] - yr[0], zr[1] - zr[0]])
    xc = sum(xr) * 0.5
    yc = sum(yr) * 0.5
    zc = sum(zr) * 0.5

    return xc, yc, zc, radius


class Allrefine():
    def __init__(self, hid, data, aexp):
        self.id = hid
        self.x = data['x_refine']
        self.y = data['y_refine']
        self.z = data['z_refine']
        self.r = data['r_refine']
        self.aexp = aexp
        
        
def match_diff_arr(x1, x2, y1, y2, z1=None, z2=None, tolerance=None, window=0.1, reverse=False):
    """
    first array is smaller. 
    """
    import numpy as np
    n1 = len(x1)
    n2 = len(x2)
    n21 = n2/n1
    if max([n1, n2]) < 200:
        window = 0.5
    # window size depends on the length of the longer array.
    window = int(n2 * window) 
    
    if tolerance is None:
        tolerance = 0.5 * x1.ptp() / n1 # half of the mean separation in 0-th axis.

    print("window:", window)
    print("tolerance:", tolerance)        
    print("array lengths:", n1, n2)

    if reverse:
        # Suppose n1 is smaller
        sort_ind1 = np.argsort(x1)[::-1]
        sort_ind2 = np.argsort(x2)[::-1]
    else:
        sort_ind1 = np.argsort(x1)
        sort_ind2 = np.argsort(x2)

    print("n1 = {}, n2 = {}".format(n1,n2))
    # initialize matched list
    match = np.zeros(n1, dtype=int)
    in21 = [np.math.ceil(i * n21) for i in np.arange(n1)]
    for i in range(n1):
        if i < window:
            ind_1 = sort_ind2[0:in21[i] + window]
        elif i > n1 - window:
            ind_1 = sort_ind2[in21[i] - window :]
        else:
            ind_1 = sort_ind2[in21[i] - window : in21[i] + window]
        xx2 = x2[ind_1]
        xx1 = x1[sort_ind1[i]] 
        d1 = np.abs(xx2 - xx1)
        indx = np.where(d1 < tolerance)[0]
#        print(indx, len(xx2))
#        print("d1", d1[indx], xx2[indx], xx1)
        if len(indx) > 0:
            if z1 is not None and z2 is not None:
                yx2 = y2[ind_1] # first N smallest elements
                yx1 = y1[sort_ind1[i]] 
                zx2 = z2[ind_1] # first N smallest elements
                zx1 = z1[sort_ind1[i]]                 
                d2 = (xx2[indx] - xx1)**2 \
                    + (yx2[indx] - yx1)**2 + (zx2[indx] - zx1)**2
                if len(np.where(d2 < tolerance**2)[0]) > 0:
                    print(np.where(d2 < tolerance**2))
                    match[i] = sort_ind2[indx[np.argmin(d2)]] # index of minimum value
                else:
                    match[i]=-1
            elif y1 is not None and y2 is not None:
                d2 = (x2[sort_ind2[indx]] - x1[sort_ind1[i]])**2 \
                            + (y2[sort_ind2[indx]] - y1[sort_ind1[i]])**2
                if len(np.where(d2 < tolerance**2)[0]) > 0:
                    print(np.where(d2 < tolerance**2))
                    match[i] = sort_ind2[indx[np.argmin(d2)]] # index of minimum value
                else:
                    match[i]=-1
                    print("Missed!", i)
        else:
            match[i]=-1
            print("Missed!", i)

    # n-th element in the list match is the index of 2nd array corresponds to
    #        the n-th smallest element of the 1st array.
    # Or, arr1[sort_ind1] = arr2[match]
    # an index to revert sort_ind1 into the original order
    reverse_arr1 = np.argsort(sort_ind1)

    return match[reverse_arr1]

def distance3d(x1,y1,z1, x2,y2,z2):
    from numpy import sqrt
    return sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)