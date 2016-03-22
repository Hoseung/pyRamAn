# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 14:51:42 2015

Routines used in sampling tasks


@author: hoseung
"""
def region_from_halo(halo, ind=None, rscale = 1.0):
    if ind == None:
        ind = range(halo.data)
    
    return set_region(xc = halo.data['x'],
                      yc = halo.data['y'],
                      zc = halo.data['z'],
                      radius = halo.data['rvir'] * rscale)

def distance_to(xc, xx):
    import numpy as np
    return np.sqrt([(xc[0] - xx[0])**2 + (xc[1] - xx[1])**2 + (xc[2] - xx[2])**2])[0]

def extract_halos_within(halos, i_center, scale=1.0, Mcut=1e5):
    import numpy as np
    '''
    Return indices of halos within SCALE * Rvir of the central halo.

    def extract_halos_within(halos, ind_center, scale=1.0)
    halos : halo finder output (single snapshot)
    ind_center : index of the central halo
    scale : multiplying factor to the Rvir of the central halo
    
    NOTE
    ----
    Additional halo selection criteria are likely to be added onward. 
    Therefore it is better to return mask array (True, True, ..., False, True)
    than a subset of indices.    
    
    Another "extract_halos_within" is in lambda_mp_xxx.py series. 
    That function accepts an explicit radius cut to work with galaxies. (rvir of galaxy is much smaller than rvir of halo.)
    
    
    To do
    -----
    TreeMaker fields are : ['p'][0:3],
    whereas HaloMaker gives ['x'], ['y'], ['z'].
    standardize it.
    '''
    xc = halos['x'][i_center]
    yc = halos['y'][i_center]
    zc = halos['z'][i_center]
    rvir= halos['rvir'][i_center]

    xx = halos['x']
    yy = halos['y']
    zz = halos['z']
    #m = np.array(halos['m'])
    dd = distance_to([xc,yc,zc], [xx,yy,zz]) 

    return (dd < (rvir * scale)) * (np.array(halos['m']) > Mcut)


def set_region_multi(xc=None, yc=None, zc=None, radius=None, **kwargs):
    """
    Return region dict.
    
    Note
    ----
    Accepts single / multiple inputs. 
    If scalar variables are given, set_region is directly called. 
    """
    
    import numpy as np
    if isinstance(xc, float):
        return set_region(xc=xc, yc=yc, zc=zc, radius = radius, **kwargs)
    else:        
        ranges = sum_region_multi(xc=xc, yc=yc, zc=zc, radius=radius)
        rr = max([np.ptp(ranges[0:2]),
                  np.ptp(ranges[2:4]),
                  np.ptp(ranges[4:6])]) * 0.5
            
        return set_region(xc=0.5 * sum(ranges[0:2]),
                            yc=0.5 * sum(ranges[2:4]),
                            zc=0.5 * sum(ranges[4:6]), radius = rr)


def sum_region_multi(xc=[], yc=[], zc=[], radius = []):
    """
    Returns minimum and maximum of given ranges. 
    
    Note
    ----
    
    """
      
    xmi=[]
    xma=[]
    ymi=[]
    yma=[]
    zmi=[]
    zma=[]
    for i in range(len(xc)):
        xmi.append(xc[i] - radius[i])
        xma.append(xc[i] + radius[i])
        ymi.append(yc[i] - radius[i])
        yma.append(yc[i] + radius[i])
        zmi.append(zc[i] - radius[i])
        zma.append(zc[i] + radius[i])
    
    return (min(xmi), max(xma), min(ymi), max(yma), min(zmi), max(zma))

    
def set_region(**region):
    """
        take any of xr, yr, zr, xc, yc, zc, ranges or centers
        and calculate the others.
        returns a dict (or defaultdict).

        example
        >>> r1 = set_region(yr=[0.45,0.733])
        >>> r2 = set_region(xc=0.7, radius=0.2)
        >>> r3 = set_region(abc=123)
        >>> c = {"centers":[0.33, 0.33, 0.64], "radius":0.1}
        >>> r4 = set_region(*c)

        To do
        1. Allow mix of xc and xr arguments.
           as xr requires more effort, xr overwrites xc.
        2. Accept both "center" and "centers"
           Also, both "range" and "ranges"
    """
    
    xxr = "xr" in region
    xyr = "yr" in region
    xzr = "zr" in region
    xxc = "xc" in region
    xyc = "yc" in region
    xzc = "zc" in region
    xrd = "radius" in region
    xcens = "centers" in region
    xrans = "ranges" in region
    
#    print(xxr, xyr, xzr, xxc,xyc,xzc,xrd,xcens,xrans)
    if not xrans:
        try:
            xr = region["xr"]
        except:
            xr = [0,1]
        try:
            yr = region["yr"]
        except:
            yr = [0,1]
        try:
            zr = region["zr"]
        except:
            zr = [0,1]
    else:
        ranges = region["ranges"]
        xr, yr, zr = ranges[0], ranges[1], ranges[2]

    if not xcens:
        if xxc:
            xc = region["xc"]
        else:
            xc = 0.5
        if xyc:
            yc = region["yc"]
        else:
            yc = 0.5
        if xzc:
            zc = region["zc"]
        else:
            zc = 0.5
        if xrd:
            radius = region["radius"]
        else:
            radius = 0.5
    else:
        centers = region["centers"]
        xc, yc, zc = centers#[0], centers[1], centers[2]
        radius = region["radius"]

    if (not xxr) and (not xyr) and (not xzr) and (not xrans):
        no_ranges = True
    else:
        no_ranges = False

    if (not xxc) and (not xyc) and (not xzc) and (not xcens) or (not xrd):
        no_centers = True
    else:
        no_centers = False

# fill empty values
    if no_ranges is True:
        if no_centers:
            print("No useful data is given.. 1")
            print("Setting to default")
            return {"ranges":[[0,1]]*3, "centers":[0.5]*3, "radius":0.5}
        else:
            # at least one value is given
            try:
                ranges=[[xc - radius, xc + radius],
                        [yc - radius, yc + radius],
                        [zc - radius, zc + radius]]
                centers=[xc, yc, zc]
            except:
                print("Something is wrong 3")
                print("Given values")
                print("xc:",xc, "yc:",yc, "zc:",zc)
                return False
    elif no_centers is True:
        # But ranges is not None
        try:
            ranges = [xr, yr, zr]
            centers=[sum(xr)*0.5, sum(yr)*0.5, sum(zr)*0.5]
            radius=max([ranges[0][1]-ranges[0][0], ranges[1][1]-ranges[1][0], ranges[2][1]-ranges[2][0]])
        except:
            print("Something is wrong")
            print("Given values")
            print("ranges:", ranges)
            return False

    region={"ranges":ranges, "centers":centers, "radius":radius}
    region2={"xr":ranges[0], "yr":ranges[1], "zr":ranges[2], \
             "xc":centers[0], "yc":centers[1], "zc":centers[2]}
    region.update(region2)
    return region
#r1 = set_region(yr=[0.45,0.733])
#r2 = set_region(xc=0.7, radius=0.2)
#r3 = set_region(center=[0.33, 0.33, 0.64], radius=0.1)
#reg={"ranges":[[0.3, 0.4], [0.3, 0.4], [0.3, 0.4]]}
#region = set_region(*reg)

#%%

class Region():
    def __init__(self):

        pass

