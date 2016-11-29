# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 23:05:02 2015

@author: hoseung
"""
import numpy as np

def link_circle_up(x, y, r, ax, finish=0):
    """
    Given two points, draw circle at the first point and link it to the second point
    without drawing the second point by default (so that it can repeat to build a long thread of bids).
    for the last point, pass the radius of the last circle to the argument 'finish'
    
    For example,

    fig = plt.figure()
    ax = fig.add_subplot(111)
    xpos = [1,1] &  ypos = [2,4]
    link_circle(xpos, ypos, 10, ax)

    xpos = [1,2] & ypos = [4,6]
    link_circle(xpos, ypos, 30, ax, finish=30)
    fig.show()
    """
    ax.plot(x[0], y[0], 'o', ms=r, lw=2,  alpha=0.7, mfc='orange')
    ax.plot(x, y, '-', c='black',alpha=0.7)
    if finish > 0:
        ax.plot(x[1], y[1], 'o', ms=20, lw=2, alpha=0.7, mfc='orange')    
        
        
def link_circle_up_lineandcircle(x, y, r, ax, first=0):
    """
    Given two points, draw circle at the first point and link it to the second point
    without drawing the second point by default (so that it can repeat to build a long thread of bids).
    for the last point, pass the radius of the last circle to the argument 'finish'
    
    For example,

    fig = plt.figure()
    ax = fig.add_subplot(111)
    xpos = [1,1] &  ypos = [2,4]
    link_circle(xpos, ypos, 10, ax)

    xpos = [1,2] & ypos = [4,6]
    link_circle(xpos, ypos, 30, ax, finish=30)
    fig.show()
    """
    ax.plot(x, y, '-', c='black',alpha=0.7)
    ax.plot(x[1], y[1], 'o', ms=r, lw=2,  alpha=0.7, mfc='orange')

    if first > 0:
        ax.plot(x[0], y[0], 'o', ms=first, lw=2, alpha=0.7, mfc='orange')    

def get_xarr(n):
    import numpy as np
    arr=[]
    a=0
    for i in range(n):
        a += (-1)**i * i
        arr.append(a)
    return np.asarray(arr)
       
        
def recursive_tree(idx, tt, nstep, ax, x0, y0, dx, mass_unit=1e10):
    prgs = tt[idx]["TREE"]
    m = np.sqrt(tt[idx]["M"] / mass_unit)
    #print("IDX:", idx, "prgs: ",prgs[:2], tt[idx]["NOUT"])
    prgs = prgs[prgs > 0]
    nprg = len(prgs)
    if nstep == 0:
        return 
    else:
        if nprg == 0:
            return
        else:
            if nprg > 1:
                dx *= 0.9
#                print("Branch!", nprg)

#            dxarr = np.arange(np.fix(-nprg/2), np.fix(nprg/2)+1)
#            xarr = dxarr[-nprg:]*dx + x0
            xarr = get_xarr(nprg) * dx + x0
            for i, x in zip(prgs, xarr):
                link_circle_up([x0, x], [y0, y0 + 1], m, ax)
                recursive_tree(i, tt, nstep - 1, ax, x, y0 + 1, dx, mass_unit=mass_unit)

def get_idx(tt, halnum, nout):
    """
    Todo: It's better to be included in the tree class.
    """
    ind_now = np.where(tt["NOUT"] == nout)
    tt_sub = tt[ind_now]
    ind_hal = np.where(tt_sub["HALNUM"] == halnum)
    return tt_sub[ind_hal]["IDX"]