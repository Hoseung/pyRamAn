# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 21:36:41 2015

Tools to make plots look better


@author: hoseung
"""

def tickInPboxFloat(xr, pboxsize, nticks=3):
    from numpy import arange
    dx = (xr[1] - xr[0]) / nticks
    return ["{:1.2f}".format(i) for i in arange(xr[0], xr[1], dx) * pboxsize]

def tickInPboxInt(xr, npix, nticks=5):
    """
    returns x_org_pos, and x_new,
    where x_org_pos is normalized position of new ticks on data point,
    and x_new is new ticks in physical unit.

    xx, xnew = tickInPboxInt(xr, s.info.pboxsize, 800, nticks=3)
    """
    from numpy import arange
    from numpy import asarray
    from math import floor, ceil
    dx_phy = (xr[1] - xr[0])
    #print(dx_phy)
    # For example, xr = [0.234, 0.345], pboxsize = 262.813870..
    # ntick=5
    # dx_phy = 29.1723396 Mpc
    multiplier, dxpn = toIntegerPart(dx_phy / nticks)
    #print(dxpn)
    # dxpn = 5.83446792...
    dxpn_good = floor(dxpn) * 10**(multiplier)
    #print(dxpn_good)
    # dxpn_good = 5

    x_new_first = ceil(xr[0] / dxpn_good) * dxpn_good
    # 65
    x_new_last = floor(xr[1] / dxpn_good) * dxpn_good
    #print(floor(xr[1] / dxpn_good))
    # 85
    x_new_num = arange(x_new_first, x_new_last * 1.00001, dxpn_good)
    # to include the end point
    x_new = asarray(["{:1.1f}".format(i) for i in x_new_num])
    x_org_pos = (x_new_num - xr[0]) / dx_phy
    return x_org_pos, x_new


def toIntegerPart(val, ndecimal=1):
    """
    Converts number to have 1 digit integer part.
    For example, 123.456 -> 1.23456

    Also returned is the original number of digit.
    , ndigit

    return ndigit, new_val

    To have 2 integer part digits,
    use the option ndecimal = 2

    """
    #from math import ceil
    nID = nIntDigit(val)
    return nID, val* 10**(-1 * nID + ndecimal - 1)

def nIntDigit(val):
    from math import floor, log10
    return floor(log10(val))