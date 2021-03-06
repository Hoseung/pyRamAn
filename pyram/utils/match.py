# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 03:08:24 2015

@author: hoseung
"""
import numpy as np

def closest(var, arr1, index=False):
    """
    Returns index of arr1 that has nearest value to var.
    ex) arr = [0, 1, 2, 3, 4], var=0.8
    > nearest(var, arr)
    > 1
    """

    return (np.abs(arr1 - var)).argmin()

def match_list_ind(arr1, arr2,
                   allow_swap=True,
                   verbose=True,
                   repeated_elements=False):
    '''
    Returns indices of matching elements in arr1.

    Parameters
    ----------
    arr1 : list or array
        longer array
    arr2 : list or array
        shorter array
    allow_swap : boolean
        If array B is longer than the array A,
        then returned indices are of array B.
        If allow_swap is False, and the array A is shorter,
        then only common elements are considered.
    repeated_elements : boolean
        False by default, and no repeated elements in arr2 arr assumed.
        If True, indices of repeated elements in arr2 are returned.
        In this case, comparing lengths of two arrays makes no sesne,
        So swaping is disabled. User must make sure that
        arr1 is the target array.

    Example,
    >>> a = [0,1,2,3,4,5,6]
    >>> b = [3,4,5,6,7]
    >>> print(a[match_list_ind(a, b)])
    >>> array([3,4,5,6])

    >>> a = [0,1,2,3,4,5,6]
    >>> b = [3,4,5,6,7]
    >>> print(a[match_list_ind(b, a, allow_swap=True)])
    >>> array([3,4,5,6])


    Supports both list and array.
    To do : Can I support any types of collections.sequence?
        does collections.sequence include types without order?


    '''

    # If arr1 is not a sequence,  return
    try:
        if len(arr1) == 0:
            return
        elif len(arr1) == 1:
            return np.where(arr1 == arr2)[0]
    except:
        raise ValueError("Array is needed")

    if len(arr2) == 1:
        return np.where(arr1 == arr2)[0]

    # If list, convert to array
    if isinstance(arr1, list):
        arr1 = np.array(arr1)
    if isinstance(arr2, list):
        arr1 = np.array(arr2)

    if len(arr1) > len(arr2) or repeated_elements:
        bigArr = arr1
        smallArr = arr2
    else:
        if allow_swap:
            bigArr = arr2
            smallArr = arr1
        else:
            bigArr = arr1
            smallArr = np.intersect1d(arr1, arr2)

    # sort big array so that we can use bisection method, which is fast.
    sortedind = np.argsort(bigArr)
    sortedbigArr = bigArr[sortedind]
    sorted_index = np.searchsorted(sortedbigArr, smallArr)
    smallindex = np.take(sortedind, sorted_index, mode="clip")
    mask = bigArr[smallindex] != smallArr

    return np.ma.array(smallindex, mask=mask).compressed()
