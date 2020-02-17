# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 17:03:40 2015

@author: hoseung

Routines borrowed from web are saved here.

"""

import numpy as np
from itertools import chain


def find(a, predicate, chunk_size=1024):
    """
    Find the indices of array elements that match the predicate.

    Parameters
    ----------
    a : array_like
        Input data, must be 1D.

    predicate : function
        A function which operates on sections of the given array, returning
        element-wise True or False for each data value.

    chunk_size : integer
        The length of the chunks to use when searching for matching indices.
        For high probability predicates, a smaller number will make this
        function quicker, similarly choose a larger number for low
        probabilities.

    Returns
    -------
    index_generator : generator
        A generator of (indices, data value) tuples which make the predicate
        True.

    See Also
    --------
    where, nonzero

    Notes
    -----
    This function is best used for finding the first, or first few, data values
    which match the predicate.

    Examples
    --------
    >>> a = np.sin(np.linspace(0, np.pi, 200))
    >>> result = find(a, lambda arr: arr > 0.9)
    >>> next(result)
    ((71, ), 0.900479032457)
    >>> np.where(a > 0.9)[0][0]
    71


    """
    if a.ndim != 1:
        raise ValueError('The array must be 1D, not {}.'.format(a.ndim))

    i0 = 0
    chunk_inds = chain(range(chunk_size, a.size, chunk_size),
                 [None])

    for i1 in chunk_inds:
        chunk = a[i0:i1]
        for inds in zip(*predicate(chunk).nonzero()):
            yield (inds[0] + i0, ), chunk[inds]
        i0 = i1
