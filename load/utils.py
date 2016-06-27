# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 21:20:50 2015

@author: hoseung
"""

import numpy as np
from utils.io import skip_fortran, read_fortran
_head_type = np.dtype('i4')


def read_header(f, dtype, check=True):
    q = np.empty(1, dtype=dtype)
    for i in range(len(dtype.fields)):
        data = read_fortran(f, dtype[i], check=check)
#        print(dtype.names[i], data)

        if np.issubdtype(dtype[i], np.string_):
            q[0][i] = data
        elif hasattr(data[0], "__len__"):
            q[0][i][:] = data[0]
        else:
            q[0][i] = data[0]
        # if string, return the whole array

    return q[0]


def read_header_string(f, dtype):
    alen = np.fromfile(f, _head_type, 1)  # == skip
    length = dtype.itemsize

    if alen != length:
        raise IOError("Unexpected FORTRAN block length %d!=%d"
                      % (alen, length))

    data = np.fromfile(f, dtype, 1)  # Actual data
    alen = np.fromfile(f, _head_type, 1)
    if alen != length:
        raise IOError("Unexpected FORTRAN block length (tail) %d!=%d"
                      % (alen, length))
    return data


# How about reading stars and dms per cpu output?
# No way to know the number of stars in individual output beforehand.
# Read all files once more only to get the size of stars? Maybe no!
# Or, append numpy array hundreds or thousands of times to construct one
# big array? Well...
# For the moment, I'll get the total number of stars and non-stars first.
# Then, allocate momeries of two populations.
#
# Question.
# Which of the three is the best for performance's sake? mask array, view,
# or index array?
#
# +
# Later, add Hilber-domain decomposition consideration. - done


