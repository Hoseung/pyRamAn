import numpy as np
cimport numpy as np

DTYPE = np.float
ctypedef np.float_t DTYPE_t

def col_over_denom(np.ndarray[np.int64_t, ndim=1] iin, \
     np.ndarray[np.float64_t, ndim=1] ixl, \
     np.ndarray[np.float64_t, ndim=1] ixr, \
     np.ndarray[np.float64_t, ndim=1] iyl, \
     np.ndarray[np.float64_t, ndim=1] iyr, \
     np.ndarray[np.float64_t, ndim=1] mass, \
     np.ndarray[np.float64_t, ndim=1] sden, \
     int nx, int ny, int column):

    cdef int percent = 0
    colden = np.zeros((nx,ny), dtype=np.float32)
    denom =  np.zeros((nx,ny), dtype=np.float32)

    if column:
        for i in iin:
#            i1 = ixl[i]
#            i2 = ixr[i]
#            j1 = iyl[i]
#            j2 = iyr[i]
#            colden[i1:i2+1, j1:j2+1] += sden[i]
            colden[ixl[i]:ixr[i]+1, iyl[i]:iyr[i]+1] += sden[i]
        return colden
    else:
        for i in iin:
            i1 = ixl[i]
            i2 = ixr[i]
            j1 = iyl[i]
            j2 = iyr[i]
            # no smoothing or max, no column
            colden[i1:i2+1, j1:j2+1] += sden[i]
            denom[i1:i2+1, j1:j2+1] += mass[i]
        return colden / denom

