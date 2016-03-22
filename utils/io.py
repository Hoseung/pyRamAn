

def read_fortran(f, dtype, n=1):
    if not isinstance(dtype, np.dtype):
        dtype = np.dtype(dtype)

    length = n * dtype.itemsize

    # A highly efficient way of reading binary data with a known data-type
    alen = np.fromfile(f, _head_type, 1) # == skip

    if alen != length:
        raise IOError("Unexpected FORTRAN block length %d!=%d" % (alen, length))

    data = np.fromfile(f, dtype, n) # Actual data

    alen = np.fromfile(f, _head_type, 1)
    if alen != length:
        raise IOError("Unexpected FORTRAN block length (tail) %d!=%d" % (alen, length))

    return data

