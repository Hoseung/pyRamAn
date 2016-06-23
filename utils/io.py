
import numpy as np
import struct
_head_type = np.dtype('i4')

def write_fortran(f, array, check=True):
    """
        Note 
        type 'i' meant to be int32. Can I explicitly tell?
    """
    f.write(struct.pack('i', array.nbytes))
    array.tofile(f)
    f.write(struct.pack('i', array.nbytes))
    
    
def read_fortran(f, dtype, n=1, check=True):
    if not isinstance(dtype, np.dtype):
        dtype = np.dtype(dtype)

    length = n * dtype.itemsize

    alen = np.fromfile(f, _head_type, 1)  # == skip

    if alen != length:
        if check:
            raise IOError("Unexpected FORTRAN block length from"
            "file {} != user given {}".format(alen, length))
        else:
            n = int(alen / dtype.itemsize)
            # Force exact read (although unintended)

    data = np.fromfile(f, dtype, n)  # Actual data

    alen = np.fromfile(f, _head_type, 1)
    if check:
        if alen != length:
            raise IOError("Unexpected FORTRAN block length (tail) %d!=%d"
                          % (alen, length))
    return data

def prettyprint(q, precise=False):
    if isinstance(q, (int, np.int32, np.int64)):
        #print("int")
        return "{:d}".format(q)
    elif abs(q) > 1e4:
        if precise:
            return "{:.6e}".format(q)
        else:
            return "{:.3e}".format(q)
    else:
        if precise:
            return "{:.6f}".format(q)
        else:
            return "{:.2f}".format(q)
