import struct
import numpy as np
_head_type = np.dtype('i4')

def write_fortran(f, array, dtype="I", check=True):
    """
        Writes an array of type dtype to a fortran file f.

        Parameters
        ----------
        f :
            fortran file stream
        array :
            numpy array
        dtype : 
            data type code
        check :
            no implemented. - what was the intention??

        Note
        ----
        dtype is for struct.pack.

        x	pad byte	no value
        c	char	bytes of length 1	1
        b	signed char	integer	1	(1),(3)
        B	unsigned char	integer	1	(3)
        ?	_Bool	bool	1	(1)
        h	short	integer	2	(3)
        H	unsigned short	integer	2	(3)
        i	int	integer	4	(3)
        I	unsigned int	integer	4	(3)
        l	long	integer	4	(3)
        L	unsigned long	integer	4	(3)
        q	long long	integer	8	(2), (3)
        Q	unsigned long long	integer	8	(2), (3)
        n	ssize_t	integer	 	(4)
        N	size_t	integer	 	(4)
        e	(7)	float	2	(5)
        f	float	float	4	(5)
        d	double	float	8	(5)
        s	char[]	bytes
        p	char[]	bytes
        P   void *

    """
    f.write(struct.pack(dtype, array.nbytes))
    array.tofile(f)
    f.write(struct.pack(dtype, array.nbytes))

def skip_fortran(f, n=1, verbose=False):
    """
        skip one record blod of the given fortran file.
    """
    alen = np.fromfile(f, _head_type, 1)
    if verbose:
        print("FORTRAN block length %d!=%d" % (alen))

    # Check
    if alen % 4 != 0:
        print("Array size is not a multiple of 4")

    n = int(alen/4)
    np.fromfile(f, _head_type, n)
    np.fromfile(f, _head_type, 1)

def read_fortran(f, dtype, n=1, check=True):
    """
    read one record block from a fortran file.

    parameters
    ----------
    f :
        fortran file stream
    dtype :
        numpy dtype
    n :
        number of data elements in the record block 
    """
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

def read_header(f, dtypes, check=True):
    """
    Read multiple SINGLE entries from a fortran file
    
    parameters
    ----------
    f :
        fortran stream
    dtypes :
        list of dtypes of header quantities

    Todo
    ----
    Handle string 
    """
    if not isinstance(dtypes, np.dtype):
        dtypes = np.dtype(dtypes)

    q = np.empty(1, dtype=dtypes)
    for i in range(len(dtypes.fields)):
        data = read_fortran(f, dtypes[i], check=check)

        if np.issubdtype(dtypes[i], np.string_):
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
