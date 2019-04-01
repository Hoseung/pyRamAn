import numpy as np
import struct
import tipsy

def load_iord(fn, swap=False):
    if swap:
        fmt = "<i"
    else:
        fmt = ">i"

    with open(fn, "rb") as f:
        #f.seek(0)
        n_id = struct.unpack(fmt, f.read(4))[0]
        return np.frombuffer(f.read(np.dtype(np.int32).itemsize * n_id),
                            dtype=fmt)




def load_partial(fn, ind_to_load, lmax=1e5, swap=False):
    """
        Load only indexed particles from the raw data IF
        the length of ind is shorter than lmax.
        Otherwise, read all at once and extract relevant items.
    """
    if len(ind_to_load) > lmax:
        pass
        """
        call load.load_data()
        """

    else:
        from tipsy.load import dtype_gas, dtype_dm, dtype_star

        #ind_to_load = np.arange(1, 100, 10)

        f = open(fn, "rb")
        header = tipsy.load.TipsyHeader(f)

        iord_dm_first = iords[header.nsph]
        iord_star_first = iords[header.nsph+header.ndark]

        gind_to_load = ind_to_load[ind_to_load < iord_dm_first]
        dind_to_load = ind_to_load[(ind_to_load > iord_dm_first) * (ind_to_load < iord_star_first)]
        sind_to_load = ind_to_load[(ind_to_load > iord_star_first) ]

        g_partial = np.zeros(len(gind_to_load), dtype=dtype_gas)
        d_partial = np.zeros(len(dind_to_load), dtype=dtype_dm)
        s_partial = np.zeros(len(sind_to_load), dtype=dtype_star)

        off_header = f.tell()

        for i, ind in enumerate(gind_to_load):
            f.seek(off_header + ind*dtype_gas.itemsize)
            g_partial[i] = np.frombuffer(f.read(dtype_gas.itemsize), dtype=dtype_gas)

        off_gas = off_header + dtype_gas.itemsize * header.nsph

        for i, ind in enumerate(dind_to_load):
            f.seek(off_gas + (ind-header.nsph)*dtype_dm.itemsize)
            d_partial[i] = np.frombuffer(f.read(dtype_dm.itemsize), dtype=dtype_dm)

        off_dm = off_header + dtype_gas.itemsize * header.nsph + dtype_dm.itemsize * header.ndark
        for i, ind in enumerate(sind_to_load):
            f.seek(off_dm + (ind-header.nsph-header.ndark)*dtype_star.itemsize)
            s_partial[i] = np.frombuffer(f.read(dtype_star.itemsize), dtype=dtype_star)

        f.close()
        if swap:
            g_partial.newbyteorder()
            d_partial.newbyteorder()
            s_partial.newbyteorder()

        return g_partial, d_partial, s_partial
