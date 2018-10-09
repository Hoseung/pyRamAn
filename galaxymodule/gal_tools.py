import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from draw import pp
from matplotlib.colors import LogNorm


def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


def _close_to(values, v_ref):
    """
        returns the index of the value which is close to the v_ref.
    """
    ind = 0
    good=values[ind]
    for i, vv in enumerate(values):
        if abs(vv - v_ref) < good:
            good = vv
            ind = i

    return ind


def save_gal(self, base='./'):

    def get_metadata(clazz):
        """
            Out of all attributes of a galaxy instance, leave only data.
        """
        return {name: attr for name, attr in clazz.__dict__.items()
                if not name.startswith("__")
                and not callable(attr)
                and not type(attr) is staticmethod}

    def get_metadata2(adict):
        return {name: attr for name, attr in adict.items()
                if not isinstance(attr, (np.ndarray, np.recarray, dict, list))}


    import h5py as hdf
    # Save data into a hdf5 file
    outfile = hdf.File(base + str(self.meta.id).zfill(6) + '_gal.hdf5',
                       'w', libver='latest')

    # Store metadata in HDF5 attributes
    attrs = get_metadata(self)
    attrs = get_metadata2(attrs)
    for name, atr in attrs.items():
        if atr != None:
            outfile.attrs.create(name, atr)
        #outfile.attrs[name] = atr

    # Store data under /selfaxy with direct assignment
    if hasattr(self, 'star'):
        #print("Saving star")
        star = outfile.create_group("star")
        for field in self.star.dtype.names:
            star.create_dataset(field, data=self.star[field])

    if hasattr(self, 'dm'):
        #print("Saving DM")
        dm = outfile.create_group("dm")
        for field in self.dm.dtype.names:
            dm.create_dataset(field, data=self.dm[field])

    if hasattr(self, 'cell'):
        #print("Saving gas")
        gas = outfile.create_group("gas")
        for field in self.cell.dtype.names:
            gas.create_dataset(field, data=self.cell[field])

    outfile.close()
