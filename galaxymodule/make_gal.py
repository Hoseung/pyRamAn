import numpy as np

def extract_data_old(halo, rscale=1.25):
    xc_tmp0 = halo['x']
    yc_tmp0 = halo['y']
    zc_tmp0 = halo['z']

    rr_tmp0 = min([halo['r'] * rscale, 0.0002])
    # it's galaxy, not halo. 
    # 'r' is the furthest stellar particle
    # 'rvir' doesn't have robust physical meaning 
    # because galaxies are hardly virialized systems.

    # arbitrary! < 20kpc
    rr_tmp0 = max([rr_tmp0, 0.000025])
    # When merger occurs, larger radius is likely to include 
    # companion galaxy resulting center to be in the middle of nowhere.
    # If you want a larger galaxy, # increase rgal_tmp instead. 
    #        
    # xx is easier to search for than x.

    if star_all is not None:
        ind_s = np.where((star_all['x'] - xc_tmp0)**2 + (star_all['y'] - yc_tmp0)**2
                        + (star_all['z'] - zc_tmp0)**2 < rr_tmp0**2)[0]
    if dm_all is not None:
        ind_d = np.where((dm_all['x'] - xc_tmp0)**2 + (dm_all['y'] - yc_tmp0)**2
                        + (dm_all['z'] - zc_tmp0)**2 < rr_tmp0**2)[0]
    if cell_all is not None:
        ind_c = np.where((cell_all['x'] - xc_tmp0)**2 + (cell_all['y'] - yc_tmp0)**2
                        + (cell_all['z'] - zc_tmp0)**2 < rr_tmp0**2)[0]
    else:
        return star_all[ind_s], dm_all[ind_d], None

#    print(len(ind_s), len(ind_d), len(ind_c))    

    return star_all[ind_s], dm_all[ind_d], cell_all[ind_c]

def extract_data(halo, rscale=1.25):
    xc_tmp0 = halo['x']
    yc_tmp0 = halo['y']
    zc_tmp0 = halo['z']

    rr_tmp0 = min([halo['r'] * rscale, 0.0002])
    # it's galaxy, not halo. 
    # 'r' is the furthest stellar particle
    # 'rvir' doesn't have robust physical meaning 
    # because galaxies are hardly virialized systems.

    # arbitrary! < 20kpc
    rr_tmp0 = max([rr_tmp0, 0.000025])
    # When merger occurs, larger radius is likely to include 
    # companion galaxy resulting center to be in the middle of nowhere.
    print(cell_all, "cell_all")
    if cell_all is not None:
        ind_c = np.where((cell_all['x'] - xc_tmp0)**2 + (cell_all['y'] - yc_tmp0)**2
                        + (cell_all['z'] - zc_tmp0)**2 < rr_tmp0**2)[0]
    else:
        return None

    return cell_all[ind_c]




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
            if not isinstance(attr, (np.ndarray, np.recarray, dict, list, str))}



def save_gal(galaxy, filename):
    import h5py as hdf

    outfile = hdf.File(filename, 'w')

    # Store metadata using HDF5 attributes
    attrs = get_metadata(galaxy.meta)
    attrs = get_metadata2(attrs)
    #outfile.attrs.create("all", attrs)
    for name, atr in attrs.items():
#        if name == "meta":
#            for name2, atr2 in atr.__dict__.items():
#                outfile.attrs.create(name2, atr2)
#        else:
        if atr != None:
            outfile.attrs.create(name, atr)
        #outfile.attrs[name] = atr

    # Store data under /galaxy with direct assignment
    if hasattr(galaxy, 'star'):
        print("Saving star")
        star = outfile.create_group("star")
        for field in galaxy.star.dtype.names:
            star.create_dataset(field, data=galaxy.star[field])

    if hasattr(galaxy, 'dm'):
        print("Saving DM")
        dm = outfile.create_group("dm")
        for field in galaxy.dm.dtype.names:
            dm.create_dataset(field, data=galaxy.dm[field])

    if hasattr(galaxy, 'cell'):
        print("Saving cell")
        gas = outfile.create_group("cell")
        for field in galaxy.cell.dtype.names:
            gas.create_dataset(field, data=galaxy.cell[field])

    outfile.close()





