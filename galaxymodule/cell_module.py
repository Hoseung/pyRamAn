def angular_momentum_shell(cells, bins = None,
    """
        parameters
        ----------
        bins=None
            If None, total angular momemtum is returned.

        example
        -------
        >>> bins, angs = angular_momentum_shell(gg.cell, bins=[0,2,5,10,20,40])
    """
    #bins = [0,2,5,10,20,40]
    if bins is None:
        # all at once
        vec_rot = np.cross(np.stack((cells["x"],
                                     cells["y"],
                                     cells["z"])).T,
                           np.stack((cells["var1"],
                                     cells["var2"],
                                     cells["var3"])).T)

        return (vec_rot.T * (cell_vec["dx"]**3 * cell_vec["var0"])).sum(axis=1)
    else:
        dd = np.sqrt(np.square(close_cell["x"])+np.square(close_cell["y"])+np.square(close_cell["z"]))
        Lns = []
        if bins[0] !=0:
            print("[cell_module.angular_momentum_shell], Warning... bin doesn't start from 0")
        for i in range(len(bins-1)):
            cell_vec = close_cell[(ds[i] < dd) * (dd < ds[i+1])]
            #print(len(cell_vec), len(close_cell))
            Mean_v = (cells["var1"].mean(),
                      cells["var2"].mean(),
                      cells["var3"].mean())
            vec_rot = np.cross(np.stack((cells["x"],
                                         cells["y"],
                                         cells["z"])).T,
                               np.stack((cells["var1"],
                                         cells["var2"],
                                         cells["var3"])).T)
                                             
            Lns.append((vec_rot.T * (cells["dx"]**3 * cells["var0"])).sum(axis=1))
        return bins, Lns
