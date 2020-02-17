def load_cat(fname):
    import pandas as pd
    import pickle 
    
    cat = pd.DataFrame(pickle.load(open(fname, "rb"))).to_records()
    return cat
    #astype([('index', '<i8'), ('Rgal_to_reff', '<f8'), ('d2t', '<f8'), ('eps', '<f8'), ('epsh', '<f8'), ('epsq', '<f8'), ('id', '<i8'), ('idx', '<i8'), ('lambda_12kpc', 'O'), ('lambda_arr', 'O'), ('lambda_arr2', 'O'), ('lambda_arrh', 'O'), ('lambda_r', '<f8'), ('lambda_r12kpc', '<f8'), #('lambda_r2', '<f8'), ('lambda_rh', '<f8'), ('mgas', '<f8'), ('mrj', '<f8'), ('mstar', '<f8'), ('nstar', '<i8'), ('pa', '<f8'), ('pah', '<f8'), ('paq', '<f8'), ('pq', 'O'), ('pt', 'O'), ('reff', '<f8'), ('reff2', '<f8'), ('rgal', '<f8'), ('rgal2', '<f8'), ('rhalo', '<f8'), ('rscale_lambda', '<f8'), ('sfr', '<f8'), ('sma', '<f8'), ('smah', '<f8'), ('smaq', '<f8'), ('smi', '<f8'), ('smih', '<f8'), ('smiq', '<f8'), ('ssfr', '<f8'), ('star', 'O'), ('vmap', 'O'), ('vx', '<f8'), ('vy', '<f8'), ('vz', '<f8'), ('x', '<f8'), ('xcen', '<f8'), ('xcenh', '<f8'), ('xcenq', '<f8'), ('y', '<f8'), ('ycen', '<f8'), ('ycenh', '<f8'), ('ycenq', '<f8'), ('z', '<f8')])



def remove_nan(x):
    import numpy as np
    return x[~np.isnan(x)]
