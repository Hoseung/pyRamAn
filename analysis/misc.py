def load_cat(fname):
    import pandas as pd
    import pickle 
    return pd.DataFrame(pickle.load(open(fname, "rb"))).to_records()



def remove_nan(x):
    import numpy as np
    return x[~np.isnan(x)]
