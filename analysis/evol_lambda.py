import numpy as np

class MainPrg():
    def __init__(self, treedata, final_gal, nout_ini=None, nout_fi=None):
        import tree.ctutils as ctu
        temp_tree = ctu.extract_main_tree(treedata, final_gal)
        if nout_ini == None:
            nout_ini = min(temp_tree['nout'])
        if nout_fi == None:
            nout_fi = max(temp_tree['nout'])            
            
        self.nouts = np.arange(nout_fi, nout_ini -1, -1)
        self.idxs = temp_tree['id'] # nout_ini, nout_fi consideration needed.
        self.ids = temp_tree['Orig_halo_id']
        self.data = None
    
    def set_data(self, cat, nout):
        """
        compile data from catalogs.
        """
        if nout in self.nouts:
            # Declare self.data first if there isn't.
            if self.data == None:
                self.data = np.zeros(len(self.nouts), dtype=cat.dtype)
            inow = self.nouts == nout
            a = np.where(cat['idx'] == self.idxs[inow])[0]
            if len(a) > 0:
                self.data[inow] = cat[a]        
            else:
                pass
                #print(self.ids[inow],cat['id'])
        else:
            pass
            #print("No {} in the catalog".format(nout))
            
    def clip_non_detection(self):
        # end of galaxy tree = last non-zero position.
        # Note that 'id' can be 0 if phantom. But phantom is a valid datapoint
        i_first_nout = max(np.where(self.data['idx'] > 0)[0])
        #print('i_first', i_first_nout)
        # then, only [0-i_first_nout] are valid.
        # earlier then 187 - 91-th are zero. so get rid of them.
        self.data = self.data[:i_first_nout].copy()
        self.nouts = self.nouts[:i_first_nout].copy()
        self.ids = self.ids[:i_first_nout].copy()
        self.idxs = self.idxs[:i_first_nout].copy()


    def fill_missing_data(self):
        assert (self.ids[-1] != 0)
        # position angles cannot be linearly interpolated.
        # skip.
        #
        # position and velocity are also not that linear..
        # but let me just interpolate them. 
        #
        # 
        # excluded=["lambda_arr2"]
        filled_fields = ["eps", "epsh", "epsq", "lambda_12kpc",
                         "lambda_arr", "lambda_arrh",
                         "lambda_r","lambda_r12kpc",
                         "lambda_r2","lambda_rh","mgas","mrj","mstar",
                         "reff","reff2","rgal","rgal2","rhalo","rscale_lambda",
                         "sfr","sma","smah","smaq","smi","smih","smiq","ssfr",
                         "vxc", "vyc", "vzc", "xc", "yc", "zc"]

        i_good_max = max(np.where(gal.data["reff"] > 0)[0])
        i_bad = np.where(gal.data['idx'] == 0)[0]
        i_bad = i_bad[i_bad < i_good_max]
        if len(i_bad) > 0:
            for field in filled_fields:
                # do not modify index and id fields.
                arr = gal.data[field] # it's a view.

                for i_b in i_bad:
                    # neighbouring array might also be empty. Search for closest valid element.
                    # left point
                    i_l = i_b - 1
                    while(i_l in i_bad):
                        i_l = i_l - 1
                    # right point
                    i_r = i_b + 1
                    while(i_r in i_bad):
                        i_r = i_r + 1
                    arr[i_b] = (arr[i_b -1] + arr[i_b +1])/2.



# In[2]:

def fixed_ind_Lr(gal):
    nnouts = len(gal.nouts)
    ind_reff_fix = np.zeros(nnouts, dtype='i4')

    #print(gal.data['rgal'])
    smooth_r = smooth(mask_outlier(gal.data['rgal'], 1.5, 1.5), 50, monotonic=False)

    # fixed Reff array
    for i in range(nnouts):
        # 1Reff = 5 points
        reff_real = smooth_r[i]
        reff = gal.data['rgal'][i]
        try:
            ind_reff_fix[i] = np.round(reff_real/reff * 5) -1
        except:
            pass
    return ind_reff_fix


def smoothed_reff(cat, nout_merger):
    """
    returns "representative" lambda at each nout by assuming monotonic change in Reff. 
    During merger, Reff can fluctuate, and if has no physical meaning to infer Labda at Reff during merger stage. 
    So Reff' is derived by linear interpolating Reffs before and after the merger. 
    
    cat is one galaxy catalog over time.
    """
    import utils.match as mtc
    i_merger = np.where(cat['nout'] == nout_merger)[0]
    ind_lower = 20
    ind_upper = 20
    
    reffs = cat['rgal']
    # left and right values chosen by sigma-clipping
    r_lefts, b, c = scipy.stats.sigmaclip(reffs[max([0,i_merger-ind_lower]):i_merger], sig_lower, sig_upper)
    #print(r_lefts)
    r_left = r_lefts[-1]
    i_left = np.where(reffs == r_left)[0]

    r_rights, b,c = scipy.stats.sigmaclip(reffs[i_merger:min([i_merger+ind_upper,len(reffs)])], sig_lower, sig_upper)
    r_right = r_rights[0]
    i_right = np.where(reffs == r_right)[0]

    r_prime = reffs
    #print("chekc")
    #print(r_prime)
    r_prime[i_left : i_right + 1] = np.linspace(r_left, r_right, i_right - i_left + 1)
    return r_prime    

