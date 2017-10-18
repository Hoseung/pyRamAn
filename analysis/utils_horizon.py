def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx


def nout2nstep(data, nout):
    return data["nstep"][np.where(data["nout"] == nout)]

def nstep2nout(data, nstep):
    try:
        len(nstep)
        import utils.match as mtc
        ind = mtc.match_list_ind(data["nstep"], nstep)
    except:
        ind = np.where(data["nstep"] == nstep)[0]
    return data["nout"][ind]

def zred2nout(data, nout):
    return data["nstep"][np.where(data["nout"] == nout)]



class MainPrg():
    def __init__(self, atree, is_root=False):
        """
            Separate tree data and lambda catalog data. 
            
            early snapshot first.
        """

        self.nsteps = atree["nstep"]
        self.idxs = atree["idx"]
        self.ids = atree["id"]
        self.nouts = nstep2nout(nnza, self.nsteps)
        self.zreds = atree["zred"]
        self.aexps = 1/(1+self.zreds)
        self.root_idx = atree["idx"][-1]
    
    def initialize_data(self, cat, force=False):
        if hasattr(self, "data"):
            if not force:
                print("self.data already exists. use force=True to re-initialize it.")
                pass
        else:
            self.data=np.zeros(len(self.nsteps),
                                  dtype=cat.dtype)
        
    def set_data(self, cat, nout):
        ind = np.where(cat["tree_root_id"] == self.root_idx)[0]
        if len(ind) == 1:
            inout = np.where(self.nouts == nout)[0]
            if len(inout) == 1:
                self.data[inout] = cat[ind]

from tree.hmutils import extract_main_tree, extract_main_tree2
