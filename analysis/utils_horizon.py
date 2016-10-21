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

def extract_main_tree(t, idx, fatherID, fatherMass):
    import numpy as np
    t_now = t[idx]
    nstep = t_now["nstep"]
    nouts = [nstep]
    atree = np.zeros(nstep + 1, dtype=t.dtype)
    atree[0]=t_now
    
    for i in range(1, nstep + 1):
        try:
            #print(t["flist_index"][idx])
            #id_father = fatherID[t["flist_index"][idx]]
            id_father = fatherID[idx]
            #print(id_father)
            #print(len(id_father))
            if len(id_father) > 1:
                #mass_father = fatherMass[t["flist_index"][idx]]
                mass_father = fatherMass[idx]
                #print(mass_father)
                id_father = id_father[np.argmax(mass_father)]
                ind_father = id_father[id_father > 0] -1
                
                nstep -= 1
                t_father = t[np.where(t["nstep"] == nstep)[0]][ind_father]
                idx = t_father["idx"]
                #print(idx)
                atree[i]=t_father
                nouts.append(nstep)
            else:
                break
        except:
            break
    return atree


def extract_main_tree2(idx, t, fatherIDx, fatherMass):
    """
    Two extract_main_tree, which one works? 
    """
    import numpy as np
    t_now = t[idx]
    nstep = t_now["nstep"]
    nouts = [nstep]
    atree = np.zeros(nstep + 1, dtype=t.dtype)
    atree[0]=t_now
    
    for i in range(1, nstep + 1):
        try:
            #id_father = fatherID[t["flist_index"][idx]]
            idx_father = fatherIDx[t["flist_index"][idx]]
            #print(idx)
            if len(idx_father) > 2:
                mass_father = fatherMass[t["flist_index"][idx]]
                #idx_father = idx_father[np.argmax(mass_father)]
                idx = idx_father[np.argmax(mass_father[mass_father > 0])]
            elif len(idx_father) == 2:
                idx = idx_father[idx_father > 0]
            else:
                break
            nstep -= 1
            atree[i]=t[idx]
            nouts.append(nstep)
        except:
            break
    return atree

