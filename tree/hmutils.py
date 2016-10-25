def extract_main_tree(t, idx, fatherID, fatherMass):
    """
    Extracts 
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
            id_father = fatherID[idx]
            if len(id_father) > 1:
                #mass_father = fatherMass[t["flist_index"][idx]]
                mass_father = fatherMass[idx]
                #print(mass_father)
                id_father = id_father[np.argmax(mass_father)]
                ind_father = id_father[id_father > 0] -1

                nstep -= 1
                t_father = t[np.where(t["nstep"] == nstep)[0]][ind_father]
                idx = t_father["idx"]
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

