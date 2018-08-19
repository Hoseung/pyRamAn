def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx


def nouts_from_zreds(zreds, base='./'):
#if True:
    """
    Look for info in wdir/snapshots/outpu*/
    """
    import glob
    from load.info import Info
    finfos = glob.glob(wdir + "snapshots/output*/info*")
    nouts_info = []
    zreds_info = []
    for fn in finfos:
        ii = Info(fn=fn)
        nouts_info.append(ii.nout)
        zreds_info.append(ii.zred)

    nouts_info = np.array(nouts_info)
    zreds_info = np.array(zreds_info)

    isort = np.argsort(zreds_info)
    nouts_info = nouts_info[isort]
    zreds_info = zreds_info[isort]

    nouts=[]
    for zz in zreds:
        nouts.append(nouts_info[find_closest(zreds_info, zz)])

    return nouts

