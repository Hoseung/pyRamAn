def mk_gal(halodata, out_q, info, i,
            save=False, rscale=0.3, verbose=True):
#            , def_param=(dm, star, cell)):
    from galaxymodule import galaxy

#    print(halodata['id'])
    print("{}-th **********************************************".format(i))

    gal = galaxy.Galaxy(halodata, radius_method='simple', info=info)
#    print(gal.halo)
    gal.mk_gal_mp(star=star, dm=dm, cell=cell,
               rscale=rscale, verbose=True)
    gal_out = {"id":0, "xc":0.0, "yc":0.0, "zc":0.0,
               "vx":0.0, "vy":0.0, "vz":0.0,
               "mstar":0.0, "nstar":0.0, "mgas":0.0,
               "lambda_arr":[], "lambda_r":0, "rgal":0}

    if gal.star is False:
        print(gal.id, " Not a good galaxy")
        out_q.put(gal_out)
    else:
        print("R eff:", gal.reff * info.pboxsize)
        gal.cal_lambda_r(npix=npix_lambda, method=1, rscale=rscale_lambda) # calculate within 1.0 * reff    
#        gal.plot_gal(save_dir= galaxy_plot_dir, ioff=True)
    #            gal.save_gal(base=wdir)

        # Instead of galaxy class, save them in a dict.    
        gal_out['mstar'] = gal.mstar
        gal_out['mgas'] = gal.mgas
        gal_out['nstar'] = gal.nstar
        gal_out['id'] = gal.id
        gal_out['xc'] = gal.xc * info.pboxsize
        gal_out['yc'] = gal.yc * info.pboxsize
        gal_out['zc'] = gal.zc * info.pboxsize
        gal_out['vx'] = gal.vxc
        gal_out['vy'] = gal.vyc
        gal_out['vz'] = gal.vzc        
        gal_out['lambda_arr'] = gal.lambda_arr
        gal_out['lambda_r'] = gal.lambda_r
        gal_out['rgal'] = gal.reff * info.pboxsize * 1000.0 # in kpc
        out_q.put(gal_out)
    print("{}-th ********************** DONE *************".format(i))
    return
    
def mk_gal3(halodata, out_q, info, i,
            save=False, rscale=0.3, verbose=True):
    print("{}-th ********************** DONE *************".format(i))
    return

def static_array2d(shape):
    import multiprocessing
    import ctypes
    import numpy as np
    shared_array_base = multiprocessing.Array(ctypes.c_double, shape[0]*shape[1], lock=False)
    # Lock = False because no process is going to modify the array. 
    # So it is safe even multiple processes access the array at the same time.
    # (But note that speed up due to simulateneous access is still limited by 
    # total memory bandwidth.)
    # Some handles to the memory? 
    shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
    # reshape => add convenience index
    return shared_array.reshape(shape[0], shape[1])

#%%

if __name__ == '__main__':
    """
    The processing pool needs to be instantiated in the main 
    thread of execution. 
    """
    
    import load
    import tree
    import numpy as np
    import utils.sampling as smp
    
    #ncore = int(input("How man cores? \n"))
    #wdir = input("Working directory \n")
    #nout = int(input("nout? \n"))
    #wdir = './'
    wdir = '/home/hoseung/Work/data/AGN2/'
    ncore = 3
    nout=132
    snout = str(nout)
    rscale = 0.8
    npix=800
    rscale_lambda = 3.0
    npix_lambda = 30
    
    m_halo_min = 1e10 # minimum halo mass above which galaxies are searched for. 
    
    fcat = wdir +"catalog" + snout + ".pickle"
    galaxy_plot_dir = wdir + 'galaxy_plot' + snout
    
    
    info = load.info.Info()
    info.setup(nout=nout, base=wdir)
    info.load()
    
    frefine= 'refine_params.txt'
    fnml = 'cosmo_200.nml'
    
    ptypes=["star id pos mass vel time metal", "dm id pos mass vel"]
    
    # Load all halo
    #hall = tree.halomodule.Halo(nout=nout, base=wdir, halofinder="RS", info=info)
    hall = tree.halomodule.Halo(nout=nout, base=wdir, halofinder="HM", info=info)
    hall.load()
    
    # subset of halos ONLY inside zoom-in region
    i_center = np.where(hall.data['np'] == max(hall.data['np']))[0]
    h_ind = smp.extract_halos_within(hall, i_center, scale=2.0)
    h = tree.halomodule.Halo()
    h.derive_from(hall, h_ind)
    #h.derive_from(hall, h_ind)#,  [4921, 5281, 5343, 5365, 5375, 5412, 5415], 5442, 5639, 5665, 6095])
    
    region = smp.set_region_multi(h.data.x, yc=h.data.y, zc=h.data.z, radius = h.data.rvir * rscale)
    print(region)
    
    hind = np.where((hall.data.x > region["xr"][0]) & (hall.data.x < region["xr"][1]) &
                    (hall.data.y > region["yr"][0]) & (hall.data.y < region["yr"][1]) &
                    (hall.data.z > region["zr"][0]) & (hall.data.z < region["zr"][1]) &
                    (hall.data.mvir > m_halo_min))[0]
    h.derive_from(hall, hind[6:30])
    region = smp.set_region_multi(h.data.x, yc=h.data.y, zc=h.data.z, radius = h.data.rvir * rscale)
    
    #%%
    s = load.sim.Sim()
    s.setup(nout, wdir)
    
    s.set_ranges(region["ranges"])
    
    s.show_cpus()
    
    s.add_part(ptypes)
    
    s.add_hydro()
    #%%
   
    # I assume there is no need to free these shared arrays.
    from load import part_shared
    
    # Set ranges
    xmi = s.info.ranges[0][0]
    xma = s.info.ranges[0][1]
    ymi = s.info.ranges[1][0]
    yma = s.info.ranges[1][1]
    zmi = s.info.ranges[2][0]
    zma = s.info.ranges[2][1]
    work_dir = s.info.base + '/snapshots/output_' + str(s.info.nout).zfill(5)
    
    print("Check 1")
    print(work_dir)
    
    #%%
    ndm_actual, nstar_actual, nsink_actual = part_shared.count_part(work_dir, xmi, xma, ymi, yma, zmi, zma)
    
    print(nstar_actual, ndm_actual, nsink_actual)
    #%%
    #star = static_array2d((nstar_actual, 10))
    #dm = static_array2d((ndm_actual, 8))
    #sink=static_array2d((8, nsink_actual))
    
    star_tmp, dm_tmp = part_shared.load_part(
                     nstar_actual, ndm_actual, nsink_actual,
                     work_dir, xmi, xma, ymi, yma, zmi, zma)
    
    #%%
    
    star = static_array2d((10,nstar_actual))
    dm = static_array2d((8,ndm_actual + nsink_actual)) # no distin
    
    # pos 3, vel 3, mass, id, time, metal
    # convert to km/s
    
    star[0] = star_tmp[:,0]
    star[1] = star_tmp[:,1]
    star[2] = star_tmp[:,2]
    star[3] = star_tmp[:,3] * s.info.kms
    star[4] = star_tmp[:,4] * s.info.kms
    star[5] = star_tmp[:,5] * s.info.kms
    star[6] = star_tmp[:,6] * s.info.msun
    star[7] = star_tmp[:,7]
    star[8] = star_tmp[:,8]
    star[9] = star_tmp[:,9]
    
    if 'dm' in s.part.pt:
        dm[0] = dm_tmp[:,0]
        dm[1] = dm_tmp[:,1]
        dm[2] = dm_tmp[:,2]
        dm[3] = dm_tmp[:,3] * s.info.kms
        dm[4] = dm_tmp[:,4] * s.info.kms
        dm[5] = dm_tmp[:,5] * s.info.kms
        dm[6] = dm_tmp[:,6] * s.info.msun
        dm[7] = dm_tmp[:,7]
    
    
    s.hydro.amr2cell(lmax=19)
    # hydro.amr2cell is writte in fortran, but does not accept 
    # shared (row-major) array. 
    # Currently, 
    # 1) read by fortran into fortran array
    # 2) COPY to numpy arrays
    # 3) Copy again to shared arrays 
    # Of course, step 2) could be excluded. 
    # But it will be impossible to keep it backward compatible. 
    # Scripts using particles must be re-written accordingly. (2015-07-25)
    
    #ngridtot, work_dir, xmi, xma, ymi, yma, zmi, zma, lmax = s.hydro.amr2cell(lmax=19, return_meta=True)
    cell = static_array2d((10,int(s.hydro.cell.shape[0])))
    
    cell[0] = s.hydro.cell['x']
    cell[1] = s.hydro.cell['y']
    cell[2] = s.hydro.cell['z']
    cell[3] = s.hydro.cell['dx']
    cell[4] = s.hydro.cell['var0']
    cell[5] = s.hydro.cell['var1']
    cell[6] = s.hydro.cell['var2']
    cell[7] = s.hydro.cell['var3']
    cell[8] = s.hydro.cell['var4']
    cell[9] = s.hydro.cell['var5']
    # Why cells are in physical unit from the debining whereas part is in code unit? 
    
    
   
    
    
    
    import time
    t1 = time.time()

    import multiprocessing as mp
    m = mp.Manager()
    out_q = m.Queue()

    nh = len(h.data)

    pool = mp.Pool(processes=ncore)
    for i in range(nh):
        print(i,"-th started")
        pool.apply(mk_gal, args=(h.data[i], out_q, s.info, i))

    time.sleep(10)
    print("----------Done---------")
    print(" Took", time.time() - t1, "Seconds")

    dictout=[]
    try:
        f = open(wdir + 'galaxies' + snout + '.txt', 'w')
    except:
        print("No filename is given.\n ")
    
    f.write(" #      ID        x          y       z[Mpc]       vx      vy     vz[km/s]")
    f.write("    Reff[kpc]      Mstar    Mgas[Msun]\n")    
    
    for i in range(nh):
#        print("dict out", i)
        dd = out_q.get(timeout=2)
        f.write("{:<4}   {:<4}   {:.5f}  {:.5f}  {:.5f}".format(i,
                dd['id'],dd['xc'],dd['yc'],dd['zc']))
        f.write("  {:.3f}  {:.3f}  {:.3f}".format(dd['vx'],dd['vy'],dd['vz']))
        f.write("  {:.6f}  {:.0f}  {:.0f}     \n".format(dd['rgal'],dd['mstar'], dd['mgas']))
        dictout.append(dd)

    f.close()    
    print("Text file written")
    import gc
    gc.collect()
    """
    import pandas as pd
    catalog = pd.DataFrame(dictout).to_records()    
    
    import pickle
    with open(fcat, 'wb') as f:
        pickle.dump(catalog, f)
    """

