def mk_gal(halodata, out_q, info, i,
            save=False, rscale=0.3, verbose=True):
#            , def_param=(dm, star, cell)):
    from galaxy import galaxy

#    print(halodata['id'])
    print("{}-th **********************************************".format(i))
    print("ID: ", id(dm))
    print("ID: ", id(star))

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
    
def static_array2d(shape):
    import multiprocessing
    import ctypes
    import numpy as np
    shared_array_base = multiprocessing.Array(ctypes.c_double, shape[0]*shape[1], lock=False)
#    shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
#    shared_array_base = mp.RawArray(ctypes.c_double, shape[0]*shape[1])

    # Some handles to the memory? 
    shared_array = np.ctypeslib.as_array(shared_array_base)
    # reshape => add convenience index
    return shared_array.reshape(shape[0], shape[1])
    
    
def extract_halos_within(halos, ind_center, scale=1.0, Mcut=1e5):
    import numpy as np
    import utils.sampling as smp
    '''
    Returns indices of halos within SCALE * Rvir of the central halo.

    def extract_halos_within(halos, ind_center, scale=1.0)
    halos : halo finder output (single snapshot)
    ind_center : index of central halo
    scale : multiplying factor to the Rvir of the central halo
    '''
    xc = halos['x'][i_center]
    yc = halos['y'][i_center]
    zc = halos['z'][i_center]
    rvir= halos['rvir'][i_center]

    xx = halos['x']
    yy = halos['y']
    zz = halos['z']
    m = np.array(halos['m'])

    dd = smp.distance_to([xc,yc,zc],[xx,yy,zz])

    i_m = m > Mcut
    i_ok = np.logical_and(dd < (rvir * scale), i_m)

    return i_ok
#%%

"""
The processing pool needs to be instantiated in the main 
thread of execution. 
"""
import multiprocessing as mp    
import load
from tree import tmtree
import numpy as np
import utils.sampling as smp
#import ctypes
import tree.halomodule as hmo 
from utils import match
#ncore = int(input("How man cores? \n"))
#wdir = input("Working directory \n")
#nout = int(input("nout? \n"))
wdir = './'
#wdir = '/home/hoseung/Work/data/AGN2/'
ncore = 2

rscale = 0.8
npix=800
rscale_lambda = 3.0
npix_lambda = 30
lmax = 19

## halo part ###

m_halo_min = 1e10 # minimum halo mass above which galaxies are searched for. 
ncore = 4
fixed_position = True
nout_ini = 37
# 27 : z=4
# 37 : z=3
# 20 : ~ z=5
nout_fi = 132
rvir=3.0
# nout_halo = 122 == nout 10
# nout_halo = 0   == nout 132
nouts = range(nout_fi, nout_ini -1, -1) 
Nnouts = len(nouts)

tt = tmtree.load(work_dir=wdir, filename="halo/TMtree.fits")
tfin = tt[np.where(tt['NOUT'] == 0)]
tini = tt[np.where(tt['NOUT'] == nout_fi - nout_ini)]
#%%
info = load.info.Info()
info.setup(nout=nout_fi, base=wdir)
info.load()
hh = hmo.Halo(base=wdir, nout=nout_fi, halofinder='HM', info=info)
hh.load()
halo = hh.data
#halo = hmu.load_data(nout_fi, work_dir=work_dir, normalize=True)
i_center = np.where(halo['np'] == max(halo['np']))
i_satellites = extract_halos_within(halo, i_center, scale=3.0)[0]
print("Total {0} halos \n{1} halos are selected".format(
      len(i_satellites),sum(i_satellites)))
#%%
# halos found inside the cluster and has tree back to nout_ini
halo_list = halo['id'][i_satellites]
#print(halo_list)
h_ind_ok, halo_ok = tmtree.check_tree_complete(tt, 0, nout_fi - nout_ini, halo_list)
print(len(halo_ok), "halos left")


for inout, nout in enumerate(nouts):
    snout = str(nout)

    print("###### ", nout, " ######\n")    
    fcat = wdir +"catalog" + snout + ".pickle"
    galaxy_plot_dir = wdir + 'galaxy_plot' + snout
    
    info = load.info.Info()
    info.setup(nout=nout, base=wdir)
    info.load()
    
    frefine= 'refine_params.txt'
    fnml = 'cosmo_200.nml'
    
    ptypes=["star id pos mass vel time metal", "dm id pos mass vel"]
    
    # Load all halo
    hh = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info)
    hh.load()  

    hind = match.match_list_ind(hh.data['id'], halo_ok[5:6,inout])
    
    h = hmo.Halo(base=wdir, nout=nout, halofinder='HM', info=info)
    h.derive_from(hh, hind)
    region = smp.set_region_multi(h.data.x, yc=h.data.y, zc=h.data.z, radius = h.data.rvir * rscale)
    #print(region)
    #%%
    s = load.sim.Sim()
    s.setup(nout, wdir)    
    s.set_ranges(region["ranges"])
    s.show_cpus()
    s.add_part(ptypes)    
    s.add_hydro()
   
    # I assume there is no need to free these shared arrays.
#    from load import part_shared
    
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
    """
    ndm_actual, nstar_actual, nsink_actual = part_shared.count_part(work_dir, xmi, xma, ymi, yma, zmi, zma)    
    print(nstar_actual, ndm_actual, nsink_actual)   
    star_tmp, dm_tmp = part_shared.load_part(
                     nstar_actual, ndm_actual, nsink_actual,
                     work_dir, xmi, xma, ymi, yma, zmi, zma)
    """
    s.part.load(fortran=True)
    print("part load")

    nstar_actual = int(s.part.nstar)
    ndm_actual = int(s.part.ndm)
    nsink_actual = int(s.part.nsink)
        
    star = static_array2d((10,nstar_actual))
    print("star allocate")
    print(nstar_actual, star.shape)
    dm = static_array2d((8,ndm_actual)) # no distin
    print("part allocate")
    print(ndm_actual, dm.shape)    
 #%%   
    # pos 3, vel 3, mass, id, time, metal
    # convert to km/s
    
    star[0] = s.part.star['x']
    star[1] = s.part.star['y']
    star[2] = s.part.star['z']
    star[3] = s.part.star['vx'] * s.info.kms
    star[4] = s.part.star['vy'] * s.info.kms
    star[5] = s.part.star['vz'] * s.info.kms
    star[6] = s.part.star['m'] * s.info.msun
    star[7] = s.part.star['id']
    star[8] = s.part.star['time']
    star[9] = s.part.star['metal']
    
    if 'dm' in s.part.pt:
        dm[0] = s.part.dm['x']
        dm[1] = s.part.dm['y']
        dm[2] = s.part.dm['z']
        dm[3] = s.part.dm['vx'] * s.info.kms
        dm[4] = s.part.dm['vy'] * s.info.kms
        dm[5] = s.part.dm['vz'] * s.info.kms
        dm[6] = s.part.dm['m'] * s.info.msun
        dm[7] = s.part.dm['id']
        s.part.dm = None
        
    s.part.star = None
    
    s.hydro.amr2cell(lmax=lmax)
    celldata = s.hydro.cell
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
    
    cell[0] = celldata['x']
    cell[1] = celldata['y']
    cell[2] = celldata['z']
    cell[3] = celldata['dx']
    cell[4] = celldata['var0']
    cell[5] = celldata['var1']
    cell[6] = celldata['var2']
    cell[7] = celldata['var3']
    cell[8] = celldata['var4']
    cell[9] = celldata['var5']
    # Why cells are in physical unit from the debining whereas part is in code unit? 
    
    import time
    t1 = time.time()

    m = mp.Manager()
    out_q = m.Queue()

    nh = len(h.data)

    pool = mp.Pool(processes=ncore)
    for i in range(nh):
        print(i,"-th started")
        pool.apply_assync(mk_gal, args=(h.data[i], out_q, s.info, i))

    time.sleep(2)
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
        print("dict out", i)
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
    import pandas as pd
    catalog = pd.DataFrame(dictout).to_records()    
    # pickling also works!
    # The problem was memory being locked!
    import pickle
    with open(fcat, 'wb') as f:
        pickle.dump(catalog, f)
