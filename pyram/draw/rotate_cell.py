import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from pyram import load, utils, tree 
import pyram.utils.samping as smp

######### constants #####################3

kpc_in_cm = 3.0857e21
msun_in_g = 1.989e33

G = 6.674e-8  # gravitational constant ; [cm^3 g^-1 s^-2]

#cell_ddx_grid_size = [2.38418579e-07,4.76837158e-07,9.53674316e-07,1.90734863e-06,3.81469727e-06,7.62939453e-06]

mingrid = 2.38418579e-07
cell_ddx_grid_size = [(2**0)*mingrid,(2**1)*mingrid,(2**2)*mingrid,(2**3)*mingrid,
                      (2**4)*mingrid,(2**5)*mingrid,(2**6)*mingrid,(2**7)*mingrid]


######################################################



def cos(theta):
    # return np.cos( theta * np.pi / 180 )
    return np.cos(theta)


def sin(theta):
    # return np.sin( theta * np.pi / 180 )
    return np.sin(theta)


#### Rotation ####
## you don't have to seperate each planes. Just for convenience ##

def rotation_xy(theta, xarray, yarray):
    x = cos(theta) * xarray - sin(theta) * yarray
    y = sin(theta) * xarray + cos(theta) * yarray
    return x, y


def rotation_xz(theta, xarray, zarray):
    x = cos(theta) * xarray - sin(theta) * zarray
    z = sin(theta) * xarray + cos(theta) * zarray
    return x, z


def rotation_yz(theta, yarray, zarray):
    y = cos(theta) * yarray - sin(theta) * zarray
    z = sin(theta) * yarray + cos(theta) * zarray
    return y, z


#### Axis Transformation ####


def cartesian_to_cylinder(xarray, yarray, zarray, vxarray, vyarray, vzarray):
    r = (xarray ** 2 + yarray ** 2) ** 0.5
    phi = np.arctan2(yarray, xarray)
    z = zarray

    v_r = np.cos(phi) * vxarray + np.sin(phi) * vyarray
    v_phi = (-np.sin(phi) * vxarray + np.cos(phi) * vyarray)
    v_z = vzarray

    return r, phi, z, v_r, v_phi, v_z


def plane_fit(x, a, b, c):
    return a * x[0] + b * x[1] + c


def auto_rotation_np(param, xarray, yarray, zarray,thetype=0):
    a, b, c = param[0], param[1], param[2]
    theta_xy = -1 * np.arctan2(a, b)
    yarray_mid, xarray_new = rotation_xy(theta_xy, yarray, xarray)
    y_abs = (a ** 2 + b ** 2) ** 0.5
    theta_yz = -1 * np.arctan2(y_abs, c)
    zarray_new, yarray_new = rotation_yz(theta_yz, zarray, yarray_mid)

    print(theta_xy * 180 / np.pi, theta_yz * 180 / np.pi)
    if thetype == 1:
        return xarray_new, yarray_new, zarray_new, theta_xy * 180 / np.pi, theta_yz * 180 / np.pi
    else:
        return xarray_new, yarray_new, zarray_new

#### Draw ####

def gm2code(arr, info):
    return (arr / info.pboxsize + 0.5)


def rotation(array1, array2, theta):
    theta *= np.pi / 180
    new_array1 = np.cos(theta) * array1 - np.sin(theta) * array2
    new_array2 = np.sin(theta) * array1 + np.cos(theta) * array2
    return new_array1, new_array2

###############################################

def get_colden(theta_xy, theta_xz, theta_yz, 
                directory=None, file_name=None, quick=False, 
                gridrate=0.5, shift=[0, 0, 0], draw=False, save=False):
    """
    Rotate gas into arbitrary direction
    """

    if gridrate < 2**(-6):
        boxsize=3*10**2
    elif gridrate < 2**(-7):
        boxsize=10**2
    else:
        boxsize = 10**4
    x = np.random.randint(1000, size=boxsize)
    y = np.random.randint(1000, size=boxsize)
    z = np.random.randint(1000, size=boxsize)
    gridsize = 1000 * gridrate  #### notice that gridsize is a half of box's side length
    x, y, z = x - 500, y - 500, z - 500

    x, y = rotation(x, y, theta_xy)
    x, z = rotation(x, z, theta_xz)
    y, z = rotation(y, z, theta_yz)

    x, y, z = x + shift[0], y + shift[1], z + shift[2]

    dsort = np.where((np.sqrt(np.square(x) + np.square(y)) < gridsize * np.sqrt(2))
                     & (abs(x) <= gridsize) & (abs(y) <= gridsize))


    if draw == True:
        plt.show()
    else:
        pass

    z_sort = np.where( (abs(z) <= gridsize) )
    X_zsorted = x[z_sort]
    Y_zsorted = y[z_sort]

    min_xshift = min(X_zsorted)/2/gridsize
    max_xshift = max(X_zsorted)/2/gridsize
    min_yshift = min(Y_zsorted)/2/gridsize
    max_yshift = max(Y_zsorted)/2/gridsize

    min_xshi = 0
    min_yshi = 0
    min_zshi = 0
    max_xshi = 0
    max_yshi = 0
    max_zshi = 0

    if quick == False:

        repeat_size = 20
        i = -repeat_size
        while i <= repeat_size:
            j = -repeat_size
            while j <= repeat_size:
                k = -repeat_size
                checknum = -1
                while k <= repeat_size:

                    if len(np.where((abs(x+2*gridsize*i) <= gridsize) & (abs(y+2*gridsize*j) <= gridsize)& (abs(z+2*gridsize*k)<=gridsize) )[0]) > 0:
                        if checknum == -1:
                            min_xshi = min(min_xshi, i)
                            min_yshi = min(min_yshi, j)
                            min_zshi = min(min_zshi, k)
                            checknum = 0
                            if min_xshi == i or min_yshi ==j or min_zshi ==k:
                                print("min_x_shift is", min_xshi)
                                print("min_y_shift is", min_yshi)
                                print("min_z_shift is", min_zshi,"\n")
                        else:
                            pass
                    else:
                        if checknum == 0:
                            max_xshi = max(max_xshi, i)
                            max_yshi = max(max_yshi, j)
                            max_zshi = max(max_zshi, k)
                            checknum = 1

                            if max_xshi == i or max_yshi == j or max_zshi == k:
                                print("max_x_shift is", max_xshi)
                                print("max_y_shift is", max_yshi)
                                print("max_z_shift is", max_zshi,"\n")
                        else:
                            pass

                    k = k + 1
                j = j +1
            print(100*abs(i+repeat_size)/repeat_size/2)
            i = i +1

            if i == repeat_size+1:
                print("############## Result #################")
                print("min_x_shift is", min_xshi)
                print("min_y_shift is", min_yshi)
                print("min_z_shift is", min_zshi, "\n")

                print("max_x_shift is", max_xshi)
                print("max_y_shift is", max_yshi)
                print("max_z_shift is", max_zshi)
                print("########################################", "\n")

        min_xshi, min_yshi, min_zshi = -1000*np.sqrt(3)/gridsize/2/2,-1000*np.sqrt(3)/gridsize/2/2,-1000*np.sqrt(3)/gridsize/2/2
        max_xshi, max_yshi, max_zshi = 1000*np.sqrt(3)/gridsize/2/2, 1000*np.sqrt(3)/gridsize/2/2, 1000*np.sqrt(3)/gridsize/2/2
    
        print("Only support False/True")
        return

    base_grid_ddx = int(max(max_xshi, abs(min_xshi)))+1
    base_grid_ddy = int(max(max_yshi, abs(min_yshi)))+1
    base_grid_ddz = int(max(max_zshi, abs(min_zshi)))+1
    print("\n","#####################################","\n","base_grid_ddx is ",base_grid_ddx,"\n","#####################################","\n")

    base_grid = np.zeros([2*base_grid_ddz+2+1 ,2*base_grid_ddy+1,2*base_grid_ddx+1])

    i = -base_grid_ddx
    while i <= base_grid_ddx:
        j = -base_grid_ddy
        while j <= base_grid_ddy:
            k = -base_grid_ddz
            while k <= base_grid_ddz:

                component_ijk = len(np.where((abs(x + 2 * gridsize * i) <= gridsize) &
                                    (abs(y + 2 * gridsize * j) <= gridsize) &
                                    (abs(z + 2 * gridsize * k) <= gridsize))[0])/boxsize

                #print(component_ijk)
                base_grid[0][j+base_grid_ddy][i+base_grid_ddx] = i
                base_grid[1][j+base_grid_ddy][i+base_grid_ddx] = j
                base_grid[k+base_grid_ddz+2][j+base_grid_ddy][i+base_grid_ddx] = component_ijk
                #base_grid[i+base_grid_ddx][j+base_grid_ddy][k+base_grid_ddz] = component_ijk
                k = k + 1
            j = j +1
        print(100*abs(i+base_grid_ddx)/base_grid_ddx/2)
        i = i +1

    print(base_grid)

    if save == True:
        save_route = directory#"/home/jangjk/PycharmProjects/PP/week3/base_grid/"
        #file_name = "first_try_basegrid.npy"
        route_name = save_route+file_name

        np.save(route_name,base_grid)

    return len(dsort[0]), base_grid  # /grid/grid/1000


def smoothing(theta_xy,theta_xz,theta_yz, directory, file_name, quick, gridrate=0.5, shift=[0, 0, 0], plt=False, save=False):
    A = []
    Base_grid = []
    i = 0
    while i < len(theta_xz):
        if i == 0:
            a, ith_base_grid = get_colden(theta_xy, theta_xz[i], theta_yz, gridrate=gridrate, shift=shift, draw=plt, quick=quick,
                                          directory =directory,file_name=file_name, save=save )
        else:
            a, ith_base_grid = get_colden(theta_xy, theta_xz[i], theta_yz, gridrate=gridrate, shift=shift, quick=quick,
                                          directory=directory, file_name=file_name ,save=save)
        A.append(a)

        if len(theta_xz) == 1:
            Base_grid = ith_base_grid
        elif len(theta_xz) > 1:
            Base_grid.append(ith_base_grid)
        else:
            pass

        i += 1
        print(np.round(100 * i / len(theta_xz), 2))

    return np.array(A), Base_grid


def load_gal(nout, idgal,theta_xz_1, theta_yz_1, save_directory,fixed_idgal,
             boxlim, celllim, boxlim_xbot, boxlim_xupp, boxlim_ybot, boxlim_yupp,
             drawlim_xbot, drawlim_xupp, drawlim_ybot, drawlim_yupp,drawpix,
             saving=True, fuv=False, BVR=False, SDSS=False, JHK=False):

    directory = '/storage1/NewHorizon/'
    directory5 = '/storage5/NewHorizon/'
    directory_galactica = '/storage5/Galactica/galaxy_14667/TREE_STARS/'
    try:
        gal = load.rd_GM.rd_gal(nout, idgal, wdir=directory)
        s = load.sim.Sim(nout, base=directory)
    except:
        try:    
            gal = load.rd_GM.rd_gal(nout, idgal, wdir=directory5)
            s = load.sim.Sim(nout, base=directory5)
        except:
            gal = load.rd_GM.rd_gal(nout, idgal, wdir=directory_galactica,fname="galactica")
            s = load.sim.Sim(nout, base=directory_galactica,data_dir='galactica')
            

    gal.header['xg'] = gm2code(gal.header['xg'], gal.info)
    gas.star = None
    """
    gal.star['x'] = gm2code(gal.star['x'], gal.info)
    gal.star['y'] = gm2code(gal.star['y'], gal.info)
    gal.star['z'] = gm2code(gal.star['z'], gal.info)

    xc, yc, zc = gal.header['xg']

    star = gal.star


    xyzsort = np.where( 
                        (abs(star["x"]-xc) <  boxlim/gal.info.boxtokpc) & 
                        (abs(star["y"]-yc) <  boxlim/gal.info.boxtokpc) & 
                        (abs(star["z"]-zc) <  boxlim/gal.info.boxtokpc)  
                      )

    print('')
    print("Loading star particles...")


    star = star[xyzsort]

    star["x"] -= xc
    star["y"] -= yc
    star["z"] -= zc


    star['vx'] -= gal.header['vg'][0]
    star['vy'] -= gal.header['vg'][1]
    star['vz'] -= gal.header['vg'][2]
    """

    ###############################################################################################3

    # radius = 0.5 * max([gal.star['x'].ptp(), gal.star['y'].ptp(), gal.star['z'].ptp()])
    radius = celllim / gal.info.boxtokpc
    Rlim_sim = radius  # region setting


    s.set_ranges([[xc - Rlim_sim, xc + Rlim_sim], [yc - Rlim_sim, yc + Rlim_sim], [zc - Rlim_sim, zc + Rlim_sim]])


    ############################################################################################
    # Load stellar particle
    """
    s.add_part(ptypes=['dm id pos vel mass'],fortran=True)

    dm = s.part.dm
    dmsort = np.where( (dm['m']>0) & (dm['id'] >=0) )
    dm = dm[dmsort]

    
    dm["x"] -= xc
    dm["y"] -= yc
    dm["z"] -= zc
    dm["m"] *= s.info.msun

    #dm["x"] *= s.info.boxtokpc
    #dm["y"] *= s.info.boxtokpc
    #dm["z"] *= s.info.boxtokpc

    """
    ################################################################################################
    # Load cell data
    s.add_hydro()

    #################### load gas cells #############################
    print('')
    print('Loading gas cells...')

    ind_g = np.where(np.square(s.hydro.cell['x'] - xc) +
                     np.square(s.hydro.cell['y'] - yc) +
                     np.square(s.hydro.cell['z'] - zc) < np.square(Rlim_sim))[0]

    cell = s.hydro.cell[ind_g]
    #print(cell["x"])
    #print(len(cell["x"]))
    
    cell['x'] -= xc
    cell['y'] -= yc
    cell['z'] -= zc


    cell['var1'] = cell['var1'] * s.info.kms - gal.header['vg'][0]
    cell['var2'] = cell['var2'] * s.info.kms - gal.header['vg'][1]
    cell['var3'] = cell['var3'] * s.info.kms - gal.header['vg'][2]



    cell_rho = cell['var0'] * s.info.unit_nH
    cell_T = cell['var4'] / cell['var0'] * s.info.unit_T2
    cell_mgas = cell['var0'] * s.info.unit_d * (cell['dx'] * s.info.boxtokpc * kpc_in_cm) ** 3 / msun_in_g

    print('complete')
    

    #######################################################################
    # Calcuate Nvec direction
    """
    print('')
    print("Loading star particles...")


    ang_x = star["m"] * (star["y"] * star["vz"] - star["z"] * star["vy"])
    ang_y = star["m"] * (star["z"] * star["vx"] - star["x"] * star["vz"])
    ang_z = star["m"] * (star["x"] * star["vy"] - star["y"] * star["vx"])

    params_ang = [np.mean(ang_x),
                  np.mean(ang_y),
                  np.mean(ang_z)]


    """
    #star["x"], star["y"], star["z"], theta_xy_0,theta_yz_0 = auto_rotation_np(params_ang,
    #                                                       star["x"], star["y"], star["z"],thetype=1)

    cell["x"], cell["y"], cell["z"], theta_xy_0, theta_yz_0 = auto_rotation_np(params_ang,
                                                           cell["x"], cell["y"], cell["z"],thetype=1)

    cell["var1"], cell["var2"], cell["var3"], theta_xy_0, theta_yz_0 = auto_rotation_np(params_ang, cell["var1"], cell["var2"],
                                                                    cell["var3"],thetype=1)

    CTC = np.vectorize(cartesian_to_cylinder)

    print('')
    print("Initiate Rotation")

    #x, y, z = star['x'],star['y'],star['z']
    #vx, vy, vz = auto_rotation_np(params_ang, star["vx"], star["vy"], star["vz"])
    #r, phi, z, v_r, v_phi, v_z = CTC(x, y, z, vx, vy, vz)

    #print("Star particles complete")
    #print("")

    #dm_x, dm_y, dm_z = auto_rotation_np(params_ang, dm["x"], dm["y"], dm["z"])
    #dm_vx, dm_vy, dm_vz = auto_rotation_np(params_ang, dm["vx"], dm["vy"], dm["vz"])
    #dm_r, dm_phi, dm_z, dm_v_r, dm_v_phi, dm_v_z = CTC(dm_x, dm_y, dm_z, dm_vx, dm_vy, dm_vz)

    #print("DM particles complete")
    #print("")

    cell_x, cell_y, cell_z = cell["x"], cell["y"], cell["z"]
    cell_vx, cell_vy, cell_vz = cell["var1"], cell["var2"], cell["var3"]
    cell_r, cell_phi, cell_z, cell_v_r, cell_v_phi, cell_v_z = CTC(cell_x, cell_y, cell_z, cell_vx, cell_vy, cell_vz)

    print("Gas Cells complete")
    print("")


    ##################################################

    #star_mass = star['m']*1e11*msun_in_g
    #dm_mass = dm['m']*msun_in_g
    gas_mass = cell_mgas * msun_in_g

    #x *= gal.info.boxtokpc
    #y *= gal.info.boxtokpc
    #z *= gal.info.boxtokpc
    #r *= gal.info.boxtokpc
    
    cell_x *= gal.info.boxtokpc
    cell_y *= gal.info.boxtokpc
    cell_z *= gal.info.boxtokpc
    cell_r *= gal.info.boxtokpc

    #dm_x *= gal.info.boxtokpc
    #dm_y *= gal.info.boxtokpc
    #dm_z *= gal.info.boxtokpc
    #dm_r *= gal.info.boxtokpc

    #x *= kpc_in_cm
    #y *= kpc_in_cm
    #z *= kpc_in_cm
    #r *= kpc_in_cm

    #dm_x *= kpc_in_cm
    #dm_y *= kpc_in_cm
    #dm_z *= kpc_in_cm
    #dm_r *= kpc_in_cm

    cell_x *= kpc_in_cm
    cell_y *= kpc_in_cm
    cell_z *= kpc_in_cm
    cell_r *= kpc_in_cm
    
    #vx *= 1e5
    #vy *= 1e5
    #vz *= 1e5

    #v_r *= 1e5
    #v_phi *= 1e5
    #v_z *= 1e5

    ##################################################
    T_seperate = 6+0.25*np.log10(cell_rho)

    cold_gas_sort = np.where(  (np.log10(cell_T) < T_seperate)   )

    cold_gas_mass = cell_mgas[cold_gas_sort]


    ang_cell = cell_r * cell_v_phi * cell_mgas
    ang_cell = np.log10(ang_cell)


    plt.hist(ang_cell[cold_gas_sort],bins=500,histtype='step',color='b')
    plt.hist(ang_cell,bins=500,histtype='step',color='orange')
    #plt.show()


    """
    sig_r = ( sum( ( v_r - sum(v_r)/len(v_r) )**2 ) / len(v_r) )**0.5
    sig_phi = ( sum( ( v_phi - sum(v_phi)/len(v_phi) )**2 ) / len(v_phi) )**0.5
    sig_z = ( sum( ( v_z - sum(v_z)/len(v_z) )**2 ) / len(v_z) )**0.5

    sig_r = np.std(v_r)
    sig_phi = np.std(v_phi)
    sig_z = np.std(v_z)
    sig_cylinder =  ( (sig_r**2 + sig_phi**2 + sig_z**2) /3 )**0.5

    VoS = abs(np.mean(v_phi)) / sig_cylinder

    print("")
    print( "V/sigma = %s" %(VoS) )
    print('')

    ##################################################
    tc = utils.cosmology.Timeconvert(info=gal.info)
    starage = tc.time2gyr(star['time'],z_now = gal.info.zred)

    r_new = list(r) + list(dm_r) + list(cell_r)
    z_new = list(z) + list(dm_z) + list(cell_z)
    mass_new = list(star_mass) + list(dm_mass) + list(gas_mass)

    r_new = np.array(r_new)
    z_new = np.array(z_new)
    mass_new = np.array(mass_new)

    radsort = np.argsort(np.square(r_new)+np.square(z_new))
    r_new = r_new[radsort]
    z_new = z_new[radsort]
    mass_new2 = mass_new[radsort]

    cum_mass = np.cumsum(mass_new2)

    dist = np.sqrt(np.square(r)+np.square(z))
    dist_bin=np.sqrt(np.square(r_new)+np.square(z_new))

    rad_index = np.digitize(dist,bins=dist_bin)-1
    rsort = np.where( (r>0) )

    KE = 0.5*star_mass[rsort]*(np.square(v_phi[rsort])+np.square(v_z[rsort])+np.square(v_r[rsort]))
    PE = -1*G*cum_mass[rad_index][rsort]*star_mass[rsort]/np.sqrt(np.square(r[rsort])+np.square(z[rsort]))
    Etot = KE + PE

    v_cir = np.sqrt(-2*Etot/star_mass[rsort])

    J_z = r[rsort]*v_phi[rsort]
    J_cir = G*cum_mass[rad_index][rsort]/v_cir

    cir_param = J_z/J_cir

    print(cir_param)

    disc_sort = np.where( (cir_param > 0.5) )
    sph_sort = np.where( (cir_param < 0.5) )


    DtT = len(disc_sort[0])/(len(sph_sort[0])+len(disc_sort[0]))
    print(DtT)

    plt.scatter(r[rsort]/kpc_in_cm,cir_param,alpha=0.5,c=starage[rsort],cmap='jet')
    #plt.savefig("./image/%s/%s" %(fixed_idgal,str(nout)+".png"))
    """
    ##################################################

    #x /= gal.info.boxtokpc
    #y /= gal.info.boxtokpc
    #z /= gal.info.boxtokpc
    #r /= gal.info.boxtokpc

    cell_x /= gal.info.boxtokpc
    cell_y /= gal.info.boxtokpc
    cell_z /= gal.info.boxtokpc
    cell_r /= gal.info.boxtokpc

    #dm_x /= gal.info.boxtokpc
    #dm_y /= gal.info.boxtokpc
    #dm_z /= gal.info.boxtokpc
    #dm_r /= gal.info.boxtokpc

    #x /= kpc_in_cm
    #y /= kpc_in_cm
    #z /= kpc_in_cm

    #dm_x /= kpc_in_cm
    #dm_y /= kpc_in_cm
    #dm_z /= kpc_in_cm

    cell_x /= kpc_in_cm
    cell_y /= kpc_in_cm
    cell_z /= kpc_in_cm

    #vx /= 1e5
    #vy /= 1e5
    #vz /= 1e5


    ##################################################
    """
    xyzsort = np.where( ( star["x"] >= boxlim_xbot/gal.info.boxtokpc ) &
                        ( star["x"] <  boxlim_xupp/gal.info.boxtokpc ) &
                        ( star["y"] >= boxlim_ybot/gal.info.boxtokpc ) &
                        ( star["y"] <  boxlim_yupp/gal.info.boxtokpc ) )

            
    star = star[xyzsort]
    """



    ###################################################

    fixed_grid_size = cell_ddx_grid_size[0] # dx_min

    changed_grid_size0 = cell_ddx_grid_size[0]
    changed_grid_size1 = cell_ddx_grid_size[1]
    changed_grid_size2 = cell_ddx_grid_size[2]
    changed_grid_size3 = cell_ddx_grid_size[3]
    changed_grid_size4 = cell_ddx_grid_size[4]
    changed_grid_size5 = cell_ddx_grid_size[5]
    changed_grid_size6 = cell_ddx_grid_size[6]
    changed_grid_size7 = cell_ddx_grid_size[7]

    grid_rate0 = fixed_grid_size / changed_grid_size0 / 2
    grid_rate1 = fixed_grid_size / changed_grid_size1 / 2
    grid_rate2 = fixed_grid_size / changed_grid_size2 / 2
    grid_rate3 = fixed_grid_size / changed_grid_size3 / 2
    grid_rate4 = fixed_grid_size / changed_grid_size4 / 2
    grid_rate5 = fixed_grid_size / changed_grid_size5 / 2
    grid_rate6 = fixed_grid_size / changed_grid_size6 / 2
    grid_rate7 = fixed_grid_size / changed_grid_size7 / 2

    grid_size0 = int(1000 * grid_rate0)
    grid_size1 = int(1000 * grid_rate1)
    grid_size2 = int(1000 * grid_rate2)
    grid_size3 = int(1000 * grid_rate3)
    grid_size4 = int(1000 * grid_rate4)
    grid_size5 = int(1000 * grid_rate5)
    grid_size6 = int(1000 * grid_rate6)
    grid_size7 = int(1000 * grid_rate7)

    shift0 = [0, 0, 0]
    shift1 = [0, 0, 0]
    shift2 = [0, 0, 0]
    shift3 = [0, 0, 0]
    shift4 = [0, 0, 0]
    shift5 = [0, 0, 0]
    shift6 = [0, 0, 0]
    shift7 = [0, 0, 0]

    #####################################################

    size = 0
    plts = True
    plts = False
    quick = True
    save = False
    directories = "/home/jangjk816/Project/Mock/week7/base_grid/"

    ######################################################

    O, rot_cell0 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                             gridrate=grid_rate0, shift=shift0, plt=plts, quick=quick,
                             directory=directories, file_name="0_basegrid_quick.npy", save=save)

    A, rot_cell1 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                             gridrate=grid_rate1, shift=shift1, plt=plts, quick=quick,
                             directory=directories, file_name="1_basegrid_quick.npy", save=save)

    B, rot_cell2 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                             gridrate=grid_rate2, shift=shift2, plt=plts, quick=quick,
                             directory=directories, file_name="2_basegrid_quick.npy", save=save)

    C, rot_cell3 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                             gridrate=grid_rate3, shift=shift3, plt=plts, quick=quick,
                             directory=directories, file_name="3_basegrid_quick.npy", save=save)

    D, rot_cell4 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                             gridrate=grid_rate4, shift=shift4, plt=plts, quick=quick,
                             directory=directories, file_name="4_basegrid_quick.npy", save=save)

    E, rot_cell5 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                             gridrate=grid_rate5, shift=shift5, plt=plts, quick=quick,
                             directory=directories, file_name="5_basegrid_quick.npy", save=save)


    F, rot_cell6 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                                         gridrate=grid_rate6, shift=shift6, plt=plts, quick=quick,
                                                                      directory=directories, file_name="6_basegrid_quick.npy", save=save)

    Gs, rot_cell7 = smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                                         gridrate=grid_rate7, shift=shift7, plt=plts, quick=quick,
                                                                      directory=directories, file_name="7_basegrid_quick.npy", save=save)



    ######################################################################
    theta_xz = theta_xz_1[0]
    #star["x"], star["z"] = rotation(star["x"], star["z"],theta_xz)
    cell["x"], cell["z"] = rotation(cell["x"], cell["z"],theta_xz)

    
    theta_yz = theta_yz_1[0]
    #star["y"], star["z"] = rotation(star["y"], star["z"],theta_yz)
    cell["y"], cell["z"] = rotation(cell["y"], cell["z"],theta_yz)



    print("\n","rotation complete","\n")

    ########################################################################


    cell_dx_min = np.min(cell["dx"])
    cell_dx_max = np.max(cell["dx"])
    
    ddx_h = cell["dx"] * 0.5
    print("min_ddx_h =",min(ddx_h),"\n","max_ddx_h =",max(ddx_h))

    dx1sort = np.where((ddx_h >= 1e-8) & (ddx_h < 1e-7))
    dx2sort = np.where((ddx_h >= 1e-7) & (ddx_h < 1e-6))
    dx3sort = np.where((ddx_h >= 1e-6) & (ddx_h < 1e-5))

    print("Number of cells in ddx:[1e-8,1e-7) = %s" % (len(ddx_h[dx1sort])))
    print("Number of cells in ddx:[1e-7,1e-6) = %s" % (len(ddx_h[dx2sort])))
    print("Number of cells in ddx:[1e-6,1e-5) = %s" % (len(ddx_h[dx3sort])))
    print("total number of cell = %s" % (len(cell)))
    print("total number of selected cell = %s" % (len(ddx_h[dx1sort]) + len(ddx_h[dx2sort]) + len(ddx_h[dx3sort])))    


    sub_cell = cell"""[((cell["x"] + ddx_h) > star["x"].min()) * (cell["x"] - ddx_h < star["x"].max()) * \
                    ((cell["y"] + ddx_h) > star["y"].min()) * (cell["y"] - ddx_h < star["y"].max()) * \
                    ((cell["z"] + ddx_h) < star["z"].max())]"""

    sub_cell_T = sub_cell['var4'] / sub_cell['var0'] * s.info.unit_T2


    print("len subcell", len(sub_cell))
    if len(sub_cell) < 1:
        print("something's wrong with sub_cell")
        return 


    ########################################################################


    sub_ddx_h = sub_cell["dx"] * 0.5
    sub_ddx_list = np.unique(sub_ddx_h)
    print("\n","sub_cell's dx/2 size =",sub_ddx_list,"\n")


    xl = sub_cell["x"] - sub_ddx_h*np.sqrt(3)
    xr = sub_cell["x"] + sub_ddx_h*np.sqrt(3)
    yl = sub_cell["y"] - sub_ddx_h*np.sqrt(3)
    yr = sub_cell["y"] + sub_ddx_h*np.sqrt(3)
    zl = sub_cell["z"] - sub_ddx_h*np.sqrt(3)
    zr = sub_cell["z"] + sub_ddx_h*np.sqrt(3)

    print("min_xl = %s, max_xl = %s" % (min(xl), max(xl)))
    print("min_xr = %s, max_xr = %s" % (min(xr), max(xr)))
    print("min_yl = %s, max_yl = %s" % (min(yl), max(yl)))
    print("min_yr = %s, max_yr = %s" % (min(yr), max(yr)))
    print("min_zl = %s, max_zl = %s" % (min(zl), max(zl)))
    print("min_zr = %s, max_zr = %s" % (min(zr), max(zr)))

    xrange = (xl.min(), xr.max())
    yrange = (yl.min(), yr.max())
    zrange = (zl.min(), zr.max())
    xspan = xrange[1] - xrange[0]
    yspan = yrange[1] - yrange[0]
    zspan = zrange[1] - zrange[0]

    ############################################################
    vmin_x, vmax_x = xrange[0], xrange[1]
    vmin_y, vmax_y = yrange[0], yrange[1]
    vmin_z, vmax_z = zrange[0], zrange[1]
    ##############################################################

    """
    ssort = np.where((star["x"] >= vmin_x) & (star["x"] <= vmax_x) &
                     (star["y"] >= vmin_y) & (star["y"] <= vmax_y))

    
    print("star_x_min and max boxsize =", min(star["x"][ssort]), max(star["x"][ssort]))
    print("star_y_min and max boxsize =", min(star["y"][ssort]), max(star["y"][ssort]))

    x_real_span = max(star["x"][ssort]) - min(star["x"][ssort])
    y_real_span = max(star["y"][ssort]) - min(star["y"][ssort])
    z_real_span = max(star["z"][ssort]) - min(star["z"][ssort])
    
    x_real_span = x_real_span*gal.info.boxtokpc
    y_real_span = y_real_span*gal.info.boxtokpc
    z_real_span = z_real_span*gal.info.boxtokpc

    print("\n","x_real_span is ",x_real_span)
    print("\n","y_real_span is ",y_real_span)
    print("\n","z_real_span is ",z_real_span)

    print(xrange)
    print(yrange)

    print("\n","xspan is ",xspan*gal.info.boxtokpc)
    print("\n","yspan is ",yspan*gal.info.boxtokpc) 
    print("\n","zspan is ",zspan*gal.info.boxtokpc)

    npixx = np.int(np.ceil(xspan / cell_dx_min))
    npixy = np.int(np.ceil(yspan / cell_dx_min))
    npixz = np.int(np.ceil(zspan / cell_dx_min))
    print("\n","(npixx,npixy,npixz) =","(",npixx,",",npixy,",",npixz,")","\n")


    h = np.histogram2d(star["x"], star["y"],
                       bins=[npixx, npixy],
                       range=[xrange, yrange])

    dxmap = dymap = cell_dx_min

    # Sort stars
    sx = np.searchsorted(h[1], star["x"]) - 1
    sy = np.searchsorted(h[2], star["y"]) - 1

    cx = np.searchsorted(h[1], sub_cell["x"]) - 1
    cy = np.searchsorted(h[2], sub_cell["y"]) - 1

    zindex = np.linspace(zrange[0],zrange[1],npixz)
    
    sz = np.searchsorted(zindex, star["z"]) -1
    cz = np.searchsorted(zindex, sub_cell["z"]) -1

    """
    rot_cell_length = [len(rot_cell0)-2,len(rot_cell1)-2,
                        len(rot_cell2)-2,len(rot_cell3)-2,
                          len(rot_cell4)-2,len(rot_cell5)-2,
                            len(rot_cell6)-2,len(rot_cell7)-2]
    
    
    
    print(rot_cell_length)
    
    x_padding_layer = int(max(rot_cell_length)/2)
    y_padding_layer = int(max(rot_cell_length)/2)
    z_padding_layer = int(max(rot_cell_length)/2) 
    
    print("... Preparing for stacking column density of cells ...")
    
    Grid = np.zeros([npixz+z_padding_layer*2+2,npixy+y_padding_layer*2,npixx+x_padding_layer*2]) # New gas distribution 
    
    ######### Notice that the order of index is not [x,y,z] but [z,y,x] ###################
    ######### Additional 2 arrays in z component is for x,y indexing when we draw maps ############ 
    
    print("x_length of empty grid = ",len(Grid[0][0]))
    print("y_length of empty grid = ",len(Grid[0]))
    print("z_length of empty grid = ",len(Grid))

    rot_cell_length = np.array(rot_cell_length)
    rot_cell_centre_num = (rot_cell_length-1)/2-1+1+2 ### -1 : indexing is start from zero, not one.
                                                        ### +2 : first and second layer is xy indexing layer.
                                                          ### +1 : length + center(1) + length
    sc = 0 #+len(sub_cell)
    while sc < len(sub_cell):
        rot_cell_index_sc = np.searchsorted(sub_ddx_list,sub_ddx_h[sc])
        if rot_cell_index_sc >= 5:
            pass

        rot_cell_length_sc = rot_cell_length[rot_cell_index_sc]
        rot_cell_centre_num_sc = int(rot_cell_centre_num[rot_cell_index_sc])
        #print(sub_ddx_h[sc],sub_ddx_list[np.searchsorted(sub_ddx_list,sub_ddx_h[sc])])
        rot_cell_sc = []
        if rot_cell_index_sc == 0:
            rot_cell_sc = rot_cell0
        elif rot_cell_index_sc == 1:
            rot_cell_sc = rot_cell1
        elif rot_cell_index_sc == 2:
            rot_cell_sc = rot_cell2
        elif rot_cell_index_sc == 3:
            rot_cell_sc = rot_cell3
        elif rot_cell_index_sc == 4:
            rot_cell_sc = rot_cell4
        elif rot_cell_index_sc == 5:
            rot_cell_sc = rot_cell5
        elif rot_cell_index_sc == 6:
            rot_cell_sc = rot_cell6
        elif rot_cell_index_sc == 7:
            rot_cell_sc = rot_cell7
        else:
            print("Error occured at calculating part of indexing sc'th sub_cell's index")
            break
        
        Gx_sc = cx[sc]+x_padding_layer
        Gy_sc = cy[sc]+y_padding_layer
        Gz_sc = cz[sc]+z_padding_layer+2 ### +2 for removing xy-indexing layer
        
        test_rcsc = np.zeros([2+rot_cell_length_sc,rot_cell_length_sc,rot_cell_length_sc])
        

        rcsc_low = rot_cell_centre_num_sc-int((rot_cell_length_sc-1)/2)
        rcsc_upp = rot_cell_centre_num_sc+int((rot_cell_length_sc-1)/2)

        test_rcsc[rcsc_low:rcsc_upp+1] = rot_cell_sc[rcsc_low:rcsc_upp+1]

        Gx_sc_low = Gx_sc - int((rot_cell_length_sc-1)/2)
        Gx_sc_upp = Gx_sc + int((rot_cell_length_sc-1)/2)

        Gy_sc_low = Gy_sc - int((rot_cell_length_sc-1)/2)
        Gy_sc_upp = Gy_sc + int((rot_cell_length_sc-1)/2)
        
        Gz_sc_low = Gz_sc - int((rot_cell_length_sc-1)/2)
        Gz_sc_upp = Gz_sc + int((rot_cell_length_sc-1)/2)
       

        
        if sub_cell_T[sc] <= 8000:
            Grid[Gz_sc_low:Gz_sc_upp+1,Gy_sc_low:Gy_sc_upp+1,Gx_sc_low:Gx_sc_upp+1] += rot_cell_sc[rcsc_low:rcsc_upp+1]*sub_cell["var0"][sc]*(sub_cell["dx"][sc]**3)/(cell_dx_min**2)

        else:
            pass

        sc += 1
        print( np.round(100*sc/len(sub_cell),2) )
    
    j = 0  
    while j < len(Grid[0]):
        i = 0
        while i < len(Grid[0][0]):
            Grid[0][j][i] = i
            Grid[1][j][i] = j
            i += 1
        j += 1
    
    
    sumup = 0
    k = 2
    while k < len(Grid):
        sumup_k = sum(sum(Grid[k]))
        sumup = sumup + sumup_k
        k = k +1

    file_name_grid = str(nout) + '_' +str(idgal) + '_' + str(theta_xz_1[0])+ '_' + str(theta_yz_1[0]) + '_test_grid.npy'
    np.save("./base_grid/%s" %(file_name_grid),Grid)
    print("Grid Saved")
    
    #Grid = np.load("./base_grid/%s" %(file_name_grid))
    """
    ############################## Calculate colden for every star ###########################
    nstar = len(star)
    colden_star = np.zeros(nstar)
    
    i = 0 #+ len(star)
    while i < len(star):
        sx_i = sx[i]
        sy_i = sy[i]
        sz_i = sz[i]
        
        #### 이 부분이 중요함. 방향이 어디이냐에 따라서 나오는 colden 결과가 완전 반대로 나올 수도 있다. ###
        colden_star[i] = sum(Grid[z_padding_layer+2:z_padding_layer+2+sz_i,y_padding_layer+sy_i,x_padding_layer+sx_i])
        #colden_star[i] = sum(Grid[z_padding_layer+2+sz_i+1:len(Grid),y_padding_layer+sy_i,x_padding_layer+sx_i])
        
        print( np.round(100*i/len(star),2) )
        i += 1
    
    colden = colden_star * gal.info.unit_nH * gal.info.unit_l #/ gal.info.boxtokpc
    
    print('')
    print('snapshot Number = %s' %(nout) )
    print('galaxy Number = %s' %(idgal) )
    print('Rotate angle along y-axis = %s' %(theta_xz_1[0]) )
    print('Rotate angle along x-axis = %s' %(theta_yz_1[0]) )
    print("Sum up all Grid componets is ",sumup)
    print('')

    file_name_colden = str(nout) + '_' +str(idgal) + '_' + str(theta_xz_1[0])+ '_' + str(theta_yz_1[0]) + '_temporary_colden.npy'
    np.save("./colden/%s" %(file_name_colden),colden)
    print("Colden Saved")
    
    #colden = np.load("./colden/%s" %(file_name_colden))
    """