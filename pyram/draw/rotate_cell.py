import numpy as np
import pickle
import pyram

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
from pyram.tree import tmtree as tmt 
from pyram.utils.cosmology import Timeconvert as TC
#from pyram.mock import ramski # requires ramski_module.so
from pyram import galaxymodule as gmo 

import pts.simulation as sm


kpc_in_cm = 3.0857e21
msun_in_g = 1.989e33
G = 6.674e-8  # gravitational constant ; [cm^3 g^-1 s^-2]


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


#### Draw ####

def gm2code(arr, info):
    return (arr / info.pboxsize + 0.5)


def rotation(array1, array2, theta):
    theta *= np.pi / 180
    new_array1 = np.cos(theta) * array1 - np.sin(theta) * array2
    new_array2 = np.sin(theta) * array1 + np.cos(theta) * array2
    return new_array1, new_array2

###############################################

def get_colden(theta_xy, theta_xz, theta_yz, n_sample_factor=1.0,
                directory=None, file_name='save.npy', quick=False, 
                gridrate=0.5, shift=[0, 0, 0], draw=False, save=False, verbose=False):
    """
    Rotate gas into arbitrary direction
    """

    if gridrate < 2**(-7):
        boxsize=10**2
    elif gridrate < 2**(-6):
        boxsize=3*10**2
    else:
        boxsize = 10**4

    x = np.random.randint(1000, size=boxsize*n_sample_factor)
    y = np.random.randint(1000, size=boxsize*n_sample_factor)
    z = np.random.randint(1000, size=boxsize*n_sample_factor)
    gridsize = 1000 * gridrate  #### notice that gridsize is a half of box's side length
    x, y, z = x - 500, y - 500, z - 500

    x, y = rotation(x, y, theta_xy)
    x, z = rotation(x, z, theta_xz)
    y, z = rotation(y, z, theta_yz)

    x, y, z = x + shift[0], y + shift[1], z + shift[2]

    dsort = np.where((np.sqrt(np.square(x) + np.square(y)) < gridsize * np.sqrt(2))
                     & (abs(x) <= gridsize) & (abs(y) <= gridsize))

    if draw:
        plt.show()
    else:
        pass

    z_sort = np.where( abs(z) <= gridsize )[0]
    X_zsorted = x[z_sort]
    Y_zsorted = y[z_sort]

    min_xshift = min(X_zsorted)/2/gridsize
    max_xshift = max(X_zsorted)/2/gridsize
    min_yshift = min(Y_zsorted)/2/gridsize
    max_yshift = max(Y_zsorted)/2/gridsize

    min_xshi, min_yshi, min_zshi = -1000*np.sqrt(3)/gridsize/2/2,-1000*np.sqrt(3)/gridsize/2/2,-1000*np.sqrt(3)/gridsize/2/2
    max_xshi, max_yshi, max_zshi = 1000*np.sqrt(3)/gridsize/2/2, 1000*np.sqrt(3)/gridsize/2/2, 1000*np.sqrt(3)/gridsize/2/2
    
    base_grid_ddx = int(max(max_xshi, abs(min_xshi)))+1
    base_grid_ddy = int(max(max_yshi, abs(min_yshi)))+1
    base_grid_ddz = int(max(max_zshi, abs(min_zshi)))+1
    print("\n","######################","\n","base_grid_ddx is ",base_grid_ddx,"\n","#####################","\n")

    base_grid = np.zeros([2*base_grid_ddz+2+1, 2*base_grid_ddy+1, 2*base_grid_ddx+1])

    i = -base_grid_ddx
    while i <= base_grid_ddx:
        j = -base_grid_ddy
        while j <= base_grid_ddy:
            k = -base_grid_ddz
            while k <= base_grid_ddz:
                component_ijk = np.sum((abs(x + 2 * gridsize * i) <= gridsize) *
                                       (abs(y + 2 * gridsize * j) <= gridsize) *
                                       (abs(z + 2 * gridsize * k) <= gridsize))/boxsize

                base_grid[0][j+base_grid_ddy][i+base_grid_ddx] = i
                base_grid[1][j+base_grid_ddy][i+base_grid_ddx] = j
                base_grid[k+base_grid_ddz+2][j+base_grid_ddy][i+base_grid_ddx] = component_ijk
                #base_grid[i+base_grid_ddx][j+base_grid_ddy][k+base_grid_ddz] = component_ijk
                k = k + 1
            j = j +1
        if i%10 == 1: print("{:.2f} % \r".format(100*abs(i+base_grid_ddx)/base_grid_ddx/2))
        i = i +1

    if verbose: print(base_grid)

    if save:
        save_route = directory
        route_name = save_route+file_name

        np.save(route_name,base_grid)

    return len(dsort[0]), base_grid  # /grid/grid/1000


def smoothing(theta_xy,theta_xz,theta_yz, 
            directory=None, file_name=None, quick=None, 
            gridrate=0.5, shift=[0, 0, 0], plt=False, 
            save=False, n_sample_factor=1):
    A = []
    Base_grid = []
    i = 0
    while i < len(theta_xz):
        if i == 0:
            a, ith_base_grid = get_colden(theta_xy, theta_xz[i], theta_yz, gridrate=gridrate, shift=shift, 
                                          draw=plt, quick=quick, n_sample_factor=n_sample_factor,
                                          directory =directory,file_name=file_name, save=save )
        else:
            a, ith_base_grid = get_colden(theta_xy, theta_xz[i], theta_yz, gridrate=gridrate, shift=shift, 
                                          quick=quick, directory=directory, file_name=file_name,
                                          save=save, n_sample_factor=n_sample_factor)
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

def rotate_vector(r, v):
    J = sum(np.cross(r,v))
    print(J)
    
    theta = np.degrees(np.arccos(J[2]))
    if J[1] > 0:
        phi = np.degrees(np.arccos(J[0]/np.sqrt(J[0]**2 + J[1]**2)))
    else:
        phi = np.degrees(-np.arccos(J[0]/np.sqrt(J[0]**2 + J[1]**2)))
    print(theta, phi)
    z = J/np.sqrt(sum(J*J))
    x = np.cross(np.array([0,0,1]),z)
    x = x/np.sqrt(sum(x*x))
    y = np.cross(z,x)
    y = y/np.sqrt(sum(y*y))
    rotate = np.vstack((x,y,z)).T
    rotate = np.matrix(rotate)
    rotate = np.linalg.inv(rotate)
    R = (rotate*r.T).T
    V = (rotate*v.T).T
    R = np.array(R)
    V = np.array(V)
    return R, V

def auto_rotation_np(param, xarray, yarray, zarray,thetype=0):
    a, b, c = param[0], param[1], param[2]
    theta_xy = -1 * np.arctan2(a, b)
    yarray_mid, xarray_new = rotation_xy(theta_xy, yarray, xarray)
    y_abs = (a ** 2 + b ** 2) ** 0.5
    theta_yz = -1 * np.arctan2(y_abs, c)
    zarray_new, yarray_new = rotation_yz(theta_yz, zarray, yarray_mid)

    if thetype == 1:
        return xarray_new, yarray_new, zarray_new, theta_xy * 180 / np.pi, theta_yz * 180 / np.pi
    else:
        return xarray_new, yarray_new, zarray_new


def make_gasmap(nout, gal_cat, theta_xz_1=[0], theta_yz_1=[0],
                wdir='./', 
                skirt_out_dir=None,
                fn_png=None, 
                nvec_frac=0.2,
                n_sample_factor=1,
                save = True,
                temp_cut = 1e8, # All cells
                Nvec_gas = False,
                plts = False,
                quick = True,
                smoothing_save_dir= "./",
                plotting="plt",
                fsave_angle=None,
                do_star=False,
                radius = 2e-4,
                shrink=False):

    centers = np.array([gal_cat['x'],gal_cat['y'],gal_cat['z']])
    gid=gal_cat['id']
    idgal=gal_cat['id']
    mingrid = 2.38418579e-07
    cell_ddx_grid_size = [2**i*mingrid for i in range(8)]
    if skirt_out_dir == None: skirt_out_dir = wdir + f"gal_{nout:05d}/"

    if fn_png == None: fn_png = f"{nout}_{idgal}_faceon_gasmap.png"
    s=pyram.load.sim.Sim(nout, base=wdir)
   # #radius = 30/(100*1e3/0.704)
    s.set_ranges([[centers[0]-radius, centers[0]+radius],
                  [centers[1]-radius, centers[1]+radius],
                  [centers[2]-radius, centers[2]+radius]])
    print(f"Galaxy {gid}, centers:{centers}, radius:{radius:.6f}")
    if do_star:
        print("loading particles")
        s.add_part(ptypes=["star id pos mass vel metal age"])
        star = s.part.star
        gal = gmo.galaxy.Galaxy()
        gal.info = s.info
        gal.star = s.part.star

        xcen = np.mean(gal.star['pos'], axis=0)
        vcen = np.mean(gal.star['vel'], axis=0)
        gal.star['pos'] -= xcen
        gal.star['vel'] -= vcen

        gal.star['x'] = gm2code(gal.star['x'], gal.info)
        gal.star['y'] = gm2code(gal.star['y'], gal.info)
        gal.star['z'] = gm2code(gal.star['z'], gal.info)
        
        print("# {} stars in total".format(len(gal.star)))

        dist = np.sqrt(np.einsum("...i,...i", gal.star["pos"],gal.star["pos"]))
        ind_close = np.argsort(dist)[:int(nvec_frac*len(gal.star))]
        close_stars = gal.star[ind_close]

        params_ang = np.mean(close_stars['m'] * np.cross(close_stars["pos"], close_stars["vel"]).T, axis=1)
        close_stars = None
        print("params_ang", params_ang)


        lvec = np.sum(np.cross(gal.star['pos'][ind_close],
                               gal.star['vel'][ind_close]), axis=0)
        nvec = lvec/np.linalg.norm(lvec)

        theta = np.degrees(np.arccos(nvec[2]))
        if nvec[1] > 0:
            phi = np.degrees(np.arccos(nvec[0]/np.sqrt(nvec[0]**2 + nvec[1]**2)))
        else:
            phi = np.degrees(-np.arccos(nvec[0]/np.sqrt(nvec[0]**2 + nvec[1]**2)))
        if fsave_angle is not None: fsave_angle.write(f"{nout} {gid} {theta:.3f}  phi={phi:.3f} \n")
        print(f"theta={theta:.3f},  phi={phi:.3f}")

        if shrink:
            radius1 = radius * shrink
            s.set_ranges([[centers[0]-radius1, centers[0]+radius1],
                          [centers[1]-radius1, centers[1]+radius1],
                          [centers[2]-radius1, centers[2]+radius1]])
            print(f"Loading again... Galaxy {gid}, centers:{centers}, radius:{radius1:.6f}")
            s.add_hydro()
            cell = s.hydro.cell
            
            #gal.star = gal.star[np.where(dist < np.max(dist) * shrink)[0]]
            #print("N stellar particle reduced to", len(gal.star))

    else:
        if shrink:
            radius1 = radius * shrink
            s.set_ranges([[centers[0]-radius1, centers[0]+radius1],
                          [centers[1]-radius1, centers[1]+radius1],
                          [centers[2]-radius1, centers[2]+radius1]])
            s.add_hydro()
            cell = s.hydro.cell

        xcen = centers
        vcen = np.mean(cell["vel"], axis=0)


    cell["pos"] -= xcen
    cell['vel'] -= vcen
    cell['vel'] *= s.info.kms

    cell_rho = cell['var0'] * s.info.unit_nH
    cell_T = cell['var4'] / cell['var0'] * s.info.unit_T2
    cell_mgas = cell['var0'] * s.info.unit_d * (cell['dx'] * s.info.boxtokpc * kpc_in_cm) ** 3 / msun_in_g

    if do_star:
        cell["x"], cell["y"], cell["z"], theta_xy_0, theta_yz_0 = auto_rotation_np(params_ang,
                                                               cell["x"], cell["y"], cell["z"],thetype=1)
    else:
        print("taking angles from .ski")
        print("Reading",  skirt_out_dir+f'{gid:05d}/{nout:05d}/g{gid}_{nout}.ski')
        # get theta from .ski
        skifile = sm.SkiFile(skirt_out_dir+f'{gid:05d}/{nout:05d}/g{gid}_{nout}.ski')
        elems = skifile._tree.xpath('//MonteCarloSimulation/instrumentSystem/InstrumentSystem/instruments/FullInstrument')
        faceon = elems[0]
        theta = float(faceon.get('inclination').split('deg')[0])
        phi   = float(faceon.get('azimuth').split('deg')[0])
        theta_xy = phi -90
        theta_yz = -theta
        # Rotation taken from auto_rotation_np
        yarray_mid, cell["x"] = rotation_xy(theta_xy, cell["y"], cell["x"])
        cell["z"], cell["y"]  = rotation_yz(theta_yz, cell["z"], yarray_mid)
        theta_xy_0 = theta_xy
        theta_yz_0 = theta_yz


    print(f"theta_xy={theta_xy_0:.3f},  theta_yz={theta_yz_0:.3f}")
    if Nvec_gas:
        #######################################################################
        # Calcuate Nvec direction
        CTC = np.vectorize(cartesian_to_cylinder)

        print("\n Initiate Rotation")

        cell_x, cell_y, cell_z = cell["x"], cell["y"], cell["z"]
        #cell_vx, cell_vy, cell_vz = cell["var1"], cell["var2"], cell["var3"]
        cell_r, cell_phi, cell_z, cell_v_r, cell_v_phi, cell_v_z = CTC(cell_x, cell_y, cell_z, cell['vx'], cell['vy'], cell['vz'])

        print("Gas Cells complete \n")

        gas_mass = cell_mgas * msun_in_g

        cell_x *= gal.info.boxtokpc
        cell_y *= gal.info.boxtokpc
        cell_z *= gal.info.boxtokpc
        cell_r *= gal.info.boxtokpc

        cell_x *= kpc_in_cm
        cell_y *= kpc_in_cm
        cell_z *= kpc_in_cm
        cell_r *= kpc_in_cm

        ##################################################
        T_seperate = 6+0.25*np.log10(cell_rho) 

        cold_gas_sort = np.where(  (np.log10(cell_T) < T_seperate)   )

        cold_gas_mass = cell_mgas[cold_gas_sort]

        ang_cell = cell_r * cell_v_phi * cell_mgas
        ang_cell = np.log10(ang_cell)

        cell_x /= gal.info.boxtokpc
        cell_y /= gal.info.boxtokpc
        cell_z /= gal.info.boxtokpc
        cell_r /= gal.info.boxtokpc

        cell_x /= kpc_in_cm
        cell_y /= kpc_in_cm
        cell_z /= kpc_in_cm
        #cell_r /= kpc_in_cm ??
    ###################################################

    fixed_grid_size = cell_ddx_grid_size[0] # dx_min

    grid_ratios = [fixed_grid_size/ddx/2 for ddx in cell_ddx_grid_size]

    rot_cells = []
    for i, grid_ratio in enumerate(grid_ratios):
        abc, rotcel= smoothing(theta_xy=theta_xy_0, theta_xz=theta_xz_1, theta_yz=theta_yz_0 + theta_yz_1[0],
                             gridrate=grid_ratio, plt=plts, quick=quick,n_sample_factor=n_sample_factor,
                             directory=smoothing_save_dir, file_name=f"{i}_basegrid_quick.npy", save=save)
        rot_cells.append(rotcel)

    ######################################################################
    cell["x"], cell["z"] = rotation(cell["x"], cell["z"],theta_xz_1[0])
    cell["y"], cell["z"] = rotation(cell["y"], cell["z"],theta_yz_1[0])

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

    sub_cell = cell
    sub_cell_T = sub_cell['var4'] / sub_cell['var0'] * s.info.unit_T2

    print("len subcell", len(sub_cell))
    if len(sub_cell) < 1:
        print("something's wrong with sub_cell")

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

    ##############################################################
    npixx = np.int(np.ceil(xspan / cell_dx_min))
    npixy = np.int(np.ceil(yspan / cell_dx_min))
    npixz = np.int(np.ceil(zspan / cell_dx_min))

    dxmap = dymap = cell_dx_min

    # Mimic stellar hist2d bins
    h = [None,
         np.linspace(xrange[0], xrange[1], npixx),
         np.linspace(yrange[0], yrange[1], npixy)]

    cx = np.searchsorted(h[1], sub_cell["x"]) - 1
    cy = np.searchsorted(h[2], sub_cell["y"]) - 1

    zindex = np.linspace(zrange[0],zrange[1],npixz)

    cz = np.searchsorted(zindex, sub_cell["z"]) -1

    rot_cell_length = [len(rotcell)-2 for rotcell in rot_cells]

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
    n_subcell = len(sub_cell)
    for sc in range(n_subcell):
        rot_cell_index_sc = np.searchsorted(sub_ddx_list,sub_ddx_h[sc])
        if rot_cell_index_sc >= 5:
            pass

        rot_cell_length_sc = rot_cell_length[rot_cell_index_sc]
        rot_cell_centre_num_sc = int(rot_cell_centre_num[rot_cell_index_sc])
        rot_cell_sc = rot_cells[rot_cell_index_sc]

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

        #if sub_cell_T[sc] <= temp_cut: # Can I impose T cut later sth like: Grid[Grid["temp"] > 1e4] = 0?
        Grid[Gz_sc_low:Gz_sc_upp+1,Gy_sc_low:Gy_sc_upp+1,Gx_sc_low:Gx_sc_upp+1] += rot_cell_sc[rcsc_low:rcsc_upp+1]*sub_cell["var0"][sc]*(sub_cell["dx"][sc]**3)/(cell_dx_min**2)

        #else:
        #    pass

        #sc += 1
        if (sc*100) % n_subcell == 0:
            print("{:.1f} % \r".format(np.round(100*sc/len(sub_cell),2) ))

    file_name_grid = str(nout) + '_' +str(idgal) + '_' + str(theta_xz_1[0])+ '_' + str(theta_yz_1[0]) + '_test_grid.npy'
    np.save("./%s_small" %(file_name_grid),Grid)
    grid_prj = np.sum(Grid, axis=0)
    np.save("./%s_prj" %(file_name_grid),grid_prj)
    if plotting == "pil":
        from PIL import Image
        im = Image.fromarray(np.uint8(cm.viridis(grid_prj)*255))
        im.save(fn_png)
    elif plotting == "plt":
        fig, ax = plt.subplots()
        fig.set_size_inches(6,6)
        ax.imshow(grid_prj, norm=LogNorm())
        plt.tight_layout()
        plt.savefig(fn_png, dpi=200)

    print("Grid Saved")


