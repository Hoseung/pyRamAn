import numpy as np
import math
import scipy.integrate
#import measure_run_time as mrt
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from multiprocessing import Pool

#import load
#import utils
#import tree
#import utils.sampling as smp
#from galaxymodule import quick_mock as qmc

import sys
import os

sys.path.append("/home/jangjk816/pymockevol")
sys.path.append("/home/jangjk816/Project/Mock/week3")
sys.path
#from galaxymodule import quick_mock_original as qmc_ori
#import try_quick_mock as tqmc

#from mload import rotation 

import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from astropy.io import fits






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
    #mingrid = 2.38418579e-07
    #cell_ddx_grid_size = [2**i*mingrid for i in range(8)]
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






















def rotation(array1, array2, theta):
    theta *= np.pi / 180
    new_array1 = np.cos(theta) * array1 - np.sin(theta) * array2
    new_array2 = np.sin(theta) * array1 + np.cos(theta) * array2
    return new_array1, new_array2

############################################################################

def init_shared(_shared, size):
    global shared_mat
 
    shared_mat = np.frombuffer(_shared.get_obj()).reshape((size, size))
 

def base_grid(theta_xy, theta_xz, theta_yz, gridrate=0.5, shift=[0, 0, 0], quick=True, draw=False, save=False, save_directory = '', save_filename = ''):

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

    min_xshi, min_yshi, min_zshi = -1000*np.sqrt(3)/gridsize/2/2,-1000*np.sqrt(3)/gridsize/2/2,-1000*np.sqrt(3)/gridsize/2/2
    max_xshi, max_yshi, max_zshi = 1000*np.sqrt(3)/gridsize/2/2, 1000*np.sqrt(3)/gridsize/2/2, 1000*np.sqrt(3)/gridsize/2/2

    base_grid_ddx = int(max(max_xshi, abs(min_xshi)))+1
    base_grid_ddy = int(max(max_yshi, abs(min_yshi)))+1
    base_grid_ddz = int(max(max_zshi, abs(min_zshi)))+1

    #print("\n","#####################################","\n","base_grid_ddx is ",base_grid_ddx,"\n","#####################################","\n")

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


                base_grid[0][j+base_grid_ddy][i+base_grid_ddx] = i
                base_grid[1][j+base_grid_ddy][i+base_grid_ddx] = j
                base_grid[k+base_grid_ddz+2][j+base_grid_ddy][i+base_grid_ddx] = component_ijk
                k = k + 1
            j = j +1
        #print(100*abs(i+base_grid_ddx)/base_grid_ddx/2)
        i = i +1


    if save == True:
        save_route = save_directory
        route_name = save_route+save_filename

        np.save(route_name,base_grid)

    return base_grid



def get_colden(gal,star, cell, cell_rho, cell_T, theta_xy0, theta_xz1, theta_yz0, theta_yz1, grid_number = 6,cfactor=2, processing_number = 8,Tcut=20000,rhocut=1, mingrid=2.38418579e-07):
    """
    index for the internal structure of input parameter 'star' and 'cell' is customized for NH data.
    If you want to use this code on another simulation output data, you should change the sub-index of input parameters.
    """

    rot_cell = []
    rot_cell_length = []

    #i = 0
    #while i < grid_number:
    for i in range(i, grid_number):
        grid_rate_i = 1/ 2**i /2
        rot_grid_i = base_grid(theta_xy = theta_xy0, theta_xz = theta_xz1, theta_yz = theta_yz0 + theta_yz1, gridrate = grid_rate_i)
        rot_cell.append(rot_grid_i)
        rot_cell_length.append(len(rot_grid_i)-2)
        #i += 1

    star["x"], star["z"] = rotation(star["x"], star["z"], theta_xz1)
    cell["x"], cell["z"] = rotation(cell["x"], cell["z"], theta_xz1)

    #theta_yz = theta_yz1
    star["y"], star["z"] = rotation(star["y"], star["z"], theta_yz1)
    cell["y"], cell["z"] = rotation(cell["y"], cell["z"], theta_yz1)

    print("\n", "rotation complete", "\n")

    
    cell_dx_min = cfactor*np.min(cell["dx"])    
    ddxsort = (cell_T < Tcut) * (cell_rho > rhocut)

    ddx_h = cell["dx"] /2
    ddx_h = max(cell["dx"][ddxsort])/2
    
    cell_ddx_list = []
    #d = 1
    #while d <= grid_number:
    for d in range(1, grid_number+1):
        cell_ddx_list.append(cell_dx_min*d/2)
        #d += 1

    ###########################################################
    # relevant cells. I may skip this part or do it beforehand.
    ind_subcell = np.where(((cell["x"] + ddx_h) > star["x"].min()) * (cell["x"] - ddx_h < star["x"].max()) * \
                           ((cell["y"] + ddx_h) > star["y"].min()) * (cell["y"] - ddx_h < star["y"].max()) * \
                           ((cell["z"] + ddx_h) > star["z"].min()) * (cell["z"] - ddx_h < star["z"].max()))[0]
                    
                    
    sub_cell = cell[ind_subcell]
    sub_cell_T = cell_T[ind_subcell]
    sub_cell_rho = cell_rho[ind_subcell]

    ############################################################

    sub_ddx_h = sub_cell["dx"] * 0.5
    sub_ddx_list = np.unique(sub_ddx_h)
    sddxsort = (sub_cell_T < Tcut) * (sub_cell_rho > rhocut)
    sub_ddx_h_max = max(sub_ddx_h[sddxsort])

    xl = sub_cell["x"] - sub_ddx_h_max*np.sqrt(3)
    xr = sub_cell["x"] + sub_ddx_h_max*np.sqrt(3)
    yl = sub_cell["y"] - sub_ddx_h_max*np.sqrt(3)
    yr = sub_cell["y"] + sub_ddx_h_max*np.sqrt(3)
    zl = sub_cell["z"] - sub_ddx_h_max*np.sqrt(3)
    zr = sub_cell["z"] + sub_ddx_h_max*np.sqrt(3)

    

    xrange = (xl.min(), xr.max())
    yrange = (yl.min(), yr.max())
    zrange = (zl.min(), zr.max())
    xspan = xrange[1] - xrange[0]
    yspan = yrange[1] - yrange[0]
    zspan = zrange[1] - zrange[0]

    ssort = (star["x"] >= xrange[0]) * (star["x"] <= xrange[1]) * (star["y"] >= yrange[0]) * (star["y"] <= yrange[1])
    x_real_span = max(star["x"][ssort]) - min(star["x"][ssort])
    y_real_span = max(star["y"][ssort]) - min(star["y"][ssort])
    z_real_span = max(star["z"][ssort]) - min(star["z"][ssort])

    #print("xspan = %s" %(xspan))
    #print("real xspan = %s" %(x_real_span))

    npixx = np.int(np.ceil(xspan / cell_dx_min))
    npixy = np.int(np.ceil(yspan / cell_dx_min))
    npixz = np.int(np.ceil(zspan / cell_dx_min))
    
    print("")
    print("---------------------")
    print("x-npix = %s" %(npixx))
    print("y-npix = %s" %(npixy))
    print("z-npix = %s" %(npixz))
    print("---------------------")
    print("")
    
    h = np.histogram2d(star["x"], star["y"],bins=[npixx, npixy],range=[xrange, yrange])

    sx = np.searchsorted(h[1], star["x"]) - 1
    sy = np.searchsorted(h[2], star["y"]) - 1

    cx = np.searchsorted(h[1], sub_cell["x"]) - 1
    cy = np.searchsorted(h[2], sub_cell["y"]) - 1

    zindex = np.linspace(zrange[0],zrange[1],npixz)
    sz = np.searchsorted(zindex, star["z"]) -1
    cz = np.searchsorted(zindex, sub_cell["z"]) -1


    x_padding_layer = int(max(rot_cell_length)/2)
    y_padding_layer = int(max(rot_cell_length)/2)
    z_padding_layer = int(max(rot_cell_length)/2)
    #print("x padding layer = %s" %(x_padding_layer) )
    #print("y padding layer = %s" %(y_padding_layer) )
    #print("z padding layer = %s" %(z_padding_layer) )
    #print("")

    Grid = np.zeros([npixz+z_padding_layer*2+2,npixy+y_padding_layer*2,npixx+x_padding_layer*2])

    from multiprocessing import Process, Lock, Value, Array #,shared_memory
    from shared_ndarray import SharedNDArray
    
    GRID = SharedNDArray((npixz+z_padding_layer*2+2,npixy+y_padding_layer*2,npixx+x_padding_layer*2))
    
    GRID = SharedNDArray.zeros_like(Grid)
    GRID.unlink()
    print("---------------------")
    print("Grid : Ready")

    rot_cell_length = np.array(rot_cell_length)
    rot_cell_centre_num = (rot_cell_length-1)/2 -1 +1 +2 ### -1 : indexing is start from zero, not one.
                                                         ### +2 : first and second layer is xy indexing layer.
                                                         ### +1 : length + center(1) + length

    nstar = len(star)
    colden_star = np.zeros(nstar)
    Colden = SharedNDArray.zeros_like(colden_star)
    #Colden = SharedNDArray.copy(colden_star)
    Colden.unlink()
    print("Colden : Ready")
    print("---------------------")
    print("")

    def repeat_grid(num_list,Grid):
        
        i = 0
        while i < len(num_list[0]):
            sc = int(num_list[0][i])

            rot_cell_index_sc = np.searchsorted(cell_ddx_list, sub_ddx_h[sc])-1 #sub_ddx_list or cell_ddx_list
            if rot_cell_index_sc == -1:
                rot_cell_index_sc = 0
            elif rot_cell_index_sc >= grid_number:
                rot_cell_index_sc = grid_number - 1

            rot_cell_length_sc = rot_cell_length[rot_cell_index_sc]
            rot_cell_centre_num_sc = int(rot_cell_centre_num[rot_cell_index_sc])
            rot_cell_sc = rot_cell[rot_cell_index_sc]

            Gx_sc = cx[sc] + x_padding_layer
            Gy_sc = cy[sc] + y_padding_layer
            Gz_sc = cz[sc] + z_padding_layer + 2  ### +2 for removing xy-indexing layer

            rcsc_low = rot_cell_centre_num_sc - int((rot_cell_length_sc - 1) / 2)
            rcsc_upp = rot_cell_centre_num_sc + int((rot_cell_length_sc - 1) / 2)

            Gx_sc_low = Gx_sc - int((rot_cell_length_sc - 1) / 2)
            Gx_sc_upp = Gx_sc + int((rot_cell_length_sc - 1) / 2)

            Gy_sc_low = Gy_sc - int((rot_cell_length_sc - 1) / 2)
            Gy_sc_upp = Gy_sc + int((rot_cell_length_sc - 1) / 2)

            Gz_sc_low = Gz_sc - int((rot_cell_length_sc - 1) / 2)
            Gz_sc_upp = Gz_sc + int((rot_cell_length_sc - 1) / 2)


            if sub_cell_T[sc] < Tcut and sub_cell_rho[sc] > rhocut:
                Grid[Gz_sc_low:Gz_sc_upp + 1, Gy_sc_low:Gy_sc_upp + 1, Gx_sc_low:Gx_sc_upp + 1] += rot_cell_sc[rcsc_low:rcsc_upp + 1] * sub_cell["var0"][sc] * (sub_cell["dx"][sc] ** 3) / (cell_dx_min ** 2)
                #print(sum(sum(sum(Grid))))
            else:
                pass
            i += 1
            if i == len(num_list[0]):
                GRID.array[:] += Grid[:]
                #print(sum(sum(sum(GRID.array))))
                #np.save("./tempo_grid%s.npy" %(num_list[1]+1),Grid)
                #print("%s process : done" %(num_list[1]+1))
        

    def repeat_colden(num_list,colden_star):
        j = 0
        while j < len(num_list):
            i = int(num_list[j])
            sx_i = sx[i]
            sy_i = sy[i]
            sz_i = sz[i]
            colden_star[i] = sum(Grid[z_padding_layer + 2:z_padding_layer + 2 + sz_i, y_padding_layer + sy_i, x_padding_layer + sx_i])
            #print("%s %s %s %s" %(sx_i,sy_i,sz_i,colden_star[i]))
            j += 1
            if j == len(num_list):
                Colden.array += colden_star
                #print("done")
                
                return

    N1 = int(len(sub_cell)/processing_number)
    N2 = int(len(star)/processing_number)

    Num_List1 = []
    Num_List2 = []
    Grid_List = []
    n = 1
    while n <= processing_number:
        Nlist_n1 = np.linspace(N1*(n-1),N1*n-1,N1)
        Nlist_n2 = np.linspace(N2*(n-1),N2*n-1,N2)
        Glist_n = np.zeros([npixz+z_padding_layer*2+2,npixy+y_padding_layer*2,npixx+x_padding_layer*2])
        #Num_List1.append(Nlist_n1)
        Num_List1.append([Nlist_n1,n-1])
        Num_List2.append(Nlist_n2)
        Grid_List.append(Glist_n)
        #print(Num_List1)
        n += 1

    rgs = []
        
    lock = Lock()
    for index, nums in enumerate(Num_List1):
        rg = Process(target=repeat_grid,args=(nums,Grid))
        rgs.append(rg)
        rg.start()

    for rg in rgs:
        rg.join()
    
    Grid += GRID.array

    print("---------------------")
    print("Grid : Done")
        

    rgs = []

    for index, nums in enumerate(Num_List2):
        rg = Process(target=repeat_colden,args=(nums,colden_star,))
        rgs.append(rg)
        rg.start()

    for rg in rgs:
        rg.join()

    colden_star += Colden.array
    colden_star = colden_star * gal.info.unit_nH * gal.info.unit_l 

    print("Colden : Done")
    print("---------------------")
    print("")

    return Grid, colden_star

