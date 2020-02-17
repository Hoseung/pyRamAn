# -*- coding: utf-8 -*
"""
Created on Thu Apr  9 15:35:43 2015

Generates a list of zoom-in target candidate.
General properties in csv format and
individual zoomin_parameters in .txt format.
Don't worry. Both are pickled.

@author: hoseung

512^3 particle data is 8.1GB.

"""

def part_2_den(part, info, region=None, proj='z', npix=800):
    """
    creates img object, calculate 2dprojection map and return the img object.
    given part, info, region, it passes x,y,z,m,npix,info arrays to pp.den2d.
    """
    import numpy as np
    from draw import img_obj, pp
    if region is None:
        region={"ranges":[[0,1]]*3, "center":[0.5]*3, "radius":0.5, \
        "xr":[0, 1], "yr":[0, 1], "zr":[0, 1]}

    ind_ok = np.where((part.x > region["xr"][0]) & (part.x < region["xr"][1]) & \
                        (part.y > region["yr"][0]) & (part.y < region["yr"][1]) & \
                        (part.z > region["zr"][0]) & (part.z < region["zr"][1]))

    if len(ind_ok[0]) > 0:
        img = img_obj.MapImg(info=info, proj='z', npix=npix, ptype=ptype)
        img.set_data(pp.den2d(part["x"][ind_ok], part["y"][ind_ok], part["z"][ind_ok], \
                    part["m"][ind_ok], npix, region=region, cic=True, norm_integer=True))
        img.set_region(region)
        return img
    else:
        return False

def get_pt_id_targets(part, h, halo_inds, r_frac=1):
    import numpy as np
    """
    Returns id list of particles in each halo in a SINGLE array.
    Also returned is the list of range of id list for each halo.

    So,
    idlist1 = idslist[ind[0] : ind[1]]
    idlist2 = idslist[ind[1] : ind[2]]

    .. note::
        You don't know how many particles will be found within halos you are interested.
        A good guess would be the sum of halo.np which is number of member particles
        defined by the halo finder.

        So, I allocate the idlist array with length of 1.2 * sum(halo.np[ind]).
    """
    assert(np.size(halo_inds) > 0)
    buffer_frac = 1.2 * np.sqrt(r_frac)
    # search radii larger than r_vir will result in more particles.
    npart_tot_guess = np.math.ceil(sum(h.data.np[halo_inds]) * buffer_frac)
    idlist = np.zeros(npart_tot_guess, dtype=np.int32)
    x = part.x
    y = part.y
    z = part.z
    print("Expected # of particles:", npart_tot_guess)
    nhalo = len(halo_inds)
    halo_range_list = np.zeros((nhalo) + 1, dtype=np.int32)
    len_idlist = len(idlist)
    for i, ihalo in enumerate(halo_inds):
        print("# part in this halo", h.data.np[ihalo])
        xr = [h.data.x[ihalo] - h.data.rvir[ihalo] * r_frac,
              h.data.x[ihalo] + h.data.rvir[ihalo] * r_frac]
        yr = [h.data.y[ihalo] - h.data.rvir[ihalo] * r_frac,
              h.data.y[ihalo] + h.data.rvir[ihalo] * r_frac]
        zr = [h.data.z[ihalo] - h.data.rvir[ihalo] * r_frac,
              h.data.z[ihalo] + h.data.rvir[ihalo] * r_frac]
        ind_x = np.where( (x > xr[0]) & (x < xr[1]))[0]
        # out of 8.1GB particle information()
        if len(ind_x) > 0:
            ind_y = np.where((y[ind_x] > yr[0]) & (y[ind_x] < yr[1]))[0]
            if len(ind_y) > 0:
                ind_z = np.where((z[ind_x[ind_y]] > zr[0]) & (z[ind_x[ind_y]] < zr[1]))[0]
                # If array is not large enough, append it.
                halo_range_list[i + 1] = halo_range_list[i] + len(ind_z)
                if halo_range_list[i + 1] > len_idlist:
                    if i + 1 == nhalo:
                        # If it's the last halo, append by exact difference.
                        npart_more = halo_range_list[i + 1] - len_idlist
                    else:
                        # Otherwise, guess from previous halos.
                        # 1.5 * mean npart so far * number of remaining halos
                        npart_more = int(1.5 * (halo_range_list[i] / (i + 1)) * (nhalo - i))
                    print("increase the array size by {:d} from {:d}".format(npart_more, len_idlist))
                    idlist = np.append(idlist, np.zeros(npart_more, dtype=np.int32))
                    len_idlist= len(idlist)

                idlist[halo_range_list[i]:halo_range_list[i+1]] = part.id[ind_x[ind_y[ind_z]]]
                print(halo_range_list[i+1], len(idlist))


    return halo_range_list, idlist

def particle_distribution(part, idlist, buffer = 0):
    from utils import match
    """
    returns center of mass of given particles, xr,yr,zr,
    and a,b,c (in the future)
    """
    ind = match.match_list_ind(part.id, idlist)
    print(" Number of matching particles: {} out of {}".format(len(ind), len(idlist)))
    xr = [part["x"][ind].min(), part["x"][ind].max()]
    yr = [part["y"][ind].min(), part["y"][ind].max()]
    zr = [part["z"][ind].min(), part["z"][ind].max()]
    radius = max([xr[1] - xr[0], yr[1] - yr[0], zr[1] - zr[0]])
    xc = sum(xr) * 0.5
    yc = sum(yr) * 0.5
    zc = sum(zr) * 0.5

    return xc, yc, zc, radius

# Load halo (final shanpshot)
import matplotlib.pyplot as plt
import numpy as np
import pickle

#import pyximport; pyximport.install(pyimport = True)
import load
import tree.halomodule as halo
import utils.sampling as smp
from draw import pp

#sys.stdout = open('stdout.txt', 'w')
#sys.stderr = open('stderr.txt', 'w')

npix = 800
work_dir = './'

decrease_region = False
show = False
denmap = True
nout_fi = 11
nout_ini = 1
nnout = nout_fi - nout_ini + 1
snout = str(nout_fi).zfill(3)
radius = 0.3
m_threshold = 1e14  # in solar mass (physical).

s = load.sim.Sim(nout=nout_fi, base=work_dir, dmo=True)
h = halo.Halo(base=work_dir, nout=nout_fi, info=s.info, halofinder="HM", load=True)

region = smp.set_region(centers=[0.5]*3, radius=radius)
s.set_ranges(ranges=region["ranges"])
s.show_ranges()

ptype = ["dm mass id pos"]
s.add_part(ptype, load=True)
part = s.part.dm
print("M ptp", part['m'].ptp())
part["m"] = part["m"] *s.info.msun

pp_halo_rscale = 10
img = part_2_den(s.part.dm, s.info, region=region, npix=npix)
#field = img.data

fname = work_dir + \
 '{}_{}_map_{}.pickle'.format(snout, ptype[0][0], "zr")
with open(fname, 'wb') as f:
    pickle.dump(img, f)

#%%
fig = plt.figure()
ax1 = fig.add_subplot(111)

#img = draw.img_obj.MapImg(info=s.info, proj='z', npix=npix, ptype=ptype)
#img.set_data(draw.pp.den2d(x, y, z, m, npix, s.info, cic=True, norm_integer=True))
img.plot_2d_den(axes=ax1, show=show, vmax=3e13, dpi=200, cname='brg', zposition=True)
#pp.plot_part_2d_den(field, show=show, vmin=1e7, vmax=3e13, dpi=500, axes=ax1, cname='brg', zposition=True)

# halo over density map

# Because pp.circle_scatter draws fixed size objects,
# you need to explicitly tune the sizes of two different plots to have the same size.
# i.e., npix is needed. (physical size of two plots are inferred from the region parameter)
# OR, you can just use plt.scatter
#
# img = fig.imshow()
# plt.scatter()
#
# Because the current axes is  maintained,
# plt.scatter automatically plots over img object.
# The downside is that the propotional size of circles to density map
# will not be pertained as you zoom in or out the image.
#
#utils.util.reimport(pp)
h_ind = np.where( (h.data.x > region["xr"][0]) & (h.data.x < region["xr"][1]) &
                (h.data.y > region["yr"][0]) & (h.data.y < region["yr"][1]) &
                (h.data.z > region["zr"][0]) & (h.data.z < region["zr"][1]) &
                (h.data.mvir > m_threshold))[0]

pp.pp_halo(h, npix, ind=h_ind, region=region, rscale=pp_halo_rscale, \
            axes=ax1, facecolor='none', edgecolor='b')

plt.savefig(work_dir + snout + "target_map.png")
# Conclusion:
# position of halos and DM particles seems to match well.

# select halos (mass cut, rank cut, and some more ,in the future.)
# opitons are: Mmin, Mmax, margin (or xr,yr,zr)

#%%
# load on momery particle IDs of all halos
range_list, id_list = get_pt_id_targets(s.part.dm, h, h_ind, r_frac=7.0)

print("++++++++++++ \n Loop over snapshots \n +++++++++++++ \n")
# Draw particle map to check.
import utils.sampling as smp

if denmap:
    scale = 8.0
    for i in h_ind:
        h_region = smp.set_region(xc=h.data.x[i], yc=h.data.y[i], zc=h.data.z[i],
                                radius = h.data.rvir[i] * scale)
    #    print(h_region)
        img = part_2_den(part, s.info, region=h_region, npix=npix)
        if img is False:
            continue
    
        sid = str(h.data.id[i]).zfill(5)
        f_save = work_dir + snout + ptype[0][0] + "halo_" + sid + img.proj + ".pickle"
        img.pickle_data(fout=f_save)
    #    with open(f_save, 'wb') as f:
    #        img = pickle.load(f)
        fout = work_dir + snout + ptype[0][0] + "halo_" + sid + img.proj + ".png"
        img.plot_2d_den(save=fout, show=show, vmin=1e7, vmax=3e10, dpi=100, zposition=True)

scale = 8.0
for i in h_ind:
    h_region = smp.set_region(xc=h.data.x[i], yc=h.data.y[i], zc=h.data.z[i],
                            radius = h.data.rvir[i] * scale)
    img = part_2_den(part, s.info, region=h_region, npix=npix)
    if img is False:
        continue

    sid = str(h.data.id[i]).zfill(5)
    f_save = work_dir + snout + ptype[0][0] + "halo_" + sid + img.proj + ".pickle"
    img.pickle_data(fout=f_save)
    fout = work_dir + snout + ptype[0][0] + "halo_" + sid + img.proj + ".png"
    img.plot_2d_den(save=fout, show=show, vmin=1e7, vmax=3e10, dpi=100, zposition=True)

# OK, it looks good! (04.12)

# Save final properties (ID, position, mass, rvir, #DM particles)
newh = halo.Halo()
newh.derive_from(h, h_ind)
# First, pickle it.
f_halo = "halo_sub.pickle"
with open(work_dir + f_halo, 'wb') as f:
    pickle.dump(newh, f)


#%%
refine_params = np.recarray((newh.data.shape[0], nnout),
                            dtype=[("x_refine", "f8"),
                                   ("y_refine", "f8"),
                                   ("z_refine", "f8"),
                                   ("r_refine", "f8")])
# for each nout
#import utils.sampling as smp
aexps=[]
for inout, thisnout in enumerate(range(nout_fi, nout_ini -1, -1)):
#    snout = str(nout).zfill(3)
    print(inout, thisnout)

    s = load.sim.Sim(nout=thisnout, base=work_dir, dmo=True, setup=True)
    if decrease_region:
        region = smp.set_region(centers=[0.5]*3, radius=radius)
        s.set_ranges(ranges=region["ranges"])

    ptype = ["dm pos id"]
    s.add_part(ptype, dmmass=False)
    s.part.load(verbose=False)
    

    for ihal, thishalo in enumerate(newh.data.id):
        xr, yr, zr, rr = particle_distribution(s.part.dm,
                            id_list[range_list[ihal]:range_list[ihal + 1]])
        refine_params[ihal][inout] = np.asarray([xr, yr, zr, rr])

    aexps.append(s.info.aexp)

    if decrease_region:
        radius = min([radius * 1.1, radius + 0.01, 0.5])
        print("\n RADIUS:", radius)


class Allrefine():
    def __init__(self, hid, data, aexp):
        self.id = hid
        self.x = data['x_refine'][::-1]
        self.y = data['y_refine'][::-1]
        self.z = data['z_refine'][::-1]
        self.r = data['r_refine'][::-1]
        self.aexp = aexp[::-1]
        

allrf = Allrefine(newh.data.id, refine_params, aexps)

with open(work_dir + "refine_params.pickle", 'wb') as f:
    pickle.dump(allrf, f)
