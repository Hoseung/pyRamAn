import os
import numpy as np
import pickle
import pyram
from pyram.tree import halomodule as hmo 
from pyram.mock import make_ramski
from pyram.utils import match
from pyram.utils.cosmology import Timeconvert as TC

def get_contam_fname(contam_dir, nout, purity="100"):
    cdir = os.path.join(contam_dir, purity+"PerCentPurity/")
    fname = cdir + f"list_gal_{nout:05d}.dat_nocontam"
    return fname

def get_contam(contam_dir, nout, gids):
    contam = np.genfromtxt(get_contam_fname(contam_dir, nout, purity="100"), 
                dtype=[("id", int), ("level", int), ("mass",float), 
                       ("x", float), ("y", float), ("z",float)])
    try:
        len(gids)
        inds = match.match_list_ind(contam['id'], gids, allow_swap=False) 
    except:
        inds = np.where(contam['id'] == gids)[0]
    return contam[inds]

def replace_name(arr, ind, name):
    temp = list(arr.dtype.names)
    temp[ind] = name
    arr.dtype.names = tuple(temp)

make_scripts=True
make_tree = True

fn_template = './template_zubko.ski'
nout_fi=868 # It's 867...
out_dir = f'gal_{nout_fi:05d}/'
if not os.path.isdir(out_dir): os.mkdir(out_dir)

sim_dir = './'
r_multiply = 2 # FoV = 120% of R90 
contam_dir = './contam/'

if make_tree:
    tree = np.load("MJ_tree/morp_1.npy" )
    galid_fi = tree['galid'][0]
    gtree = np.zeros(len(tree), dtype=[('nout', int), ("id", int), ("level", int), ("mass",float),
                           ("r", float), ("x", float), ("y", float), ("z",float)])

    for i, tt in enumerate(tree):
        nout = tt['nout']
        galid = tt['galid']
        contam = get_contam(contam_dir, nout, galid)
        for nn in contam.dtype.names:
            gtree[i][nn] = contam[nn]
        gtree[i]['r'] = tt['R90'] * r_multiply
        gtree[i]['nout'] = nout

    pickle.dump(gtree, open(f"gtree_{galid_fi}_{nout_fi}.pickle", "wb"))
else:
    galid_fi = int(input("input galid_fi:"))
    gtree = pickle.load(open(f"gtree_{galid_fi}_{nout_fi}.pickle", "rb"))

ftext_name = out_dir+f"centers_{galid_fi}_{nout_fi}.txt"
# Check for partial outputs
# 
try:
    fsaved = np.genfromtxt(ftext_name)
    nout_last = int(fsaved[-1][0])
except:
    nout_last = nout_fi

#nout_last = 409
print("Nout of the last output :", nout_last)
print("Continuing from there")


# For testing purpose
#gtree = gtree[::100]


#####################################################
# Prepare scripts for all outputs first, and start running SKIRT from the last a few snapshots.
# Dont' wait for ~a day for all the outputs are ready.
if make_scripts:
    all_nouts = gtree['nout']

    nchunk = 24 
    for cnt_script in range(nchunk):
        nouts_chunk = gtree['nout'][cnt_script::nchunk]
        fname_script = out_dir + f"run{cnt_script}.sh"
        with open(fname_script, 'w') as f:
            f.write("""#!/bin/sh 
    #PBS -S /bin/sh\n""")
            f.write(f"#PBS -N SK{galid_fi}_{cnt_script}\n")
            f.write("""#PBS -j oe
    #PBS -l nodes=1:ppn=24,walltime=24:00:00

    module () {
      eval $(/usr/bin/modulecmd bash $*)
    }

    module purge
    module load gcc/8.3.0 cmake/3.14.7 git 
    export CCX=g++
    export CC=gcc
    # Then, you can build SKIRT.
    export OMP_NUM_THREADS=24
    export PATH="${HOME}/SKIRT/release/SKIRT/main:${PATH}" \n""")
            for no in nouts_chunk:
                f.write(f"cd /scratchb01/hoseung/gal_{nout_fi:05d}/{galid_fi:05d}/{no:05d}/ && skirt g*_*.ski\n")



#################################
# ver 1
r_fixed = 25
nvec_frac=0 # Don't calculate rotation vector of the galaxy
# Make a function to load it from the last gxx.ski file or take it as a return from make_ram_input
ski_params =[dict(name='face_on',  pixel_scale=50, nphoton=1e8, incli_off = 160.23863,  azim_off = 52.22508, roll_off = 0),
             dict(name='d30',      pixel_scale=50, nphoton=1e8, incli_off = 130.23863,  azim_off = 52.22508, roll_off = 0),
             dict(name='d60',      pixel_scale=50, nphoton=1e8, incli_off = 100.23863,  azim_off = 52.22508, roll_off = 0),
             dict(name='edge_on',  pixel_scale=50, nphoton=1e8, incli_off = 70.23863,  azim_off = 52.22508, roll_off = 0),
             dict(name='z',        pixel_scale=50, nphoton=1e8, incli_off = 0, azim_off = 0, roll_off = 0),
             dict(name='y',        pixel_scale=50, nphoton=1e8, incli_off = 0, azim_off = 0, roll_off = 90)]


if True:
    with open(ftext_name, 'a+') as fsave:
        for gdat in gtree:
            print(gdat)
            nout = gdat['nout']
            if nout > nout_last:
                continue
            # x,y,z in code unit, mass in Msun, r in kpc
            s=pyram.load.sim.Sim(nout, base=sim_dir, sim_type='nh')
            tc = TC(s.info)
            make_ramski.make_ram_input(s, tc, gdat, fsave,
                                           dir_out=out_dir,
                                           fn_template=fn_template,
                                           ski_params=ski_params,
                                           smooth = 50, 
                                           nvec_frac=nvec_frac,
                                           plot_stellar=True,
                                           more_props=True,
                                           r_fixed=r_fixed)
