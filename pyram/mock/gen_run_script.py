import pickle
import numpy as np

galid_fi = 1
nout_fi = 925
gtree = pickle.load(open(f"gtree_{galid_fi}_{nout_fi}.pickle", "rb"))

all_nouts = gtree['nout']

nchunk = 4
for cnt_script in range(nchunk):
    nouts_chunk = gtree['nout'][cnt_script:,:,nchunk]
    fname_script = out_dir + f"run{cnt_script}.sh"
    with open(fname_script, 'w') as f:
        f.write("""#!/bin/sh 
        #PBS -S /bin/sh
        #PBS -N SK867_1
        #PBS -j oe
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
        export PATH="${HOME}/SKIRT/release/SKIRT/main:${PATH}" """)
        for no in nouts_chunk:
            f.wrtie(f"cd /scratchb01/hoseung/gal_{nout:05d}_{galid_if}/{galid_fi:05d}/{nout:05d}/ && skirt g*_*.ski")

