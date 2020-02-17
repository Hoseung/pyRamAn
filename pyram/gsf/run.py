import sys, os, numexpr
import gsf

print('This script calls gsf given a valid simulation output file, snaphot_file (required input)')
print('gsf.py is a collection of functions that cover the computation of the gravitational potential,')
print('a wrapper for the Gaussian Mixture Models of scikit-learn and various plotting options.')
print('To have an idea of the various arguments that can be set use verbose=True.')
print('If no optional argument is given, the code will look for 2 clusters in the kinematic stellar space of')
print('(jz/jc,jp/jc,binding_energy) of the most massive halo in snaphot_file.')

print('The runtime for the computation of the potential scales as ~N^2')
print('If you want to run this part of gsf in parallel, set the system variable OMP_NUM_THREADS first.')
       
snaphot_file = '/scratch/database/nihao/gasoline2.1/g8.26e11/g8.26e11.01024'
out_dir='./'

print('Running gsf for the sim in %s'%snaphot_file)
gsf.wgsf(snaphot_file, out_dir=out_dir, number_of_clusters=2, covariance_type='full', whiten_data=None, 
                halo_id=1, radius_align=0.1, align_with='baryon', n_init=1, plot=True, verbose=True)