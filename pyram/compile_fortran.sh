#/bin/bash

cd gsf
./setup

cd ../tree
f2py -c tmtree.f90 -m cnt_tree
f2py -c tmtree_dp.f90 -m cnt_tree_dp
f2py -c readhtm.f90 -m readhtm

cd ../load
f2py -c amr2cell_fun.f90 -m a2c
f2py -c part_cpu_module.f90 -m part_load
gfortran -x f95-cpp-input -c readr.f90
f2py -c readr.f90 -m readr --opt='-O3 -x f95-cpp-input'

