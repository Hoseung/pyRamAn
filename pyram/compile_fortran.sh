#/bin/bash

hspy
cd gsf
./setup

cd ../tree
f2py -c tmtree.f90 -m cnt_tree
f2py -c tmtree_dp.f90 -m cnt_tree_dp

cd ../load
f2py -c amr2cell_fun.f90 -m a2c
f2py -c part_cpu_module.f90 -m part_load
f2py -c readr.f90 -m readr
