cd load
python setup.py build_ext --inplace 
f2py -c amr2cell_fun.f90 -m a2c
f2py -c part_cpu_module.f90 -m part_load
cd ../tree
python setup.py build_ext --inplace 
f2py -c tmtree_d.f90 -m cnt_tree
cd ../draw
python setup.py build_ext --inplace 
