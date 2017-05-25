f2py -c --fcompiler=gnu95 -m a2c amr2cell_fun.f90
f2py -c --fcompiler=gnu95 -m part_shared part_cpu.f90
cd ../tree
f2py -c --fcompiler=gnu95 -m cnt_tree tmtree.f90

