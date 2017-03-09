f2py -c --fcompiler=gnu95 -m a2c amr2cell_fun.f90
f2py -c --fcompiler=gnu95 -m part_shared part_cpu.f90
f2py -c --fcompiler=gnu95 --build-dir ../tree -m cnt_tree ../tree/tmtree.f90
mv cnt_tree.cpython* ../tree

