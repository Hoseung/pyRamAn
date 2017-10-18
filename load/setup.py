from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

# compile the fortran modules without linking
fortran_mod_comp = 'gfortran part_cpu_module.f90 -c -o part_cpu_module.o -O3 -fPIC'
#print fortran_mod_comp
system(fortran_mod_comp)
shared_obj_comp = 'gfortran py_part_cpu.f90 -c -o py_part_cpu.o -O3 -fPIC'

system(shared_obj_comp)


ext_modules = [Extension(# module name:
                         'py_part_cpu',
                         # source file:
                         ['py_part_cpu.pyx'],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3'],
                         # other files to link to
                         extra_link_args=['part_cpu_module.o', 'py_part_cpu.o'])]

setup(name = 'py_part_cpu',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)
