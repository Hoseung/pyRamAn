from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("rd_hal", ["c_utils_caller.pyx", "c_rd_halo.cpp"],
          language='c++',)]
setup(cmdclass = {'build_ext': build_ext},
              ext_modules = ext_modules)

