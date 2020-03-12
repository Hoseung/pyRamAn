from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("rd_hal", ["c_rd_halo.cpp", "load_hal.pyx"],
          language='c++',)]
setup(cmdclass = {'build_ext': build_ext},
              ext_modules = ext_modules)

