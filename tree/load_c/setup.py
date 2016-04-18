from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("rd_hal", ["load_hal.pyx", "c_rd_halo.cpp"],
          language='c++',)]
setup(cmdclass = {'build_ext': build_ext},
              ext_modules = ext_modules)

