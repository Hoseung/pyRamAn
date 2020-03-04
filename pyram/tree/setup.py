from setuptools import setup, Extension
#from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

ext_modules = [Extension("rd_hal",
                        sources=["load_hal.pyx"],
                        libraries=["c_rd_halo"],
                        library_dirs=["./"],
                        runtime_library_dirs=["./"],
                        language='c++',)]
setup(cmdclass = {'build_ext': build_ext},
              ext_modules = cythonize(ext_modules))

