from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("rd_hal",
                        sources=["load_hal.pyx"],
                        libraries=["c_rd_halo"],
                        library_dirs=["/home/hoseung/Work/pyclusterevol/tree/"],
                        runtime_library_dirs=["/home/hoseung/Work/pyclusterevol/tree/"],
                         language='c++',)]
setup(cmdclass = {'build_ext': build_ext},
              ext_modules = ext_modules)

