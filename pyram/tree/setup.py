from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import Cythonize
#from Cython.Distutils import build_ext
import numpy as np

ext_modules = [Extension("rd_hal",
                        sources=["rd_hal.pyx"],
                        libraries=["c_rd_halo"],
                        library_dirs=["./"],
                        runtime_library_dirs=["./"],
                         language='c++',)]
setup(ext_modules = Cythonize.cythonize(ext_modules),
      include_dirs = [np.get_include()],)

