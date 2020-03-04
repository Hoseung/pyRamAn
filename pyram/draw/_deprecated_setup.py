from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

extensions=[Extension("ppc",["ppc.pyx"])]
setup(
    ext_modules = cythonize(extensions),
    include_dirs = [np.get_include()],
	extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-fopenmp" ],
    extra_link_args=['-fopenmp']
)
