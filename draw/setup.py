from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

extensions=[Extension("ppc",["ppc.pyx"])]
setup(
    ext_modules = cythonize(extensions),
    include_dirs=[numpy.get_include()]
)
