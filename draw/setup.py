from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize


extensions=[Extension("ppc",["ppc.pyx"])]
setup(
    ext_modules = cythonize(extensions)
)
