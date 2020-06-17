from distutils.core import setup, Extension
from Cython.Build import cythonize


ext = Extension(name='rd_hal', sources=["rd_hal.pyx"], language="c++")

setup(ext_modules=cythonize(ext))
