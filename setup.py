import setuptools

from distutils.core import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy
from setuptools import find_packages

"""
 packages=setuptools.find_packages()
 -> list of packages needed in the package.
    find_packages() automatically finds all used packages.

 classifier: according to pypi classification scheme.
 https://pypi.org/classifiers/

"""

def main():
    with open("README.md", "r") as fh:
        long_description = fh.read()

    ext_modules = [Extension("pyram.draw.ppc",["pyram/draw/ppc.pyx"],
                            include_dirs=[numpy.get_include()])]
                #Extension( name='pyram/load/part_load',
                #        sources= ['pyram/load/part_cpu_module.f90'] )]

    setup(
        name="pyram",
        version="0.0.1",
        author="Hoesung Choi",
        author_email="hopung@gmail.com",
        description="Testing",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="github.com",
        packages=find_packages(),
        classifiers=[
            "Programming Language :: Python 3",
            "License ::OSI Approved :: MIT License",
            "Operating System :: Linux"
        ],
        python_requires='>=3.7',
        ext_modules = cythonize(ext_modules),
        #extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-fopenmp" ],
        #extra_link_args=['-fopenmp']
    )
    #setup(cmdclass = {'build_ext': build_ext},

    

        

if __name__ == "__main__":
    main()
    
#from numpy import distutils
#from numpy.distutils.core import setup
#from numpy.distutils.extension import Extension
# Fortran extension
#distutils.core.setup(cmdclass = {'build_ext': build_ext},
#        ext_modules = ([distutils.extension.Extension( 'pyram/load/part_cpu', ['pyram/load/part_cpu_module.f90'])]))
    
