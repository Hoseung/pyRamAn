from setuptools import setup, Extension, find_packages
from Cython.Distutils import build_ext
from Cython.Build import cythonize


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


    ext_modules = [Extension("pyram/tree/rd_hal",
                            sources=["pyram/tree/c_rd_halo.cpp",
                             "pyram/tree/load_hal.pyx"],
                            libraries=["pyram/tree/c_rd_halo"],
    #                        library_dirs=[path.join(here, 'tree/')],
    #                        runtime_library_dirs=[path.join(here, 'tree/')],
                            language='c++',)]

    # Fortran extension
    #ext_modules.append(Extension( 'load/part_cpu', ['load/part_cpu.f90'] ))

    #ext_modules = [Extension("pyram/tree/rd_hal",
    #["pyram/tree/c_rd_halo.cpp", "pyram/tree/load_hal.pyx"],
    #        language='c++',)]

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
        ext_modules = cythonize(ext_modules)
    )
    #setup(cmdclass = {'build_ext': build_ext},
        

if __name__ == "__main__":
    main()
