#!/bin/bash

# in the hspy env.


# https://stackoverflow.com/questions/54201961/how-to-tell-f2py-module-to-look-in-current-directory-for-shared-object-dependenc
#
# Although you cannot pass linker flags from the command line in f2py, it will read the LDFLAGS environment variable.
# However, the default behavior for numpy is to overwrite the flags used in compiling rather than appending them which will cause failure in compiling if the required flags are not present in LDFLAGS.
# Support was added in numpy version 1.16.0 for optionally appending these linker flags by setting the environment variable NPY_DISTUTILS_APPEND_FLAGS=1
#

export LDFLAGS=-Wl,-rpath=.
export NPY_DISTUTILS_APPEND_FLAGS=1

gfortran -shared -fPIC -o ramski_module.so ramski_module.f90

python -m numpy.f2py -c ramski_gas.f90 ramski_module.so -m ramski
#f2py -c -m ramski ramski_gas.f90 ramski_module.so
