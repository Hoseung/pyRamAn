#!/bin/bash

# This allows Python to access Fortran's internal memory!
gfortran -x f95-cpp-input -c readr.f90
f2py -c readr.f90 -m readr --opt='-O3 -x f95-cpp-input'

gfortran -x f95-cpp-input -c amr2cell_rur.f90
f2py -c amr2cell_rur.f90 -m a2c --opt='-O3 -x f95-cpp-input'

