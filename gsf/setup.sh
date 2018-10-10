#!/bin/bash

## Compile the fortran module with omp to get the shared lib twobody.so and the interface twobody.pyf
##

f2py -m twobody --fcompiler=gfortran --f90flags='-fopenmp' -lgomp -c twobody.f95

f2py -h twobody.pyf -m twobody --overwrite-signature twobody.f95

## A few lines in twobody.pyf need a bit of editing, which is done by edit_pyf.py

python edit_pyf.py
