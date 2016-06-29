Getting Started
============================================================================

Prerequisites
---------------------------------------------------------------------------
Some python packages are needed.

- I recommend using virtualenv. ::
  
  $ virtualenv -p  python3 venv
  $ source venv/bin/activate


- Then install packages ::

  (venv) $ pip3 install numpy scipy matplotlib pandas astropy pyfits ipython

  ** There are a few scripts using h5py, or vispy. But let's just ignore them. 

- Manual compilation

  Small part of codes are written in Fortran, Cython, or C++. 
  I don't know how to automatically compile them, yet. ::

   $ load/compile_f2py.sh
   $ load/compile_part.sh
   $ cd draw
   $ python3 setup.py build_ext --inplace
   $ cd ../tree/load_c/
   $ python3 setup.py build_ext --inplace
   $ cd ../uitls/
   $ python3 setup.py build_ext --inplace

