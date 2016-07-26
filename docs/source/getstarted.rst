Getting Started
============================================================================

Prerequisites
---------------------------------------------------------------------------
Some python packages are needed.

I recommend using virtualenv. If there is no virtualenv, install it by::

   $ pip3 install virtualenv

Then, ::
  
   $ virtualenv -p python3 --always-copy venv
   $ source venv/bin/activate


Install packages ::

  (venv) $ pip3 install numpy scipy matplotlib pandas astropy pyfits ipython cython

   There are a few scripts using h5py, or vispy.
   But let's just ignore them. 
   sphinx, numpydoc and autodoc are needed to build this document.


Manual compilation

  Some parts of the code are written in Fortran, Cython, or C++. 
  I don't know how to automatically compile them, yet.
  To compile external funtions::

   $ load/compile_f2py.sh
   $ load/compile_part.sh
   $ cd draw
   $ python3 setup.py build_ext --inplace
   $ cd ../tree/load_c/
   $ python3 setup.py build_ext --inplace
   $ cd ../uitls/
   $ python3 setup.py build_ext --inplace


