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

  ** There are a few scripts using h5py, or vispy.
   But let's just ignore them. 


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


Loading simulation raw data
^^^^^^^^^^^^^^^^^^^^^^^^^^^
By default, wdir='./' or base='./' for all methods/funcionts. 
So it is convenient to start analysis under a cluster directory (for example, /data1/good/01605/)::

   >>> import load
   >>> s = load.sim.Sim(base='./', nout=187)
   >>> s.add_part(ptypes=["star id pos vel", "dm id pos"])
   >>> s.add_hydro()
   >>> s.part.star['x']
   array([ 0.47147833,  0.48642134,  0.49925411, ...,  0.52071744,
        0.51837657,  0.50797539])
   >>> s.hydro.cell['var0'] 
   array([ 0.03143702,  0.04156503,  0.15768691, ...,  0.49851362,
        0.06658726,  0.0774503 ])
   
   * No real field names are given to the hydro variables.

Units
^^^^^
To be updated.


Loading Halo/GalaxyMaker catalogues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Each cluster simulation has it's halo/galaxy catalog at each snapshot.
I use HaloMaker/GalaxyMaker for structure finding, and data are already there.
To load the cataloguse::

   >>> import tree.halomodule as hmo
   >>> hals = hmo.Halo(nout=187, is_gal=False)
   >>> gals = hmo.Halo(nout=187, is_gal=True)
   >>> gals.data['mvir']
   array([  4.22624002e+11,   2.18115220e+07,   3.20211712e+08,
         4.12808602e+10,   6.30470042e+10,   4.08385952e+08,


Units
^^^^^
To be updated.



Loading Halo/Galaxy tree
^^^^^^^^^^^^^^^^^^^^^^^^
...



Basic plottings
^^^^^^^^^^^^^^^
particle density
hydro variables
halo 
merger tree...
