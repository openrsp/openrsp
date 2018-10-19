

==================
Compiling the code
==================

Dalton on stallo
----------------

First we make sure that we have all required
modules loaded::

  $ module load python/2.7.3 intel/13.0 mkl/11.0.0 cmake/2.8.9 git/1.8.3.4 openmpi/1.6.2

We verify that mpif90 really resolves to Intel::

  $ mpif90 --version

Now we can configure the code::

  $ ./setup --mpi --mkl=parallel

And build it::

  $ cd build
  $ make dalton.x

BLAS/LAPACK errors
------------------

On some occasions (at least experienced when compiled sequentially on Ubuntu), there might be errors with missing references to BLAS or LAPACK libraries (the reason is not fully understood, but it could just be that the default version of these libraries on Ubuntu (14.04) is incomplete). This was resolved by following the steps by user 'treerink' on http://ubuntuforums.org/showthread.php?t=1505249​ ::

  "Go to: System -> Synaptic -> Administration -> Package Manager -> 
  search on lapack (and/or blas), and mark for installation:

  libblas3gf
  libblas-doc
  libblas-dev
  liblapack3gf
  liblapack-doc
  liblapack-dev

  -> Apply "

