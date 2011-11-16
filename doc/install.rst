

Installation
============

XCint is configured using CMake, typically via the setup script,
and subsequently compiled using make or gmake.


Dependencies
------------

The configuration/compilation step requires:
CMake  >= 2.6.3

The configuration/compilation step using setup script requires:
CMake  >= 2.6.3
python >= 2.4.0


What to do if CMake is not available or too old?
------------------------------------------------

It is your machine and you have Ubuntu or Debian::

  $ sudo apt-get install cmake cmake-curses-gui cmake-qt-gui

Similar mechanisms exist for other distributions or
operating systems. Please consult Google.

If it is a cluster, please ask the Administrator to install/upgrade CMake.

If it is a cluster, but you prefer to install it yourself (it's easy):

1. Download the latest tarball from http://www.cmake.org/cmake/resources/software.html
2. Extract the tarball to, say, ~/cmake-2.8.5-Linux-i386

Now you can enjoy your CMake instead of the old one::

  $ ./setup [other flags] --cmake=~/cmake-2.8.5-Linux-i386/bin/cmake


Installation for normal people
------------------------------

XCint is configured using CMake.
The setup script is a useful front-end to CMake.
It does nothing else than creating the directory "build" and calling
CMake with appropriate environment variables and flags::

  $ ./setup [--flags]
  $ cd build
  $ make

Call setup without flags to see all available options::

  $ ./setup

You can see the CMake command using::

  $ ./setup [--flags] --show

and use it directly to call CMake without setup.


Most typical examples
---------------------

In order to get familiar with the new configuration setup, let us demonstrate
some of the most typical configuration scenarios.

Configure for parallel compilation using MPI (make sure to properly export MPI
paths)::

  $ ./setup --fc=mpif90 --cc=mpicc

Configure for parallel compilation using MPI and 64-bit integers::

  $ ./setup --fc=mpif90 --cc=mpicc --int64

Configure for sequential compilation using ifort/icc::

  $ ./setup --fc=ifort --cc=icc

Configure for sequential compilation using gfortran/gcc::

  $ ./setup --fc=gfortran --cc=gcc

You get the idea. The configuration is pretty good at detecting math libraries
automatically, provided you export proper paths, see "Linking to external math
libraries".


Out of source compilation
-------------------------

By default CMake builds out of source.  This means that all object files and
the final binary are generated outside of the source root.  Typically the build
directory is called "build", but you can modify this default behavior using the
--build flag. This strategy offers several advantages. One obvious advantage is
that you can now build several binaries with the same source::

  $ cd /sourcepath
  $ ./setup --fc=gfortran --cc=gcc --build=/gfortran-buildpath
  $ cd /gfortran-buildpath
  $ make
  $ cd /sourcepath
  $ ./setup --fc=ifort --cc=icc --build=/ifort-buildpath
  $ cd /ifort-buildpath
  $ make


Basic installation without the setup script
-------------------------------------------

If you are familiar with CMake you don't have to use the setup script.
The setup script does nothing else than calling CMake with appropriate
environment variables and flags, it is a convenient front-end.

Minimal example::

  $ mkdir build
  $ cd build
  $ cmake ..
  $ make

The two following strategies are completely
equivalent:

Using CMake directly::

  $ mkdir build
  $ cd build
  $ FC=mpif90 CC=mpicc cmake -DENABLE_MPI=1 -DCMAKE_BUILD_TYPE=Release ..
  $ make

Using setup::

  $ ./setup --fc=mpif90 --cc=mpicc (--mpi)
  $ cd build
  $ make

If the compiler contains "mpi", then you can omit the flag --mpi, setup will set
it in this case automatically.

Please note that the defaults for performance optimization are different for
setup and direct CMake: by default setup configures for optimization, whereas
direct CMake commands configure code without optimization. Both defaults can be
changed.

There is nothing special about the directory "build".
You can do this instead::

  $ mkdir /buildpath
  $ cd /buildpath
  $ cmake /sourcepath
  $ make


Linking to external math libraries
----------------------------------

Typically you will want to link to external math (BLAS and LAPACK) libraries,
for instance provided by MKL or Atlas.

The CMake configuration script will automatically find them if you define MATH_ROOT::

  $ export MATH_ROOT='/opt/intel/mkl'

Do not use full path MATH_ROOT='/opt/intel/mkl/lib/ia32'. CMake will append the
correct paths depending on the processor and the default integer type.  If the
MKL libraries that you want to use reside in
/opt/intel/mkl/10.0.3.020/lib/em64t, then MATH_ROOT is defined as::

  $ export MATH_ROOT='/opt/intel/mkl/10.0.3.020'

Then::

  $ ./setup [--flags]                 # do not need to define --math
  $ cd build
  $ make

Alternatively::

  $ cd build
  $ [FC=gfortran CC=gcc] MATH_ROOT='/opt/intel/mkl' cmake ..
  $ make

Exporting MATH_ROOT is equivalent to calling setup with --math-dir::

  $ ./setup --math-dir=/opt/intel/mkl

If automatic detection of math libraries fails for whatever reason, you can
always call the libraries explicitly like here::

  $ ./setup --math="-L/path -lfoo -lbar"


Running CMake using GUI
-----------------------

You prefer GUI? No problem. You can configure with GUI::

  $ cd build
  $ cmake ..
  $ cmake-gui ..

You may have to install cmake-gui for it, on debian/ubuntu::

  $ sudo apt-get install cmake cmake-curses-gui cmake-qt-gui


Running tests
-------------

You can run the test suite with::

  $ make test

It is HIGHLY recommended to run the test set after you have compiled
XCint to make sure that your binary correctly reproduces reference results.


Make install
------------

Make install is very useful to make XCint available to other users on the same
machine::

  $ ./setup [--flags] --install=/path
  $ cd build
  $ make
  $ make install


Where should $PATH point to? Source directory or build directory?
-----------------------------------------------------------------

We recommend to let $PATH point to the install directory::

  $ ./setup [--flags] --install=/install/path
  $ cd build
  $ make
  $ make install

This way everything (binary, scripts, basis sets) will be at the right place
under /install/path and $PATH should contain /install/path.


Compiling in verbose mode
-------------------------

Sometimes you want to see the actual compiler flags and definitions::

  $ make VERBOSE=1


Compiling on many cores
-----------------------

Yes, it works. Try::

  $ make -j4


How can I change optimization flags?
------------------------------------

If you want to turn optimization off (debug mode), there are several ways to do that.

Either use setup::

  $ ./setup --debug [other flags]
  $ cd build
  $ make

Or use Cmake directly (default here is debug mode)::

  $ mkdir build
  $ cd build
  $ [FC=ifort CC=icc] cmake ..
  $ make

If you want to modify compiler flags, edit cmake/FCompilers.cmake and/or cmake/CCompilers.cmake.
