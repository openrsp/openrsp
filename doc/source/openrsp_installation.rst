.. _chapter-installation:

Installation of OpenRSP Library
===============================

Before installing OpenRSP, you need to make sure the following programs are
installed on your computer:

#. Git,

#. CMake (:math:`\ge2.8`),

#. C, C++ (if C++ APIs built) and/or Fortran 2003 (if Fortran APIs built) compilers,

#. HDF 5 (:math:`\ge1.8`) if it is enabled in
   `QcMatrix library <https://gitlab.com/bingao/qcmatrix>`_,

#. BLAS and LAPACK libraries.

CMake
-----

For the time being, only CMake could be used to compile OpenRSP. In general,
OpenRSP should be compiled together with the host programs. See for example
the Dalton program (``qcmatrix`` branch).

You could also compile OpenRSP alone to be familiar with how it works. But
no real calculations will be performed, all the callback functions in the
OpenRSP unit testing only return pre-defined data or read data from file.
Let us assume that you want to compile the library in directory ``build``,
you could invoke the following commands::

    mkdir build
    cd build
    ccmake ..
    make

During the step ``ccmake``, you need to set some parameters appropriately
for the compilation. For instance, if you enable ``OPENRSP_TEST_EXECUTABLE``, some
executables for the test suite will be built and can run after compilation. So
that you are able to check if OpenRSP has been successfully compiled. A detailed
list of the parameters controlling the compilation is given in the following table:

.. tabularcolumns:: |L|p{0.5\textwidth}|L|
.. list-table:: OpenRSP CMake parameters
   :header-rows: 1

   * - CMake parameters
     - Description
     - Default
   * - ``OPENRSP_USER_CONTEXT``
     - Enable user context in callback functions.
     - ``OFF``
   * - ``LIB_OPENRSP_NAME``
     - Sets the name of the OpenRSP library.
     - ``None``
   * - ``OPENRSP_CXX_API``
     - Build C++ API.
     - ``OFF``
   * - ``OPENRSP_F03_API``
     - Build Fortran 2003 API.
     - ``OFF``
   * - ``OPENRSP_TEST_EXECUTABLE``
     - Build the test suite excutables.
     - ``ON``
   * - ``QCMATRIX_LIB``
     - Sets the QcMatrix library (``xx/libqcmatrix.a``).
     - ``None``
   * - ``QCMATRIX_INCLUDE_DIR``
     - Sets the include directory of QcMatrix library.
     - ``None``
   * - ``QCMATRIX_MODULE_DIR``
     - Sets the module directory of QcMatrix library (if Fortran 2003 API built).
     - ``None``

..   * - ``OPENRSP_PERTURBATION_FREE``
       - Enable perturbation free.
       - ``ON``
