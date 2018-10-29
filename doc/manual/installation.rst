.. _chapter_installation:

Installation
============

Before installing OpenRSP, you need to make sure the following programs are
installed on your computer:

#. Git,
#. CMake (:math:`\ge2.8`),
#. C, C++ (if C++ APIs built) and/or Fortran 2003 (if Fortran APIs built) compilers,
#. HDF 5 (:math:`\ge1.8`) if it is enabled in QcMatrix library,
#. BLAS and LAPACK libraries, and
#. `QcMatrix library <https://gitlab.com/bingao/qcmatrix>`_.

For the time being, only CMake can be used to compile OpenRSP. In general,
OpenRSP should be compiled together with the host programs. See for example the
LSDalton program.

You can also compile OpenRSP alone to be familiar with how it works. But no
real calculations will be performed, all the callback functions in the OpenRSP
unit testing only return pre-defined data or read data from file. Let us
assume that you want to compile the library in directory ``build``, you could
invoke the following commands::

    mkdir build
    cd build
    ccmake ..
    make

During the step ``ccmake``, you need to set some parameters appropriately for
the compilation. For instance, if you enable ``OPENRSP_TEST_EXECUTABLE``, some
executables for the test suite will be built and can run after compilation. So
that you are able to check if OpenRSP has been successfully compiled. A
detailed list of the parameters controlling the compilation is given in the
following table:

.. tabularcolumns:: |L|p{0.5\textwidth}|L|
.. list-table:: OpenRSP CMake parameters
   :header-rows: 1

   * - CMake parameters
     - Description
     - Default
   * - ``OPENRSP_BUILD_WEB``
     - Build OpenRSP from WEB files (only useful for developers)
     - ``OFF``
   * - ``OPENRSP_FORTRAN_API``
     - Build Fortran 2003 API
     - ``OFF``
   * - ``OPENRSP_PERT_LABEL_BIT``
     - Number of bits for a perturbation label (used for perturbation free scheme)
     - 10
   * - ``OPENRSP_TEST_EXECUTABLE``
     - Build test suite as excutables (otherwise, as functions in the library)
     - ``ON``
   * - ``OPENRSP_USER_CONTEXT``
     - Enable user context in callback functions
     - ``OFF``
   * - ``OPENRSP_ZERO_BASED``
     - Zero-based numbering
     - ``ON``
   * - ``QCMATRIX_HEADER_DIR``
     - Directory of header files of QcMatrix library
     - ``None``
   * - ``QCMATRIX_LIB``
     - Name of QcMatrix library with absolute path
     - ``None``
   * - ``QCMATRIX_MODULE_DIR``
     - Directory of Fortran modules of QcMatrix library
     - ``None``

..   * - ``OPENRSP_PERTURBATION_FREE``
       - Enable perturbation free.
       - ``ON``
