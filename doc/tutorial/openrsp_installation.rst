============
Installation
============
.. include:: background.rst

CMake parameters
================
.. include:: background.rst

The OpenRSP Manual has described how to compile OpenRSP using the CMake,
see Chapter 2 "**INSTALLATION**". There are some parameters affecting
the compiling of OpenRSP:

* ``OPENRSP_USER_CONTEXT`` will enable the user defined context in
  callback functions. OpenRSP only talks to callback functions with
  minimal and necessary numbers of arguments. If one needs to pass
  additional but necessary information to their callback functions,
  they can pass it through the user defined context.

  In the ``test/c`` we only show a simple use of the user defined
  context to account the times of callback function being invoked.

.. nextslide::
   :increment:
.. include:: background.rst

* ``LIB_OPENRSP_NAME`` specifies the name of compiled OpenRSP
  library so that the OpenRSP library will be ``lib|LIB_OPENRSP_NAME|.a``,
  where the specified name will substitute ``|LIB_OPENRSP_NAME|``.

* ``OPENRSP_CXX_API`` not implemented yet.

* ``OPENRSP_F03_API`` will build OpenRSP Fortran 2003 APIs, so that one
  can use OpenRSP in their Fortran codes.

* ``OPENRSP_TEST_EXECUTABLE`` will build the test suite excutables that
  can be run directly for OpenRSP unit testing; otherwise, the unit testing
  will still be compiled but put into the OpenRSP library that can be invoked
  by user to do test.

.. nextslide::
   :increment:
.. include:: background.rst

The following CMake parameters provide the information of the QcMatrix
library:

* ``QCMATRIX_LIB`` specifies the QcMatrix library as ``xx/libqcmatrix.a``
  that will be linked if test suite excutables built.

* ``QCMATRIX_INCLUDE_DIR`` gives the include directory of the QcMatrix
  library, which should contains its ``include`` and ``build`` directories.

* ``QCMATRIX_MODULE_DIR`` sets the module directory of the QcMatrix library,
  which is normally its ``build`` directory; this is needed only if the OpenRSP
  Fortran 2003 APIs are built.

OpenRSP testing
===============
.. include:: background.rst

After successfully compiled the OpenRSP library, one is recommended to perform
the OpenRSP unit testing. If the test suite excutables built, one can run them
by ``./|LIB_OPENRSP_NAME|_c_test`` and ``./|LIB_OPENRSP_NAME|_f_test`` if the
OpenRSP Fortran 2003 APIs are built.

Moreover, one should also implement integration testing if they use the OpenRSP
library in their codes.

The authors of OpenRSP can be contacted if there is any error during compiling,
testing and/or using the OpenRSP library.

