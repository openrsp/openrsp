===============
OpenRSP Context
===============
.. include:: background.rst

OpenRSP APIs
============
.. include:: background.rst

In order to use OpenRSP, C users should first include the header file
of OpenRSP in their codes::

  #inclde "openrsp.h"

while Fortran users should use the OpenRSP module::

  use openrsp_f

All the OpenRSP APIs (application programming interface) can be invoked as::

  OpenRSP open_rsp;
  QErrorCode ierr;
  ierr = OpenRSP...(&openrsp, ...);

or for Fortran users as::

  type(OpenRSP) open_rsp
  integer(kind=4) ierr
  ierr = OpenRSP..._f(open_rsp, ...)

where ``open_rsp`` contains the context of calculations by the OpenRSP,
and is always of the first argument for all the APIs.

The ``ierr`` contains error information that one should check if it
equals to ``QSUCCESS`` (constant defined in
`QcMatrix library <https://gitlab.com/bingao/qcmatrix>`_). If not,
there was error happened in the invoked OpenRSP API, and one should
stop the calculations and check the error message.

Basic OpenRSP APIs
==================
.. include:: background.rst

There are several OpenRSP APIs only take the OpenRSP context
``open_rsp`` as the argument and must be invoked by the users
during calculations:

.. c:function:: QErrorCode OpenRSPCreate(OpenRSP *open_rsp)

.. c:function:: QErrorCode OpenRSPAssemble(OpenRSP *open_rsp)

.. c:function:: QErrorCode OpenRSPDestroy(OpenRSP *open_rsp)

in which :c:func:`OpenRSPCreate` and :c:func:`OpenRSPDestroy`
must be called respectively **at the beginning** and **at the end**
of the calculations, to create and destroy the context of the
OpenRSP library. They should be called **only once**.

The API :c:func:`OpenRSPAssemble` should be called **after** all
ingredients (see :ref:`slide-ingredients`) have been set, and
**before** any response function or residue calculation.

This API will examine if the context of OpenRSP has been set
correctly, and must be called **at least once**.

Check the OpenRSP context
=========================
.. include:: background.rst

Often users would like to see how the OpenRSP context has been set
in a readable manner, that can be done by calling:

.. c:function:: QErrorCode OpenRSPWrite(OpenRSP *open_rsp, QChar *file_name)

The OpenRSP context will be written (or more exactly **appended**)
into the end of the file ``|file_name|`` that can be read and sent
to the OpenRSP authors if there is anything wrong.

This API can be called many times as you want.

Fortran users
=============
.. include:: background.rst

Functions of OpenRSP API (Fortran) are similar to those of the C version,
except that an extra ``_f`` should be appended to each function.

Other differences have been described in the OpenRSP Manual (Section
"**3.2 Functions of OpenRSP API (Fortran version)**").

