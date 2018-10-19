=================
Complete Examples
=================
.. include:: background.rst

Typical procedure of using OpenRSP
==================================
.. include:: background.rst

After reading previous parts, you have already known that you first need to
declare the OpenRSP context and error handler (C code)::

  #inclde "openrsp.h"
  OpenRSP open_rsp;
  QErrorCode ierr;

or (Fortran code)::

  use openrsp_f
  type(OpenRSP) open_rsp
  integer(kind=4) ierr

Afterwards, you could creat the OpenRSP context::

  ierr = OpenRSPCreate(&open_rsp);
  if (ierr!=QSUCCESS) {
      /* error handling */
  }

We will only show the C code here, because the difference between C and
Fortran is not much. Users can refer to Section
"**3.2 Functions of OpenRSP API (Fortran version)**" of the OpenRSP Manual
for the use of OpenRSP Fortran APIs.

After creating the OpenRSP context, users could set: 

#. Perturbations involved in calculations by calling
   :c:func:`OpenRSPSetPerturbations`;
#. Electronic Hamiltonian, by calling [#]_

   #. :c:func:`OpenRSPSetPDBS`,
   #. :c:func:`OpenRSPAddOneOper`,
   #. :c:func:`OpenRSPAddTwoOper`,
   #. :c:func:`OpenRSPAddXCFun`;

.. [#] Users may not need all these 4 APIs. For instance, Hartree-Fock
       calculations do not need to call :c:func:`OpenRSPAddXCFun`.

3. nuclear Hamiltonian by calling :c:func:`OpenRSPSetNucContributions`;
4. Linear response equation solver by calling
   :c:func:`OpenRSPSetLinearRSPSolver`.

After setting the above information, users **must** call :c:func:`OpenRSPAssemble`
to examine if the context of OpenRSP has been set correctly. Otherwise,
calculations could have problems during running.

After calling :c:func:`OpenRSPAssemble`, users could consider using
:c:func:`OpenRSPWrite` to write the OpenRSP context (in a readable
format) into a file. If the file exists, the OpenRSP will append its
context to the file.

This file can be read and sent to the OpenRSP authors if there is anything
wrong during calculations.

If everything is OK, users can then:

#. call :c:func:`OpenRSPGetRSPFun` to calculate response functions, and/or
#. call :c:func:`OpenRSPGetResidue` to calculate residues.

After all calculations performed, users should call :c:func:`OpenRSPDestroy`
to release the memory used by the OpenRSP context.

The above is a typical procedure of using OpenRSP library. Users could
also refer to the tests codes in ``tests/c`` (C version) and ``tests/f90``
(Fortran version), to learn how the callback functions can be prepared.

