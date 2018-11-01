.. _chapter_getting_started:

Getting Started
===============

OpenRSP is a computer library that uses recursive routines [Ringholm2014]_ to
identify and assemble contributions to molecular properties ("response
functions" or "residues") based on the density matrix-based response theory
[Thorvaldsen2008]_.

Therefore, OpenRSP extensively bases on the matrix operations, which are built
on top of the `QcMatrix library <https://gitlab.com/bingao/qcmatrix>`_. Please
refer to the manual and tutorial of this library if you are not familiar with
it.

Briefly, to use OpenRSP, one has to provide:

#. perturbations,
#. functions for evaluating overlap integrals,
#. different one-electron operators,
#. different two-electron operators,
#. different exchange-correlation functionals,
#. different zero-electorn operators, like nuclear Hamiltonian,
#. linear response equation solver.

The reference state (usually the ground state), and excited states (if
calculating residues) are sent to OpenRSP APIs as input arguments.

Within the above ingredients, one does not need to provide both one-electron
operators, two-electron operators and exchange-correlation functionals.

Indeed only one of them is needed for OpenRSP to calculate the electron
contributions to the response functions or residues. For instance, only one-
and two-electron operators are needed for the Hartree-Fock calculations.

One therefore only needs to provide OpenRSP the necessary operators for their
interested Hamiltonian, which is done by providing appropriate callback
functions to OpenRSP:

#. rank of a perturbation component (**not invoked by the current release**),
#. (perturbed) overlap integrals,
#. (perturbed) one-electron operators,
#. (perturbed) two-electron operators,
#. (perturbed) exchange-correlation functionals,
#. (perturbed) zero-eletron operators,
#. response parameters solved from the linear response equation.

The use of callback functions makes one freely choose the appropriate functions
during runtime to calculate molecular properties.

After the necessary ingredients properly provided, one can invoke OpenRSP APIs
:c:func:`OpenRSPGetRSPFun` or :c:func:`OpenRSPGetResidue` to calculate response
functions or residues.

Perturbation free
-----------------

Although it is not available in the current release, we would like to mention
that, the use of callback function to determine the rank of a perturbation
component (see definition of perturbation component and rank in
:ref:`chapter_notations_and_conventions`) can make OpenRSP a **perturbation
free** library in the future. That is:

#. OpenRSP does not need to know the meaning of each perturbation.
   All perturbations are treated equally as symbols/variables in
   OpenRSP for differentiation.
#. Only the order of perturbation labels matters, that OpenRSP will
   follow to generate necessary perturbation tuples (see definition of
   perturbation label and tuple in :ref:`chapter_notations_and_conventions`)
   during calculations.
#. These labels are sent to OpenRSP by the APIs :c:func:`OpenRSPGetRSPFun`
   or :c:func:`OpenRSPGetResidue` during runtime.

Requirements on callback functions
----------------------------------

OpenRSP further has the following requirements on the callback functions that
users should be keep in mind (more detail can be found at the beginning of
:ref:`chapter_callback_functions`):

1. OpenRSP always ask for **complex expectation values** for different zero-,
   one- and two-electron operators, and exchange-correlation functionals, and
   these values are presented in memory that the real and imaginary parts of
   each value are consecutive.

2. OpenRSP requires that calculated integral matrices and expectation values
   should **be added to the returned argument**. OpenRSP will zero the entries
   of these matrices and expectation values at first.

   This requirement affects the callback functions of zero-, one- and
   two-electron operators, and exchange-correlation functionals.

Typical procedure of using OpenRSP
----------------------------------

To summarize, you first need to declare the OpenRSP context and error handler
(C code) for using OpenRSP::

  #inclde "OpenRSP.h"
  OpenRSP open_rsp;
  QErrorCode ierr;

or (Fortran code)::

  use OpenRSP_f
  type(OpenRSP) open_rsp
  integer(kind=4) ierr

Afterwards, you can creat the OpenRSP context (we only show the C code here,
because the difference between C and Fortran is not much)::

  ierr = OpenRSPCreate(&open_rsp, num_atoms);
  if (ierr!=QSUCCESS) {
      /* error handling */
  }

**NOTE**: the last argument ``num_atoms`` in the API :c:func:`OpenRSPCreate` is
the number of atoms, which **will be removed** after the perturbation free
scheme implemented in OpenRSP.

After creating the OpenRSP context, users could set:

#. Perturbations involved in calculations by calling
   :c:func:`OpenRSPSetPerturbations`;

#. Electronic Hamiltonian, by calling

   #. :c:func:`OpenRSPSetOverlap`,
   #. :c:func:`OpenRSPAddOneOper`,
   #. :c:func:`OpenRSPAddTwoOper`,
   #. :c:func:`OpenRSPAddXCFun`;

   Note that users may not need all the above 4 APIs. For instance,
   Hartree-Fock calculations do not need to call :c:func:`OpenRSPAddXCFun`.

#. Zero-electron operator, like nuclear Hamiltonian by calling
   :c:func:`OpenRSPAddZeroOper`;

#. Linear response equation solver by calling
   :c:func:`OpenRSPSetLinearRSPSolver`.

After setting the above information, users **must** call
:c:func:`OpenRSPAssemble` to examine if the context of OpenRSP has been set
correctly. Otherwise, calculations could have problems during running.

Afterwards, users could use :c:func:`OpenRSPWrite` to write the OpenRSP context
(in a readable format) into a file. If the file exists, the OpenRSP will append
its context to the file.

This file can be read and sent to the OpenRSP authors if there is anything
wrong during calculations.

If everything is OK, users can then:

#. call :c:func:`OpenRSPGetRSPFun` to calculate response functions, and/or
#. call :c:func:`OpenRSPGetResidue` to calculate residues.

After all calculations performed, users should call :c:func:`OpenRSPDestroy` to
release the memory used by the OpenRSP context.

The above is a typical procedure of using OpenRSP. Users can also refer to the
unit testing codes in the directory ``tests`` (C version), to learn how the
callback functions can be prepared.

In the following, we will describe how to prepare the above ingredients
respectively, and perform the calculations step by step.

Before proceeding, to make yourself familiar with OpenRSP, please refer to
:ref:`chapter_notations_and_conventions` for the notations and conventions used
through the OpenRSP and this tutorial.
