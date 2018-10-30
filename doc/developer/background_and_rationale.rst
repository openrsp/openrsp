.. _chapter_background_and_rationale:

Background and Rationale behind OpenRSP
=======================================

In this chapter, we will discuss the idea behind the implementation, that will
be useful for further development and maintenance, and be useful for new
developers to understand the library and to start their work on top of the
current development.

Theoretical Background
----------------------

The density matrix-based quasienergy formulation of the Kohn-Sham density
functional response theory using perturbation- and time-dependent basis sets
can be found in [Thorvaldsen2008]_ and

.. [Bast2011] Radovan Bast, Ulf Ekström, Bin Gao, Trygve Helgaker, Kenneth Ruud
   and Andreas J. Thorvaldsena, Phys. Chem. Chem. Phys. 13, 2627-2651 (2011).

A relativistic implementation can be found in:

.. [Bast2009] Radovan Bast, Andreas J. Thorvaldsen, Magnus Ringholm and Kenneth Ruud,
   Chem. Phys. 356(1-3), 177-186 (2009).

The recursive programming techniques implemented in OpenRSP can be found in
[Ringholm2014]_.

The recursive programming techniques used for the first order residues can be
found in [Friese2015]_.

Rationale behind OpenRSP
------------------------

The name OpenRSP stands for "open-ended response theory", that is, the library
is:

#. open-ended for different levels of theory, i.e., one-, two- and
   four-component levels;
#. open-ended for different wave functions, e.g., atomic-orbital (AO) based
   density matrix, molecular orbital (MO) cofficients and coupled cluster (CC);
#. open-ended for different kinds of perturbations; and
#. open-ended for different host programs.

For the time being, OpenRSP has implemented:

#. AO based density matrix response theory (source codes in ``src/ao_dens``),

which works for one-, two- and four-component levels by simply setting the
appropriate Hamiltonian. We are now planning to implement the MO and CC based
response theories.

**NOTE**: The codes in ``src/ao_dens`` are written in Fortran, but OpenRSP APIs
are implemented using C language. Therefore, adapter codes between them are
implemented in ``src/ao_dens/adapter``, for OpenRSP APIs calling the codes of
AO based density matrix response theory, also for the AO based density matrix
response theory codes calling the callback functions (as function pointers
saved by OpenRSP APIs).

To make OpenRSP work for any perturbation, we are now trying to implement the
so called **perturbation free scheme**, see :ref:`section_perturbation_free`.

In order to make it easy for implementing OpenRSP into different host programs
(written in different programming languages), we agree to use the **callback
function scheme** in OpenRSP in the 2015 Skibotn meeting.  The callback
functions are specified by host programs by calling the OpenRSP APIs (both C
and Fortran APIs implemented) during run time, and will be used by OpenRSP
during calculations, to get contributions from electronic and nuclear
Hamiltonian, and to get response parameters from solving the linear response
equation.

Another important issue affects the implementation of OpenRSP into different
host programs is the matrix and its different operations that OpenRSP
extensively depends on. Different host programs can have different types of
matrices (dense and sparse, sequential and parallel) and written by different
programming languages (e.g. C and Fortran).

To best utilize the host program's developed matrix routines (if there is), and
also to remove this complexity of matrix problem from OpenRSP, we also agree to
build OpenRSP on top of the `QcMatrix library <https://gitlab.com/bingao/qcmatrix>`_
in the 2015 Skibotn meeting.

The QcMatrix library works as an adapter between OpenRSP and different matrix
routines (implemented in different host programs) that can be written in C and
Fortran. If there is no matrix routines implemented in a host program, it can
fully use the QcMatrix library that will invoke BLAS and LAPACK libraries for
matrix operations.

Therefore, a full picture of OpenRSP used in a C host program can be
(the use of OpenRSP in a Fortran host program can be found in Secion
:ref:`section_fortran_api_impl`):

.. image:: /_static/openrsp_framework.png
   :scale: 30 %
   :align: center

As shown in the above picture, the OpenRSP library is divided into three parts:

#. The "OpenRSP C APIs" have been described in :ref:`chapter_api_reference`
   which work mostly between the host program driver routine and other parts of
   the OpenRSP library;
#. The "OpenRSP response" is the core part in which the response theory
   calculations will be performed;
#. The "OpenRSP C struct" will be described in :ref:`chapter_openrsp_design`,
   which saves the information of perturbations, Hamiltonian and linear
   response equation solver, and will be used by the "OpenRSP response" part
   during calculating response functions and residues.
