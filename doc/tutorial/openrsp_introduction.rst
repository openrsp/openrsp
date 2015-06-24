============
Introduction
============
.. include:: background.rst

What is OpenRSP?
================
.. include:: background.rst

OpenRSP is a computer library that uses recursive routines [#]_
to identify and assemble contributions to molecular properties
("response functions" or "residues") based on the density
matrix-based response theory [#]_.

Therefore, OpenRSP extensively bases on the matrix operations,
which are built on top of the `QcMatrix library <https://gitlab.com/bingao/qcmatrix>`_.
Please refer to the manual and tutorial of this library if you
are not familiar with it.

.. [#] Magnus Ringholm, Dan Jonsson, and Kenneth Ruud,
   J. Comput. Chem., 35, 622 (2014).
.. [#] Andreas J. Thorvaldsen, Kenneth Ruud, Kasper Kristensen,
   Poul Jørgensen and Sonia Coriani, J. Chem. Phys., 129, 214108 (2008).

.. _slide-ingredients:

Ingredients needed by OpenRSP
=============================
.. include:: background.rst

Briefly, to use OpenRSP, one has to provide:

#. functions for evaluating overlap integrals,
#. different one-electron operators,
#. different two-electron operators,
#. different exchange-correlation functionals,
#. nuclear Hamiltonian,
#. linear response equation solver,
#. perturbations.

The reference state (usually the ground state), and excited states
(if calculating residues) are sent to OpenRSP APIs as input arguments.

.. nextslide::
   :increment:
.. include:: background.rst

Within the above ingredients, one does not need to provide both
one-electron operators, two-electron operators and exchange-correlation
functionals.

Indeed only one of them is needed for OpenRSP to calculate the electron
contributions to the response functions or residues. For instance, only
one- and two-electron operators are needed for the Hartree-Fock calculations.

One therefore only needs to provide OpenRSP the necessary operators
for their interested Hamiltonian.

Callback function scheme
========================
.. include:: background.rst

For the first 6 ingredients, one has to provide OpenRSP appropriate callback
functions that will be invoked by OpenRSP during calculations, whenever the
following quantities needed

#. (perturbed) overlap integrals,
#. (perturbed) one-electron operators,
#. (perturbed) two-electron operators,
#. (perturbed) exchange-correlation functionals,
#. (perturbed) nuclear Hamiltonian,
#. response parameters,

in which the response parameters will be solved from the callback function of
the linear response equation solver.

The use of callback functions make one freely choose the appropriate functions
during runtime.

OpenRSP Notations and Conventions
=================================
.. include:: background.rst

After the necessary ingredients properly provided, one can invoke OpenRSP
APIs :c:func:`OpenRSPGetRSPFun` or :c:func:`OpenRSPGetResidue` to calculate
the response functions or residues.

In the following, we will first describe how to install OpenRSP and then how to
prepare the above ingredients respectively, and perform the calculations.

Before proceeding, to make yourself familiar with OpenRSP, please refer to
Section "**1.2 OpenRSP Notations and Conventions**" in the OpenRSP Manual
for the notations and conventions used through the OpenRSP library and
this tutorial.

