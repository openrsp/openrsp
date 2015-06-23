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

.. [#] Magnus Ringholm, Dan Jonsson, and Kenneth Ruud,
   J. Comput. Chem., 35, 622 (2014).
.. [#] Andreas J. Thorvaldsen, Kenneth Ruud, Kasper Kristensen,
   Poul Jørgensen and Sonia Coriani, J. Chem. Phys., 129, 214108 (2008).

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
#. reference state (usually the ground state),
   and excited states if calculating residues,
#. perturbations.

OpenRSP Notations and Conventions
=================================
.. include:: background.rst

After the necessary ingredients properly provided, one can invoke OpenRSP
APIs :c:func:`OpenRSPGetRSPFun` or :c:func:`OpenRSPGetResidue` to calculate
the response functions or residues. In the following, we will first describe
how to install OpenRSP and then how to prepare the above ingredients respectively.

Before proceeding, to make yourself familiar with OpenRSP, please refer to
Section "1.2 OpenRSP Notations and Conventions" in OpenRSP Manual for the
notations and conventions used through the OpenRSP library and this tutorial.
