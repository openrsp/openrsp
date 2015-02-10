.. _chapter-openrsp-limitations:

Limitations of OpenRSP
======================

#. Interface of callback functions for (1) linear response equation solver,
   (2) exchange-correlation, and (3) nuclear contributions.

#. Documentation and unit testings.

#. Implemented in DALTON and LSDALTON, integration testings.

#. OpenRSP core part (``ao_dens``) will access C pointer so that host
   program developers do not need to add all operators in a Hamiltonian,
   instead they can choose one-, two-electron operators or exchange-correlation
   functionals.

#. OpenRSPSetRSPEigenSolver(): will implement or change in the future.

#. OpenRSPGetResidue():  will implement or change in the future.

#. OpenRSPSetElecEOM(): will implement or change in the future.

#. OpenRSPSetPerturbations(): will implement or change in the future.
