.. _chapter_limitations:

Limitations or Known Problems
=============================

**Are the following two issues also true for OpenRSP?**

* Not safe to use "regular" DFT on Dalton's openrsp branch.
* vibgamma test broken on GNU.

* Currently we use ``QcPertInt`` (defined as ``QInt`` type in
  ``include/RSPPerturbation.h``, and ``src/fortran/RSPPertBasicTypes.F90`` for
  Fortran APIs) to reprenset several perturbation labels (see
  :ref:`chapter_notations_and_conventions`), in which one label is described by
  ``OPENRSP_PERT_LABEL_BIT`` bits (that can be modified during the step
  ``ccmake``, see :ref:`chapter_installation`).

  For the time being, we do not suggest that users change the type of
  ``QcPertInt``, because other integer types are not supported by OpenRSP yet.

**To Daniel: what is the limitation for residue calculations? Only single and double residues supported?**
