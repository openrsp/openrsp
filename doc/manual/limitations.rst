.. _chapter_limitations:

Limitations or Known Problems
=============================

* "T matrix contributions" - i.e. contributions from the perturbed
  "half-time-differentiated" overlap matrix - are not yet supported in the
  newest version of the code. These contributions are only nonzero for
  perturbations that both a) affect the basis set and b) have frequencies other
  than zero. The most relevant such kind of perturbation is the magnetic dipole
  perturbations using London atomic orbitals. Properties consisting of only
  other kinds of perturbations - such as geometric displacement of the nuclei
  or electric dipole perturbations - are unaffected by the lack of T matrix
  contributions.

* Currently we use ``QcPertInt`` (defined as ``QInt`` type in
  ``include/RSPPerturbation.h``, and ``src/fortran/RSPPertBasicTypes.F90`` for
  Fortran APIs) to reprenset several perturbation labels (see
  :ref:`chapter_notations_and_conventions`), in which one label is described by
  ``OPENRSP_PERT_LABEL_BIT`` bits (that can be modified during the step
  ``ccmake``, see :ref:`chapter_installation`).

  For the time being, we do not suggest that users change the type of
  ``QcPertInt``, because other integer types are not supported by OpenRSP yet.

* The current implementation of residues is just tested for electric field
  perturbations and single residues.
