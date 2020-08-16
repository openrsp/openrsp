.. _subsection_limitations:

Limitations or Known Problems
-----------------------------

* "T matrix contributions" - i.e. contributions from the perturbed
  "half-time-differentiated" overlap matrix - are not yet supported.
  These contributions are only nonzero for
  perturbations that both a) affect the basis set and b) have frequencies other
  than zero. The most relevant such kind of perturbation is the magnetic dipole
  perturbations using London atomic orbitals. Properties consisting of only
  other kinds of perturbations - such as geometric displacement of the nuclei
  or electric dipole perturbations - are unaffected by the lack of T matrix
  contributions.

* Currently we use ``QcPertInt`` (defined as ``QInt`` type in
  ``include/RSPPerturbation.h``, and ``src/fortran/RSPPertBasicTypes.F90`` for
  Fortran APIs) to reprenset several perturbation labels (see
  :ref:`subsection_notations_and_conventions`), in which one label is described
  by ``OPENRSP_PERT_LABEL_BIT`` bits (that can be modified during the step
  ``ccmake``, see :ref:`subsection_compile`).

  For the time being, we do not suggest that users change the type of
  ``QcPertInt``, because other integer types are not supported by OpenRSP yet.

* The current implementation for calculation of residues of response functions
  is significantly limited in generality. Currently, only electric dipole perturbations
  and single residues are possible; furthermore, there are significant limitations for
  the calculation setup. These limitations are described in further detail in the manual
  of LSDalton in its (at the time of writing unreleased) 2020 version.
