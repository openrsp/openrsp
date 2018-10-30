.. _chapter_limitations:

Limitations or Known Problems
=============================

* "T matrix contributions" - i.e. contributions from the perturbed "half-time-differentiated" overlap matrix - are not yet supported in the newest version of the code. These contributions are only nonzero for perturbations that both a) affect the basis set and b) have frequencies other than zero. The most relevant such kind of perturbation is the magnetic dipole perturbations using London atomic orbitals. Properties consisting of only other kinds of perturbations - such as geometric displacement of the nuclei or electric dipole perturbations - are unaffected by the lack of T matrix contributions.

**To Daniel: what is the limitation for residue calculations? Only single and double residues supported?**
