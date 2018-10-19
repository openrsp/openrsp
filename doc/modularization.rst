

======================
OpenRSP modularization
======================

This document will contain some topics pertaining to the modularization of
OpenRSP. Currently just a stub to check out the documentation framework.


How to specify the perturbations in a perturbation tuple
--------------------------------------------------------

The perturbation tuple datatype of the OpenRSP library (henceforth openrsp-lib)
is somewhat differently organized from how the perturbation tuples should be
passed between openrsp-lib and the host program. We pass as simple datatypes as
possible to avoid having to handle complicated datatypes.

A perturbation tuple is completely specified by the following::

  integer :: num_pert
  integer, dimension(num_pert) :: pert_dims
  integer, dimension(num_pert) :: pert_first_comp
  integer, dimension(num_pert) :: pert_labels
  complex(8), dimension(num_pert) :: pert_freqs

The integer ``num_pert`` describes the number of perturbations in the
perturbation tuple. The array ``pert_dims`` describes how many compontents each
perturbation has. 

The array ``pert_first_comp`` describes which component of the full-component
perturbation each perturbation starts at. The full-component perturbation is in
this context taken to be the one contained in the host program. The
corresponding number in ``pert_dims`` is therefore the number of components
requested from openrsp-lib. For example, if a perturbation had 12 components in
the host program and one wanted a calculation for components 4 to 9, the
corresponding element of ``pert_first_comp`` would be 4, and the corresponding
element of ``pert_dims`` would be 6.

The array ``pert_labels`` describes which kind of perturbation each
perturbation is. This will simply be passed as integers, since openrsp-lib
doesn't need to know to which kind of perturbation each perturbation is - it
only needs to know which of them are possibly equivalent. For example, the host
program could manage that the perturbations labelled as perturbation no. 1 and
2 are a geometrical perturbation and an electric dipole perturbation,
respectively. Then, ``pert_labels`` for e.g. the (electric) dipole moment
Hessian would be ``(1,1,2)``.

The array ``pert_freq`` contains the values of the frequencies associated with
each perturbation.
