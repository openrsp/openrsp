! Copyright 2012      Magnus Ringholm
!
! This source code form is subject to the terms of the
! GNU Lesser General Public License, version 2.1.
! If a copy of the GNU LGPL v2.1 was not distributed with this
! code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

! Contains routines/functions and definitions related to the SDF datatype
! in which (perturbed) overlap, density and Fock matrices are stored.
!
! UPDATE (2016): Now defunct after introduction of new, general 
! 'contribution cache' system

module rsp_sdf_caching

  use rsp_field_tuple
  use rsp_indices_and_addressing
!   use matrix_defop, matrix => openrsp_matrix
!   use matrix_lowlevel, only: mat_init, mat_zero_like
  use qcmatrix_f

  implicit none

  contains

end module
