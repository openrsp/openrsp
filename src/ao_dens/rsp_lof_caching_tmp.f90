! Copyright 2012      Magnus Ringholm
!
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

! Contains routines/functions and definitions related to fock_lowerorder datatype
! in which (perturbed) overlap, density and Fock matrices are stored.

module rsp_lof_caching

  use rsp_field_tuple, only: p_tuple,                &
                             p_tuple_standardorder,  &
                             p_tuples_compare,       &
                             p_tuples_standardorder, &
                             merge_p_tuple,          &
                             p1_cloneto_p2
  use rsp_indices_and_addressing
!   use matrix_defop, matrix => openrsp_matrix
!   use matrix_lowlevel, only: mat_init
  use qcmatrix_f

  implicit none

  contains


end module
