! Copyright 2012 Magnus Ringholm
!
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

! Contains routines and functions for the calculation of (perturbed) matrices,
! arising in the response equations, that use (perturbed) S, D, and F
! as building blocks. These matrices are W, Y, Z, Lambda, and Zeta.

module rsp_perturbed_matrices

  use matrix_defop
  use rsp_field_tuple
  use rsp_sdf_caching
!   use rsp_general, only: rsp_cfg

  implicit none

  public derivative_superstructure_getsize
  public derivative_superstructure
  public rsp_get_matrix_w
  public rsp_get_matrix_y
  public rsp_get_matrix_z
  public rsp_get_matrix_lambda
  public rsp_get_matrix_zeta
  public get_fds_data_index
  public frequency_zero_or_sum


  type rsp_cfg

     type(matrix) :: zeromat

  end type
  contains

  recursive function derivative_superstructure_getsize(mol, pert, kn, &
                     primed, current_derivative_term) result(superstructure_size)

    implicit none

    logical :: primed
    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    type(p_tuple), dimension(3) :: current_derivative_term
    integer, dimension(2) :: kn
    integer :: i, superstructure_size

    superstructure_size = 0

    if (pert%n_perturbations > 0) then

       superstructure_size = superstructure_size + derivative_superstructure_getsize( &
                             mol,p_tuple_remove_first(pert), kn, primed, &
                             (/p_tuple_extend(current_derivative_term(1), &
                             p_tuple_getone(pert, 1)), current_derivative_term(2:3)/))

       superstructure_size = superstructure_size + derivative_superstructure_getsize( &
                             mol, p_tuple_remove_first(pert), kn, primed, &
                             (/current_derivative_term(1), &
                             p_tuple_extend(current_derivative_term(2), &
                             p_tuple_getone(pert,1)), current_derivative_term(3)/))

       superstructure_size = superstructure_size + derivative_superstructure_getsize( &
                             mol, p_tuple_remove_first(pert), kn, primed, &
                             (/current_derivative_term(1:2), &
                             p_tuple_extend(current_derivative_term(3), &
                             p_tuple_getone(pert, 1))/))

    else

       if (primed .EQV. .TRUE.) then

          if ( ( ( (current_derivative_term(1)%n_perturbations <= kn(2)) .AND.&
              current_derivative_term(2)%n_perturbations <= kn(2) ) .AND. &
              current_derivative_term(3)%n_perturbations <= kn(2) ) .eqv. .TRUE.) then

             superstructure_size = 1

          else

             superstructure_size = 0

          end if

       else

          if ( ( (kn_skip(current_derivative_term(1)%n_perturbations, current_derivative_term(1)%pid, kn) .OR. &
                  kn_skip(current_derivative_term(2)%n_perturbations, current_derivative_term(2)%pid, kn) ) .OR. &
                  kn_skip(current_derivative_term(3)%n_perturbations, current_derivative_term(3)%pid, kn) ) .eqv. &
                  .FALSE.) then

             superstructure_size = 1

          else

             superstructure_size = 0

          end if

       end if

    end if

  end function


  recursive subroutine derivative_superstructure(mol, pert, kn, primed, &
                       current_derivative_term, superstructure_size, & 
                       new_element_position, derivative_structure)

    implicit none

    logical :: primed
    integer :: i, superstructure_size, new_element_position
    integer, dimension(2) :: kn    
    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    type(p_tuple), dimension(3) :: current_derivative_term
    type(p_tuple), dimension(superstructure_size, 3) :: derivative_structure

    if (pert%n_perturbations > 0) then

       call derivative_superstructure(mol, p_tuple_remove_first(pert), kn, primed, &
            (/p_tuple_extend(current_derivative_term(1), p_tuple_getone(pert, 1)), &
            current_derivative_term(2:3)/), superstructure_size, new_element_position, &
            derivative_structure)

       call derivative_superstructure(mol, p_tuple_remove_first(pert), kn, primed, &
            (/current_derivative_term(1), p_tuple_extend(current_derivative_term(2), &
            p_tuple_getone(pert, 1)), current_derivative_term(3)/), &
            superstructure_size, new_element_position, derivative_structure)

       call derivative_superstructure(mol, p_tuple_remove_first(pert), kn, primed, &
            (/current_derivative_term(1:2), p_tuple_extend(current_derivative_term(3), &
            p_tuple_getone(pert, 1))/), superstructure_size, new_element_position, &
            derivative_structure)

    else


       if (primed .EQV. .TRUE.) then

          if ( ( ( (current_derivative_term(1)%n_perturbations <= kn(2)) .AND.&
              current_derivative_term(2)%n_perturbations <= kn(2) ) .AND. &
              current_derivative_term(3)%n_perturbations <= kn(2) ) .eqv. .TRUE.) then

             new_element_position = new_element_position + 1
             derivative_structure(new_element_position, :) = current_derivative_term
 
          end if

       else

          if ( ( (kn_skip(current_derivative_term(1)%n_perturbations,  &
                          current_derivative_term(1)%pid, kn) .OR. &
                  kn_skip(current_derivative_term(2)%n_perturbations,  &
                          current_derivative_term(2)%pid, kn) ) .OR. &
                  kn_skip(current_derivative_term(3)%n_perturbations, &
                          current_derivative_term(3)%pid, kn) ) .eqv. .FALSE.) then

             new_element_position = new_element_position + 1 
             derivative_structure(new_element_position, :) = current_derivative_term(:)

          end if

       end if

    end if

  end subroutine


  subroutine rsp_get_matrix_w(mol, superstructure_size, &
           deriv_struct, total_num_perturbations, which_index_is_pid, &
           indices_len, ind, F, D, S, W)

    implicit none

    integer :: i, total_num_perturbations, superstructure_size, indices_len
    type(rsp_cfg) :: mol
    type(p_tuple), dimension(superstructure_size, 3) :: deriv_struct
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: ind
    type(sdf) :: F, D, S
    type(matrix) :: W, A, B, C

    W%elms_alpha = 0.0
!     W = mat_zero_like(mol%zeromat)

    A = mat_alloc_like(mol%zeromat)
    A%elms_alpha = 0.0
!     A = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(A)

    B = mat_alloc_like(mol%zeromat)
    B%elms_alpha = 0.0
!     B = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(B)

    C = mat_alloc_like(mol%zeromat)
    C%elms_alpha = 0.0
!     C = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(C)

    do i = 1, superstructure_size

!     A = mat_zero_like(mol%zeromat)
!     B = mat_zero_like(mol%zeromat)
!     C = mat_zero_like(mol%zeromat)

       call sdf_getdata_s(D, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(F, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(D, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       W = W + A * B * C


if (.not.(frequency_zero_or_sum(deriv_struct(i,1)) == 0.0) .and. &
    .not.(frequency_zero_or_sum(deriv_struct(i,3)) == 0.0)) then

       call sdf_getdata_s(D, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(S, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(D, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       W = W + ((1.0)/(2.0)) * (frequency_zero_or_sum(deriv_struct(i,3)) - &
                                frequency_zero_or_sum(deriv_struct(i,1))) * &
               A * B * C

elseif (.not.(frequency_zero_or_sum(deriv_struct(i,1)) == 0.0) .and. &
             (frequency_zero_or_sum(deriv_struct(i,3)) == 0.0)) then

       call sdf_getdata_s(D, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(S, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(D, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       W = W + ((-1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,1)) *  &
               A * B * C

elseif (.not.(frequency_zero_or_sum(deriv_struct(i,3)) == 0.0) .and. &
             (frequency_zero_or_sum(deriv_struct(i,1)) == 0.0)) then


       call sdf_getdata_s(D, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(S, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(D, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       W = W + ((1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,3))  * &
               A * B * C

end if


    end do

    A = 0
    B = 0
    C = 0

  end subroutine


  subroutine rsp_get_matrix_y(mol, superstructure_size, deriv_struct, &
           total_num_perturbations, which_index_is_pid, indices_len, &
           ind, F, D, S, Y)

    implicit none

    integer :: i, total_num_perturbations, superstructure_size, indices_len
    type(rsp_cfg) :: mol
    type(p_tuple), dimension(superstructure_size, 3) :: deriv_struct
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: ind
    type(sdf) :: F, D, S
    type(matrix) :: Y, A, B, C

    Y%elms_alpha = 0.0
!     Y = mat_zero_like(mol%zeromat)

    A = mat_alloc_like(mol%zeromat)
    A%elms_alpha = 0.0
!     A = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(A)

    B = mat_alloc_like(mol%zeromat)
    B%elms_alpha = 0.0
!     B = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(B)

    C = mat_alloc_like(mol%zeromat)
    C%elms_alpha = 0.0
!     C = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(C)

! write(*,*) 'ind is', ind
    
    do i = 1, superstructure_size

! write(*,*) 'i is', i

       call sdf_getdata_s(F, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(S, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)



       Y = Y + A*B*C



       call sdf_getdata_s(S, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(F, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)


       Y = Y - A*B*C



if (.not.(frequency_zero_or_sum(deriv_struct(i,1)) == 0.0) .and. &
    .not.(frequency_zero_or_sum(deriv_struct(i,3)) == 0.0)) then

! MR: MAKE SURE THAT THESE (AND B) ARE ACTUALLY THE CORRECT MATRICES TO USE HERE AND BELOW

       call sdf_getdata_s(S, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(S, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Y = Y + ((-1.0)/(2.0)) * (frequency_zero_or_sum(deriv_struct(i,3)) + &
                                frequency_zero_or_sum(deriv_struct(i,1))) * &
               A * B * C

elseif (.not.(frequency_zero_or_sum(deriv_struct(i,1)) == 0.0) .and. &
             (frequency_zero_or_sum(deriv_struct(i,3)) == 0.0)) then

       call sdf_getdata_s(S, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(S, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Y = Y + ((-1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,1)) *  &
               A * B * C

elseif (.not.(frequency_zero_or_sum(deriv_struct(i,3)) == 0.0) .and. &
             (frequency_zero_or_sum(deriv_struct(i,1)) == 0.0)) then

       call sdf_getdata_s(S, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(S, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Y = Y + ((-1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,3))  * &
               A * B * C

end if




    end do

    A = 0
    B = 0
    C = 0

  end subroutine


  subroutine rsp_get_matrix_z(mol, superstructure_size, deriv_struct, kn, &
           total_num_perturbations, which_index_is_pid, indices_len, &
           ind, F, D, S, Z)

    implicit none

    integer :: i, total_num_perturbations, superstructure_size, indices_len
    type(rsp_cfg) :: mol
    type(p_tuple), dimension(superstructure_size, 3) :: deriv_struct
    type(p_tuple) :: merged_p_tuple
    integer, dimension(2) :: kn
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: ind
    type(sdf) :: F, D, S
    type(matrix) :: Z, A, B, C

    Z%elms_alpha = 0.0
!     Z = mat_zero_like(mol%zeromat)

    A = mat_alloc_like(mol%zeromat)
    A%elms_alpha = 0.0
!     A = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(A)

    B = mat_alloc_like(mol%zeromat)
    B%elms_alpha = 0.0
!     B = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(B)

    C = mat_alloc_like(mol%zeromat)
    C%elms_alpha = 0.0
!     C = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(C)


    do i = 1, superstructure_size

       call sdf_getdata_s(D, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(S, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(D, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)


       Z = Z + A*B*C

    end do

    merged_p_tuple = merge_p_tuple(deriv_struct(1,1), merge_p_tuple(deriv_struct(1,2), deriv_struct(1,3)))

    if (kn_skip(total_num_perturbations, merged_p_tuple%pid, kn) .eqv. .FALSE.) then

    A = mat_zero_like(mol%zeromat)

        call sdf_getdata_s(D, merged_p_tuple, get_fds_data_index(merged_p_tuple, &
        total_num_perturbations, which_index_is_pid, indices_len, ind), A)

        Z = Z - A

    end if

call p_tuple_deallocate(merged_p_tuple)

    A = 0
    B = 0
    C = 0

!  write(*,*) 'end of function, returning', Z%elms_alpha

  end subroutine


  subroutine rsp_get_matrix_lambda(mol, p_tuple_a, superstructure_size, deriv_struct, &
           total_num_perturbations, which_index_is_pid, indices_len, ind, D, S, L)

    implicit none

    integer :: i, total_num_perturbations, superstructure_size, indices_len
    type(rsp_cfg) :: mol
    type(p_tuple) :: p_tuple_a, merged_A, merged_B
    type(p_tuple), dimension(superstructure_size, 3) :: deriv_struct
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: ind
    type(sdf) :: D, S
    type(matrix) :: L, A, B, C

    L%elms_alpha = 0.0
!     L = mat_zero_like(mol%zeromat)

    A = mat_alloc_like(mol%zeromat)
    A%elms_alpha = 0.0
!     A = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(A)

    B = mat_alloc_like(mol%zeromat)
    B%elms_alpha = 0.0
!     B = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(B)

    C = mat_alloc_like(mol%zeromat)
    C%elms_alpha = 0.0
!     C = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(C)


    do i = 1, superstructure_size

       merged_A = merge_p_tuple(p_tuple_a, deriv_struct(i,1))
       merged_B = merge_p_tuple(p_tuple_a, deriv_struct(i,3))

       call sdf_getdata_s(D, merged_A, get_fds_data_index(merged_A, &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(S, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(D, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       L = L + A * B * C
        
       call sdf_getdata_s(D, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, merged_B, get_fds_data_index(merged_B, &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       L = L - A * B * C

call p_tuple_deallocate(merged_A)
call p_tuple_deallocate(merged_B)
       
    end do

    A = 0
    B = 0
    C = 0

  end subroutine


!MR: ZETA ROUTINE STILL NOT OPTIMIZED


  subroutine rsp_get_matrix_zeta(mol, p_tuple_a, kn, superstructure_size, deriv_struct, &
           total_num_perturbations, which_index_is_pid, indices_len, &
           ind, F, D, S, Zeta)

    implicit none

    integer :: i, total_num_perturbations, superstructure_size, indices_len
    type(rsp_cfg) :: mol
    type(p_tuple) :: p_tuple_a, merged_p_tuple, merged_A, merged_B
    type(p_tuple), dimension(superstructure_size, 3) :: deriv_struct
    integer, dimension(2) :: kn
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: ind
    type(sdf) :: F, D, S
    type(matrix) :: Zeta, A, B, C

    Zeta%elms_alpha = 0.0
!     Zeta = mat_zero_like(mol%zeromat)

    A = mat_alloc_like(mol%zeromat)
    A%elms_alpha = 0.0
!     A = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(A)

    B = mat_alloc_like(mol%zeromat)
    B%elms_alpha = 0.0
!     B = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(B)

    C = mat_alloc_like(mol%zeromat)
    C%elms_alpha = 0.0
!     C = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(C)

    do i = 1, superstructure_size

       merged_A = merge_p_tuple(p_tuple_a, deriv_struct(i,1))
       merged_B = merge_p_tuple(p_tuple_a, deriv_struct(i,3))

       call sdf_getdata_s(F, merged_A, get_fds_data_index(merged_A, &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(S, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Zeta = Zeta + A * B * C

       call sdf_getdata_s(F, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(S, merged_B, get_fds_data_index(merged_B, &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Zeta = Zeta - A * B * C


if (.not.(frequency_zero_or_sum(deriv_struct(i,1)) == 0.0) .and. &
    .not.(frequency_zero_or_sum(deriv_struct(i,2)) == 0.0)) then

       call sdf_getdata_s(S, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)

       Zeta = Zeta + ( ((1.0)/(2.0))*frequency_zero_or_sum(deriv_struct(i,1)) + &
                       frequency_zero_or_sum(deriv_struct(i,2)) ) * A * B * C


elseif (.not.(frequency_zero_or_sum(deriv_struct(i,1)) == 0.0) .and. &
             (frequency_zero_or_sum(deriv_struct(i,2)) == 0.0)) then

       call sdf_getdata_s(S, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)

       Zeta = Zeta + ((1.0)/(2.0))*frequency_zero_or_sum(deriv_struct(i,1)) * A * B * C

elseif (.not.(frequency_zero_or_sum(deriv_struct(i,2)) == 0.0) .and. &
             (frequency_zero_or_sum(deriv_struct(i,1)) == 0.0)) then

       call sdf_getdata_s(S, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)

       Zeta = Zeta + frequency_zero_or_sum(deriv_struct(i,2)) * A * B * C


end if


       call sdf_getdata_s(S, deriv_struct(i,1), get_fds_data_index(deriv_struct(i,1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(D, deriv_struct(i,2), get_fds_data_index(deriv_struct(i,2), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), B)
       call sdf_getdata_s(F, merged_B, get_fds_data_index(merged_B, &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Zeta = Zeta +  A * B * C



       call sdf_getdata_s(S, merged_A, get_fds_data_index(merged_A, &
            total_num_perturbations, which_index_is_pid, indices_len, ind), A)
       call sdf_getdata_s(F, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Zeta = Zeta - A * B * C





if (.not.(frequency_zero_or_sum(deriv_struct(i,2)) == 0.0) .and. &
    .not.(frequency_zero_or_sum(deriv_struct(i,3)) == 0.0)) then

       call sdf_getdata_s(S, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Zeta = Zeta + ( ((1.0)/(2.0))*frequency_zero_or_sum(deriv_struct(i,3)) + &
                     frequency_zero_or_sum(deriv_struct(i,2)) ) * A * B * C


elseif (.not.(frequency_zero_or_sum(deriv_struct(i,2)) == 0.0) .and. &
             (frequency_zero_or_sum(deriv_struct(i,3)) == 0.0)) then

       call sdf_getdata_s(S, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Zeta = Zeta + ((1.0)/(2.0))*frequency_zero_or_sum(deriv_struct(i,3)) * A * B * C


elseif (.not.(frequency_zero_or_sum(deriv_struct(i,3)) == 0.0) .and. &
             (frequency_zero_or_sum(deriv_struct(i,2)) == 0.0)) then

       call sdf_getdata_s(S, deriv_struct(i,3), get_fds_data_index(deriv_struct(i,3), &
            total_num_perturbations, which_index_is_pid, indices_len, ind), C)

       Zeta = Zeta + frequency_zero_or_sum(deriv_struct(i,2)) * A * B * C



end if


call p_tuple_deallocate(merged_A)
call p_tuple_deallocate(merged_B)
         
    end do

    merged_p_tuple = merge_p_tuple(p_tuple_a, merge_p_tuple(deriv_struct(1,1), &
                     merge_p_tuple(deriv_struct(1,2), deriv_struct(1,3))))

    if (kn_skip(merged_p_tuple%n_perturbations, &
        merged_p_tuple%pid, kn) .eqv. .FALSE.) then

    A = mat_zero_like(mol%zeromat)

       call sdf_getdata_s(F, merged_p_tuple, get_fds_data_index(merged_p_tuple, & 
       total_num_perturbations, which_index_is_pid, indices_len, ind), A)

       Zeta = Zeta - A

    end if

call p_tuple_deallocate(merged_p_tuple)

    A = 0
    B = 0
    C = 0

  end subroutine


  function get_fds_data_index(pert_tuple, total_num_perturbations, which_index_is_pid, &
                              indices_len, indices)

    implicit none

    type(p_tuple) :: pert_tuple
    integer :: i, total_num_perturbations, indices_len
    integer, allocatable, dimension(:) :: get_fds_data_index
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: indices

    if (pert_tuple%n_perturbations == 0) then

       allocate(get_fds_data_index(1))
       get_fds_data_index(1) = 1

    else

       allocate(get_fds_data_index(pert_tuple%n_perturbations))

       do i = 1, pert_tuple%n_perturbations

          get_fds_data_index(i) = indices(which_index_is_pid(pert_tuple%pid(i)))

       end do

    end if

  end function


  function frequency_zero_or_sum(pert_tuple)

    implicit none

    type(p_tuple) :: pert_tuple
    complex(8) :: frequency_zero_or_sum
    integer :: i

    frequency_zero_or_sum = 0.0

    if (pert_tuple%n_perturbations > 0) then

       do i = 1, pert_tuple%n_perturbations

          frequency_zero_or_sum = frequency_zero_or_sum + pert_tuple%freq(i)

       end do

    end if

  end function

end module