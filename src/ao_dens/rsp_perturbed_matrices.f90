! Copyright 2012 Magnus Ringholm
!
! This source code form is subject to the terms of the
! GNU Lesser General Public License, version 2.1.
! If a copy of the GNU LGPL v2.1 was not distributed with this
! code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

! Contains routines and functions for the calculation of (perturbed) matrices,
! arising in the response equations, that use (perturbed) S, D, and F
! as building blocks. These matrices are W, Y, Z, Lambda, and Zeta.

module rsp_perturbed_matrices

!  use matrix_defop, matrix => openrsp_matrix
!  use matrix_lowlevel, only: mat_init
  use rsp_field_tuple, only: p_tuple,              &
                             p_tuple_extend,       &
                             p_tuple_getone,       &
                             merge_p_tuple,        &
                             p_tuple_remove_first, &
                             p_tuple_deallocate
  use rsp_sdf_caching
  use rsp_indices_and_addressing
  use rsp_property_caching
  use qcmatrix_f

  implicit none

  public derivative_superstructure_getsize
  public derivative_superstructure
  public get_fds_data_index
  public frequency_zero_or_sum
  
  public rsp_get_matrix_w
  public rsp_get_matrix_y
  public rsp_get_matrix_z
  public rsp_get_matrix_lambda
  public rsp_get_matrix_zeta

  contains

  ! Get size of derivative superstructure for perturbed matrices call
  recursive function derivative_superstructure_getsize(pert, kn, &
                     primed, current_derivative_term) result(superstructure_size)

    implicit none

    logical :: primed
    type(p_tuple) :: pert
    type(p_tuple), dimension(3) :: current_derivative_term, next_deriv_term
    integer, dimension(2) :: kn
    integer :: i, superstructure_size

    superstructure_size = 0

    ! Recurse
    if (pert%npert > 0) then

       next_deriv_term = (/p_tuple_extend(current_derivative_term(1), &
                             p_tuple_getone(pert, 1)), current_derivative_term(2:3)/)

       superstructure_size = superstructure_size + derivative_superstructure_getsize( &
                             p_tuple_remove_first(pert), kn, primed, next_deriv_term)

       next_deriv_term = (/current_derivative_term(1), &
                             p_tuple_extend(current_derivative_term(2), &
                             p_tuple_getone(pert,1)), current_derivative_term(3)/)

       superstructure_size = superstructure_size + derivative_superstructure_getsize( &
                             p_tuple_remove_first(pert), kn, primed, next_deriv_term)

       next_deriv_term = (/current_derivative_term(1:2), &
                             p_tuple_extend(current_derivative_term(3), &
                             p_tuple_getone(pert, 1))/)

       superstructure_size = superstructure_size + derivative_superstructure_getsize( &
                             p_tuple_remove_first(pert), kn, primed, next_deriv_term)

    ! End of recursion: Determine if term is kept and if so, increment counter
    else
   
       if (primed .EQV. .TRUE.) then

          if ( ( ( (current_derivative_term(1)%npert <= kn(2)) .AND.&
              current_derivative_term(2)%npert <= kn(2) ) .AND. &
              current_derivative_term(3)%npert <= kn(2) ) .eqv. .TRUE.) then

             superstructure_size = 1

          else

             superstructure_size = 0

          end if

       else

          if ( ( (kn_skip(current_derivative_term(1)%npert, current_derivative_term(1)%pid, kn) .OR. &
                  kn_skip(current_derivative_term(2)%npert, current_derivative_term(2)%pid, kn) ) .OR. &
                  kn_skip(current_derivative_term(3)%npert, current_derivative_term(3)%pid, kn) ) .eqv. &
                  .FALSE.) then
                  
             superstructure_size = 1

          else
          
             superstructure_size = 0

          end if

       end if

    end if

  end function

  ! Get derivative superstructure for perturbed matrices call
  recursive subroutine derivative_superstructure(pert, kn, primed, &
                       current_derivative_term, superstructure_size, & 
                       new_element_position, derivative_structure)

    implicit none

    logical :: primed
    integer :: i, superstructure_size, new_element_position
    integer, dimension(2) :: kn    
    type(p_tuple) :: pert
    type(p_tuple), dimension(3) :: current_derivative_term, next_deriv_term
    type(p_tuple), dimension(superstructure_size, 3) :: derivative_structure

    ! Recurse
    if (pert%npert > 0) then

       next_deriv_term = (/p_tuple_extend(current_derivative_term(1), &
                           p_tuple_getone(pert, 1)), current_derivative_term(2:3)/)


       call derivative_superstructure(p_tuple_remove_first(pert), kn, primed, &
            next_deriv_term, superstructure_size, new_element_position, &
            derivative_structure)

       next_deriv_term = (/current_derivative_term(1), &
                             p_tuple_extend(current_derivative_term(2), &
                             p_tuple_getone(pert,1)), current_derivative_term(3)/)


       call derivative_superstructure(p_tuple_remove_first(pert), kn, primed, &
            next_deriv_term, &
            superstructure_size, new_element_position, derivative_structure)

       next_deriv_term = (/current_derivative_term(1:2), &
                             p_tuple_extend(current_derivative_term(3), &
                             p_tuple_getone(pert, 1))/)


       call derivative_superstructure(p_tuple_remove_first(pert), kn, primed, &
            next_deriv_term, superstructure_size, new_element_position, &
            derivative_structure)

    ! End of recursion: Determine if term is kept and if so, store the superstructure element
    else

       next_deriv_term = current_derivative_term

       if (primed .EQV. .TRUE.) then

          if ( ( ( (current_derivative_term(1)%npert <= kn(2)) .AND.&
              current_derivative_term(2)%npert <= kn(2) ) .AND. &
              current_derivative_term(3)%npert <= kn(2) ) .eqv. .TRUE.) then

             new_element_position = new_element_position + 1
             derivative_structure(new_element_position, :) = next_deriv_term
 
          end if

       else

          if ( ( (kn_skip(current_derivative_term(1)%npert,  &
                          current_derivative_term(1)%pid, kn) .OR. &
                  kn_skip(current_derivative_term(2)%npert,  &
                          current_derivative_term(2)%pid, kn) ) .OR. &
                  kn_skip(current_derivative_term(3)%npert, &
                          current_derivative_term(3)%pid, kn) ) .eqv. .FALSE.) then

             new_element_position = new_element_position + 1 
             derivative_structure(new_element_position, :) = next_deriv_term(:)

          end if

       end if

    end if

  end subroutine


  ! Calculate a perturbed W matrix
  subroutine rsp_get_matrix_w(superstructure_size, &
           deriv_struct, total_num_perturbations, which_index_is_pid, &
           indices_len, ind, len_f, F, len_d, D, len_s, S, W)

    implicit none

    integer :: i, total_num_perturbations, superstructure_size, indices_len
    type(p_tuple), dimension(superstructure_size, 3) :: deriv_struct
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: ind
    type(p_tuple), dimension(1) :: ds1, ds2, ds3
    integer, dimension(:), allocatable :: iu1, iu2, iu3
    integer :: len_f, len_d, len_s
    type(contrib_cache_outer), dimension(len_f) :: F
    type(contrib_cache_outer), dimension(len_d) :: D
    type(contrib_cache_outer), dimension(len_S) :: S
    type(QcMat) :: W, A, B, C, T
    logical :: calc_contrib

    call QcMatInit(A)
    call QcMatInit(B)
    call QcMatInit(C)

    call QcMatInit(T)
    
    do i = 1, superstructure_size

       ds1(1) = deriv_struct(i,1)
       ds2(1) = deriv_struct(i,2)
       ds3(1) = deriv_struct(i,3)

       allocate(iu1(ds1(1)%npert))
       allocate(iu2(ds2(1)%npert))
       allocate(iu3(ds3(1)%npert))

       iu1 = get_fds_data_index(ds1(1), total_num_perturbations, which_index_is_pid, &
            indices_len, ind)
       iu2 = get_fds_data_index(ds2(1), total_num_perturbations, which_index_is_pid, &
            indices_len, ind)
       iu3 = get_fds_data_index(ds3(1), total_num_perturbations, which_index_is_pid, &
            indices_len, ind)

       call contrib_cache_getdata_outer(len_d, D, 1, ds1, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu1, mat_sing=A) 
       call contrib_cache_getdata_outer(len_f, F, 1, ds2, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu2, mat_sing=B)                
       call contrib_cache_getdata_outer(len_d, D, 1, ds3, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu3, mat_sing=C)       

!      W = W + A * B * C
       call QcMatkABC(1.0d0, A, B, C, T)
       call QcMatRAXPY(1.0d0, T, W)
            

       calc_contrib = .not.find_residue_info(deriv_struct(i,2))

       if (calc_contrib) then
         if (.not.(frequency_zero_or_sum(deriv_struct(i,1)) == 0.0) .and. &
             .not.(frequency_zero_or_sum(deriv_struct(i,3)) == 0.0)) then

          call contrib_cache_getdata_outer(len_d, D, 1, ds1, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu1, mat_sing=A) 
          call contrib_cache_getdata_outer(len_s, S, 1, ds2, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu2, mat_sing=B)                
          call contrib_cache_getdata_outer(len_d, D, 1, ds3, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu3, mat_sing=C)               
               

!            W = W + ((1.0)/(2.0)) * (frequency_zero_or_sum(deriv_struct(i,1)) - &
!                                    frequency_zero_or_sum(deriv_struct(i,3))) * A * B * C
               
          call QcMatcABC(((1.0)/(2.0)) * (frequency_zero_or_sum(deriv_struct(i,1)) - &
                        frequency_zero_or_sum(deriv_struct(i,3))), A, B, C, T)
          call QcMatRAXPY(1.0d0, T, W)
               
         
         elseif (.not.(frequency_zero_or_sum(deriv_struct(i,1)) == 0.0) .and. &
                      (frequency_zero_or_sum(deriv_struct(i,3)) == 0.0)) then

          call contrib_cache_getdata_outer(len_d, D, 1, ds1, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu1, mat_sing=A) 
          call contrib_cache_getdata_outer(len_s, S, 1, ds2, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu2, mat_sing=B)                
          call contrib_cache_getdata_outer(len_d, D, 1, ds3, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu3, mat_sing=C)               

!           W = W + ((1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,1)) * A * B * C
               
          call QcMatcABC(((1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,1)), A, B, C, T)
          call QcMatRAXPY(1.0d0, T, W)
               


         elseif (.not.(frequency_zero_or_sum(deriv_struct(i,3)) == 0.0) .and. &
                      (frequency_zero_or_sum(deriv_struct(i,1)) == 0.0)) then
               
          call contrib_cache_getdata_outer(len_d, D, 1, ds1, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu1, mat_sing=A) 
          call contrib_cache_getdata_outer(len_s, S, 1, ds2, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu2, mat_sing=B)                
          call contrib_cache_getdata_outer(len_d, D, 1, ds3, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu3, mat_sing=C)                

!           W = W + ((-1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,3))  * A * B * C

          call QcMatcABC(((-1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,3)), A, B, C, T)
          call QcMatRAXPY(1.0d0, T, W)   
          
         end if
       end if

       deallocate(iu1)
       deallocate(iu2)
       deallocate(iu3)


    end do

    call QcMatDst(A)
    call QcMatDst(B)
    call QcMatDst(C)
    call QcMatDst(T)

  end subroutine

  ! Calculate a perturbed Y matrix
  subroutine rsp_get_matrix_y(superstructure_size, deriv_struct, &
           total_num_perturbations, which_index_is_pid, indices_len, &
           ind, len_f, F, len_d, D, len_s, S, Y, select_terms_arg)

    implicit none

    logical, optional :: select_terms_arg
    logical :: select_terms, calc_contrib
    integer :: i, total_num_perturbations, superstructure_size, indices_len, j
    type(p_tuple), dimension(superstructure_size, 3) :: deriv_struct
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    type(p_tuple), dimension(1) :: ds1, ds2, ds3
    integer, dimension(:), allocatable :: iu1, iu2, iu3
    integer, dimension(indices_len) :: ind
    integer :: len_f, len_d, len_s
    type(contrib_cache_outer), dimension(len_f) :: F
    type(contrib_cache_outer), dimension(len_d) :: D
    type(contrib_cache_outer), dimension(len_S) :: S
    type(QcMat) :: Y, A, B, C, T

    call QcMatInit(A)
    call QcMatInit(B)
    call QcMatInit(C)

    call QcMatInit(T)
    
    select_terms = .FALSE.
    
    if (present(select_terms_arg)) then
    
       select_terms = select_terms_arg
    
    end if
    
    
    do i = 1, superstructure_size

       if (.not.select_terms) then
         calc_contrib = .true.
       else 
         calc_contrib = .not.find_residue_info(deriv_struct(i,3))
       end if

       ds1(1) = deriv_struct(i,1)
       ds2(1) = deriv_struct(i,2)
       ds3(1) = deriv_struct(i,3)

       allocate(iu1(ds1(1)%npert))
       allocate(iu2(ds2(1)%npert))
       allocate(iu3(ds3(1)%npert))

       iu1 = get_fds_data_index(ds1(1), total_num_perturbations, which_index_is_pid, &
            indices_len, ind)
       iu2 = get_fds_data_index(ds2(1), total_num_perturbations, which_index_is_pid, &
            indices_len, ind)
       iu3 = get_fds_data_index(ds3(1), total_num_perturbations, which_index_is_pid, &
            indices_len, ind)



       call contrib_cache_getdata_outer(len_f, F, 1, ds1, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu1, mat_sing=A)             
       call contrib_cache_getdata_outer(len_d, D, 1, ds2, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu2, mat_sing=B)             
       call contrib_cache_getdata_outer(len_s, S, 1, ds3, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu3, mat_sing=C) 

!        Y = Y + A*B*C
       if (calc_contrib) then
         call QcMatkABC(1.0d0, A, B, C, T)
         call QcMatRAXPY(1.0d0, T, Y)
       end if

       if (.not.(frequency_zero_or_sum(deriv_struct(i,2)) == 0.0)) then
             
          if (.not.select_terms) then 
            calc_contrib = .true.
          else
            calc_contrib = .not.(find_residue_info(deriv_struct(i,1)).or.find_residue_info(deriv_struct(i,3)))
          end if
          call contrib_cache_getdata_outer(len_s, S, 1, ds1, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu1, mat_sing=A)

!           Y = Y - frequency_zero_or_sum(deriv_struct(i,2)) * A * B * C

          if (calc_contrib) then
            call QcMatcABC(-1.0 * frequency_zero_or_sum(deriv_struct(i,2)), A, B, C, T)
            call QcMatRAXPY(1.0d0, T, Y)
          end if
          
       end if

       if (.not.select_terms) then
          calc_contrib = .true.
       else 
          calc_contrib = .not.find_residue_info(deriv_struct(i,1))
       end if
       call contrib_cache_getdata_outer(len_s, S, 1, ds1, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu1, mat_sing=A)  
       call contrib_cache_getdata_outer(len_f, F, 1, ds3, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu3, mat_sing=C)  
                       

!        Y = Y - A * B * C
       if (calc_contrib) then 
         call QcMatkABC(-1.0d0, A, B, C, T)
         call QcMatRAXPY(1.0d0, T, Y)
       end if

       if (.not.select_terms) then
         calc_contrib = .true.
       else
         calc_contrib = .not.(find_residue_info(deriv_struct(i,1)).or.find_residue_info(deriv_struct(i,3)))
       end if
       if (.not.(frequency_zero_or_sum(deriv_struct(i,1)) == 0.0) .and. &
           .not.(frequency_zero_or_sum(deriv_struct(i,3)) == 0.0)) then
          ! MaR: MAKE SURE THAT THESE (AND B) ARE ACTUALLY THE CORRECT 
          ! MATRICES TO USE HERE AND BELOW

          call contrib_cache_getdata_outer(len_s, S, 1, ds1, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu1, mat_sing=A)  
          call contrib_cache_getdata_outer(len_s, S, 1, ds3, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu3, mat_sing=C)  
       
               
!           Y = Y + ((-1.0)/(2.0)) * (frequency_zero_or_sum(deriv_struct(i,1)) + &
!                                     frequency_zero_or_sum(deriv_struct(i,3))) * A * B * C
                                    
          if (calc_contrib) then
            call QcMatcABC(((-1.0)/(2.0)) * (frequency_zero_or_sum(deriv_struct(i,1)) + &
                                    frequency_zero_or_sum(deriv_struct(i,3))), A, B, C, T)
            call QcMatRAXPY(1.0d0, T, Y)
          end if

       elseif (.not.(frequency_zero_or_sum(deriv_struct(i,1)) == 0.0) .and. &
             (frequency_zero_or_sum(deriv_struct(i,3)) == 0.0)) then
               
          call contrib_cache_getdata_outer(len_s, S, 1, ds1, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu1, mat_sing=A)  
          call contrib_cache_getdata_outer(len_s, S, 1, ds3, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu3, mat_sing=C)  

!           Y = Y + ((-1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,1)) * A * B * C
          if (calc_contrib) then 
            call QcMatcABC(((-1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,1)), A, B, C, T)
            call QcMatRAXPY(1.0d0, T, Y)
          end if

       elseif (.not.(frequency_zero_or_sum(deriv_struct(i,3)) == 0.0) .and. &
                    (frequency_zero_or_sum(deriv_struct(i,1)) == 0.0)) then
               
          call contrib_cache_getdata_outer(len_s, S, 1, ds1, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu1, mat_sing=A)  
          call contrib_cache_getdata_outer(len_s, S, 1, ds3, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu3, mat_sing=C)  
               
!           Y = Y + ((-1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,3)) * A * B * C
          if (calc_contrib) then
            call QcMatcABC(((-1.0)/(2.0)) * frequency_zero_or_sum(deriv_struct(i,3)), A, B, C, T)
            call QcMatRAXPY(1.0d0, T, Y)
          end if

       end if

       deallocate(iu1)
       deallocate(iu2)
       deallocate(iu3)

    end do

    call QcMatDst(A)
    call QcMatDst(B)
    call QcMatDst(C)
    call QcMatDst(T)


  end subroutine

  ! Calculate a perturbed Z matrix
  subroutine rsp_get_matrix_z(superstructure_size, deriv_struct, kn, &
           total_num_perturbations, which_index_is_pid, indices_len, &
           ind, len_f, F, len_d, D, len_s, S, Z, select_terms_arg)

    implicit none

    logical, optional :: select_terms_arg
    logical :: select_terms, calc_contrib
    integer :: i, total_num_perturbations, superstructure_size, indices_len
    type(p_tuple), dimension(superstructure_size, 3) :: deriv_struct
    type(p_tuple), dimension(1) :: merged_p_tuple
    type(p_tuple), dimension(1) :: ds1, ds2, ds3
    integer, dimension(:), allocatable :: iu1, iu2, iu3, ium
    integer, dimension(2) :: kn
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: ind
    integer :: len_f, len_d, len_s
    type(contrib_cache_outer), dimension(len_f) :: F
    type(contrib_cache_outer), dimension(len_d) :: D
    type(contrib_cache_outer), dimension(len_S) :: S
    type(QcMat) :: Z, A, B, C, T

    call QcMatInit(A)
    call QcMatInit(B)
    call QcMatInit(C)

    call QcMatInit(T)
 
    select_terms = .FALSE.
    
    if (present(select_terms_arg)) then
    
       select_terms = select_terms_arg
    
    end if
 
 

    do i = 1, superstructure_size

       ds1(1) = deriv_struct(i,1)
       ds2(1) = deriv_struct(i,2)
       ds3(1) = deriv_struct(i,3)

       allocate(iu1(ds1(1)%npert))
       allocate(iu2(ds2(1)%npert))
       allocate(iu3(ds3(1)%npert))

       iu1 = get_fds_data_index(ds1(1), total_num_perturbations, which_index_is_pid, &
            indices_len, ind)
       iu2 = get_fds_data_index(ds2(1), total_num_perturbations, which_index_is_pid, &
            indices_len, ind)
       iu3 = get_fds_data_index(ds3(1), total_num_perturbations, which_index_is_pid, &
            indices_len, ind)


       if (.not.select_terms) then
         calc_contrib = .true.
       else
         calc_contrib = .not.find_residue_info(deriv_struct(i,2))
       end if

       call contrib_cache_getdata_outer(len_d, D, 1, ds1, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu1, mat_sing=A) 
       call contrib_cache_getdata_outer(len_s, S, 1, ds2, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu2, mat_sing=B)             
       call contrib_cache_getdata_outer(len_d, D, 1, ds3, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu3, mat_sing=C)
            
!        Z = Z + A*B*C
       if (calc_contrib) then 
         call QcMatkABC(1.0d0, A, B, C, T)
         call QcMatRAXPY(1.0d0, T, Z)
       end if

       deallocate(iu1)
       deallocate(iu2)
       deallocate(iu3)

    end do

    merged_p_tuple(1) = merge_p_tuple(deriv_struct(1,1), merge_p_tuple(deriv_struct(1,2), deriv_struct(1,3)))

    allocate(ium(total_num_perturbations))

    ium = (/get_fds_data_index(merged_p_tuple(1), &
        total_num_perturbations, which_index_is_pid, indices_len, ind)/)

    if (kn_skip(total_num_perturbations, merged_p_tuple(1)%pid, kn) .eqv. .FALSE.) then

       call contrib_cache_getdata_outer(len_d, D, 1, merged_p_tuple, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=ium, mat_sing=A) 
            
!         Z = Z - A
        call QcMatRAXPY(-1.0d0, A, Z)

    end if

    deallocate(ium)

    call p_tuple_deallocate(merged_p_tuple(1))

    call QcMatDst(A)
    call QcMatDst(B)
    call QcMatDst(C)
    call QcMatDst(T)


  end subroutine

  ! Calculate a perturbed Lambda matrix
  subroutine rsp_get_matrix_lambda(p_tuple_a, superstructure_size, deriv_struct, &
           total_num_perturbations, which_index_is_pid, indices_len, &
           ind, len_d, D, len_s, S, L,&
           select_terms_arg)

    implicit none

    integer :: i, total_num_perturbations, superstructure_size, indices_len
    type(p_tuple) :: p_tuple_a
    type(p_tuple), dimension(1) :: merged_A, merged_B
    type(p_tuple), dimension(superstructure_size, 3) :: deriv_struct
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: ind
    type(p_tuple), dimension(1) :: ds1, ds2, ds3
    integer, dimension(:), allocatable :: iu1, iu2, iu3, iuma, iumb
    integer :: len_d, len_s
    type(contrib_cache_outer), dimension(len_d) :: D
    type(contrib_cache_outer), dimension(len_S) :: S
    type(QcMat) :: L, A, B, C, T
    logical :: calc_contrib, select_terms
    logical, optional :: select_terms_arg

    call QcMatInit(A)
    call QcMatInit(B)
    call QcMatInit(C)

    call QcMatInit(T)

    calc_contrib = .true.
    
    select_terms = .FALSE.
    
    if (present(select_terms_arg)) then
    
       select_terms = select_terms_arg
    
    end if
    

    do i = 1, superstructure_size

       ds1(1) = deriv_struct(i,1)
       ds2(1) = deriv_struct(i,2)
       ds3(1) = deriv_struct(i,3)

       allocate(iu1(ds1(1)%npert))
       allocate(iu2(ds2(1)%npert))
       allocate(iu3(ds3(1)%npert))

       iu1 = get_fds_data_index(ds1(1), total_num_perturbations, which_index_is_pid, &
            indices_len, ind)
       iu2 = get_fds_data_index(ds2(1), total_num_perturbations, which_index_is_pid, &
            indices_len, ind)
       iu3 = get_fds_data_index(ds3(1), total_num_perturbations, which_index_is_pid, &
            indices_len, ind)

            
       if (select_terms) then
      
          calc_contrib = .not.find_residue_info(deriv_struct(i,2))
         
       end if
           

      if (calc_contrib) then
      

       merged_A(1) = merge_p_tuple(p_tuple_a, deriv_struct(i,1))
       merged_B(1) = merge_p_tuple(p_tuple_a, deriv_struct(i,3))

       allocate(iuma(merged_A(1)%npert))
       allocate(iumb(merged_B(1)%npert))

       iuma = (/get_fds_data_index(merged_A(1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind)/)

       iumb = (/get_fds_data_index(merged_B(1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind)/)

       call contrib_cache_getdata_outer(len_d, D, 1, merged_A, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iuma, mat_sing=A) 
       call contrib_cache_getdata_outer(len_s, S, 1, ds2, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu2, mat_sing=B) 
       call contrib_cache_getdata_outer(len_d, D, 1, ds3, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu3, mat_sing=C)    
            
!        L = L + A * B * C
       call QcMatkABC(1.0d0, A, B, C, T)
       call QcMatRAXPY(1.0d0, T, L)            
            
       call contrib_cache_getdata_outer(len_d, D, 1, ds1, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu1, mat_sing=A) 
       call contrib_cache_getdata_outer(len_d, D, 1, merged_B, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iumb, mat_sing=C)             

!        L = L - A * B * C
       call QcMatkABC(-1.0d0, A, B, C, T)
       call QcMatRAXPY(1.0d0, T, L)    

       call p_tuple_deallocate(merged_A(1))
       call p_tuple_deallocate(merged_B(1))

       deallocate(iuma)
       deallocate(iumb)

      end if

       deallocate(iu1)
       deallocate(iu2)
       deallocate(iu3)
       
    end do
    
    
    call QcMatDst(A)
    call QcMatDst(B)
    call QcMatDst(C)
    call QcMatDst(T)


  end subroutine

  ! Calculate a perturbed Zeta matrix
  subroutine rsp_get_matrix_zeta(p_tuple_a, kn, superstructure_size, deriv_struct, &
           total_num_perturbations, which_index_is_pid, indices_len, &
           ind, len_f, F, len_d, D, len_s, S, Zeta, select_terms_arg)

    implicit none

    integer :: i, total_num_perturbations, superstructure_size, indices_len, j
    type(p_tuple) :: p_tuple_a
    type(p_tuple), dimension(1) :: merged_p_tuple, merged_A, merged_B
    type(p_tuple), dimension(superstructure_size, 3) :: deriv_struct
    integer, dimension(2) :: kn
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: ind
    integer :: len_f, len_d, len_s
    type(p_tuple), dimension(1) :: ds1, ds2, ds3
    integer, dimension(:), allocatable :: iu1, iu2, iu3, ium, iuma, iumb
    type(contrib_cache_outer), dimension(len_f) :: F
    type(contrib_cache_outer), dimension(len_d) :: D
    type(contrib_cache_outer), dimension(len_S) :: S
    type(QcMat) :: Zeta, A, B, C, T
    logical :: calc_contrib, select_terms
    logical, optional :: select_terms_arg

    call QcMatInit(A)
    call QcMatInit(B)
    call QcMatInit(C)

    call QcMatInit(T)
    
    calc_contrib = .true.
    
    select_terms = .FALSE.
    
    if (present(select_terms_arg)) then
    
       select_terms = select_terms_arg
    
    end if

    do i = 1, superstructure_size

       ds1(1) = deriv_struct(i,1)
       ds2(1) = deriv_struct(i,2)
       ds3(1) = deriv_struct(i,3)

       allocate(iu1(ds1(1)%npert))
       allocate(iu2(ds2(1)%npert))
       allocate(iu3(ds3(1)%npert))

       merged_A(1) = merge_p_tuple(p_tuple_a, deriv_struct(i,1))
       merged_B(1) = merge_p_tuple(p_tuple_a, deriv_struct(i,3))

       allocate(iuma(merged_A(1)%npert))
       allocate(iumb(merged_B(1)%npert))
       
       iuma = (/get_fds_data_index(merged_A(1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind)/)

       iumb = (/get_fds_data_index(merged_B(1), &
            total_num_perturbations, which_index_is_pid, indices_len, ind)/)

       iu1 = get_fds_data_index(ds1(1), total_num_perturbations, which_index_is_pid, &
            indices_len, ind)
       iu2 = get_fds_data_index(ds2(1), total_num_perturbations, which_index_is_pid, &
            indices_len, ind)
       iu3 = get_fds_data_index(ds3(1), total_num_perturbations, which_index_is_pid, &
            indices_len, ind)

     if (select_terms) then 
     
        calc_contrib = .not.find_residue_info(deriv_struct(i,3))
        
     end if

     if (calc_contrib) then     
    

       call contrib_cache_getdata_outer(len_f, F, 1, merged_A, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iuma, mat_sing=A)   
       call contrib_cache_getdata_outer(len_d, D, 1, ds2, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu2, mat_sing=B)
       call contrib_cache_getdata_outer(len_s, S, 1, ds3, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu3, mat_sing=C)            
     
!        Zeta = Zeta + A * B * C
       call QcMatkABC(1.0d0, A, B, C, T)
       call QcMatRAXPY(1.0d0, T, Zeta)   

     end if
     
     if (select_terms) then
     
        calc_contrib = .not.find_residue_info(merged_B(1))
        
     end if

     if (calc_contrib) then
 
       call contrib_cache_getdata_outer(len_f, F, 1, ds1, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu1, mat_sing=A)          
       call contrib_cache_getdata_outer(len_s, S, 1, merged_B, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iumb, mat_sing=C)       
            
!        Zeta = Zeta - A * B * C

       call QcMatkABC(-1.0d0, A, B, C, T)
       call QcMatRAXPY(1.0d0, T, Zeta)   

     end if

     if (select_terms) then 
     
        calc_contrib = find_residue_info(deriv_struct(i,2))
        
     end if

     if (calc_contrib) then

       if (.not.(frequency_zero_or_sum(deriv_struct(i,1)) == 0.0) .and. &
           .not.(frequency_zero_or_sum(deriv_struct(i,2)) == 0.0)) then
               
          call contrib_cache_getdata_outer(len_s, S, 1, ds1, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu1, mat_sing=A)                     

!           Zeta = Zeta + ( ((1.0)/(2.0))*frequency_zero_or_sum(deriv_struct(i,1)) + &
!                            frequency_zero_or_sum(deriv_struct(i,2)) ) * A * B * C
          call QcMatcABC(((1.0d0)/(2.0d0))*frequency_zero_or_sum(deriv_struct(i,1)) + &
                           frequency_zero_or_sum(deriv_struct(i,2)), A, B, C, T)
          call QcMatRAXPY(1.0d0, T, Zeta)   

       elseif (.not.(frequency_zero_or_sum(deriv_struct(i,1)) == 0.0) .and. &
                    (frequency_zero_or_sum(deriv_struct(i,2)) == 0.0)) then

          call contrib_cache_getdata_outer(len_s, S, 1, ds1, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu1, mat_sing=A)  

!           Zeta = Zeta + ((1.0)/(2.0))*frequency_zero_or_sum(deriv_struct(i,1)) * A * B * C
          call QcMatcABC(((1.0d0)/(2.0d0))*frequency_zero_or_sum(deriv_struct(i,1)), A, B, C, T)
          call QcMatRAXPY(1.0d0, T, Zeta)   

       elseif (.not.(frequency_zero_or_sum(deriv_struct(i,2)) == 0.0) .and. &
                    (frequency_zero_or_sum(deriv_struct(i,1)) == 0.0)) then

          call contrib_cache_getdata_outer(len_s, S, 1, ds1, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu1, mat_sing=A)                 

!           Zeta = Zeta + frequency_zero_or_sum(deriv_struct(i,2)) * A * B * C
          call QcMatcABC(frequency_zero_or_sum(deriv_struct(i,2)), A, B, C, T)
          call QcMatRAXPY(1.0d0, T, Zeta)  

       end if

     end if

     if (select_terms) then
     
        calc_contrib = .not.find_residue_info(deriv_struct(i,1))
    
     end if

     if (calc_contrib) then

       call contrib_cache_getdata_outer(len_s, S, 1, ds1, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu1, mat_sing=A) 
       call contrib_cache_getdata_outer(len_d, D, 1, ds2, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu2, mat_sing=B) 
       call contrib_cache_getdata_outer(len_f, F, 1, merged_B, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iumb, mat_sing=C) 
            
!        Zeta = Zeta +  A * B * C
       call QcMatkABC(1.0d0, A, B, C, T)
       call QcMatRAXPY(1.0d0, T, Zeta)  

     end if

     if (select_terms) then 
     
        calc_contrib = .not.find_residue_info(merged_A(1)) 
        
     end if

      if (calc_contrib) then
                 
       call contrib_cache_getdata_outer(len_s, S, 1, merged_A, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iuma, mat_sing=A) 
       call contrib_cache_getdata_outer(len_f, F, 1, ds3, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=iu3, mat_sing=C) 

!        Zeta = Zeta - A * B * C
       call QcMatkABC(1.0d0, A, B, C, T)
       call QcMatRAXPY(-1.0d0, T, Zeta)  

      end if

      if (select_terms) then
      
         calc_contrib = .not.find_residue_info(deriv_struct(i,2))
         
      end if

      if (calc_contrib) then

       if (.not.(frequency_zero_or_sum(deriv_struct(i,2)) == 0.0) .and. &
           .not.(frequency_zero_or_sum(deriv_struct(i,3)) == 0.0)) then
                    
          call contrib_cache_getdata_outer(len_s, S, 1, ds3, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu3, mat_sing=C) 

!           Zeta = Zeta - ( ((1.0)/(2.0))*frequency_zero_or_sum(deriv_struct(i,3)) + &
!                         frequency_zero_or_sum(deriv_struct(i,2)) ) * A * B * C
          call QcMatcABC(((-1.0d0)/(2.0d0))*frequency_zero_or_sum(deriv_struct(i,3)) + &
                        frequency_zero_or_sum(deriv_struct(i,2)), A, B, C, T)
          call QcMatRAXPY(1.0d0, T, Zeta)  

       elseif (.not.(frequency_zero_or_sum(deriv_struct(i,3)) == 0.0) .and. &
                    (frequency_zero_or_sum(deriv_struct(i,2)) == 0.0)) then

          call contrib_cache_getdata_outer(len_s, S, 1, ds3, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu3, mat_sing=C) 

!           Zeta = Zeta - ((1.0)/(2.0))*frequency_zero_or_sum(deriv_struct(i,3)) * A * B * C
          call QcMatcABC(((-1.0d0)/(2.0d0))*frequency_zero_or_sum(deriv_struct(i,3)), A, B, C, T)
          call QcMatRAXPY(1.0d0, T, Zeta)  

       elseif (.not.(frequency_zero_or_sum(deriv_struct(i,2)) == 0.0) .and. &
                    (frequency_zero_or_sum(deriv_struct(i,3)) == 0.0)) then

          call contrib_cache_getdata_outer(len_s, S, 1, ds3, .FALSE., contrib_size=1, &
               ind_len=indices_len, ind_unsorted=iu3, mat_sing=C) 

!           Zeta = Zeta - frequency_zero_or_sum(deriv_struct(i,2)) * A * B * C
          call QcMatcABC(-1.0d0 * frequency_zero_or_sum(deriv_struct(i,2)), A, B, C, T)
          call QcMatRAXPY(1.0d0, T, Zeta) 

       end if

      end if

       call p_tuple_deallocate(merged_A(1))
       call p_tuple_deallocate(merged_B(1))

       deallocate(iu1)
       deallocate(iu2)
       deallocate(iu3)

       deallocate(iuma)
       deallocate(iumb)
         
    end do

    merged_p_tuple(1) = merge_p_tuple(p_tuple_a, merge_p_tuple(deriv_struct(1,1), &
                     merge_p_tuple(deriv_struct(1,2), deriv_struct(1,3))))

    allocate(ium(total_num_perturbations))

    ium = (/get_fds_data_index(merged_p_tuple(1), &
       total_num_perturbations, which_index_is_pid, indices_len, ind)/)


    ! NOTE JAN 16: Does not seem likely that the skip condition will ever be met here: Any calculation
    ! of Zeta should already be "not skip" for this merged tuple
    if (kn_skip(merged_p_tuple(1)%npert, &
        merged_p_tuple(1)%pid, kn) .eqv. .FALSE.) then

       call contrib_cache_getdata_outer(len_f, F, 1, merged_p_tuple, .FALSE., contrib_size=1, &
            ind_len=indices_len, ind_unsorted=ium, mat_sing=A) 
            
!        Zeta = Zeta - A
       call QcMatRAXPY(-1.0d0, A, Zeta)         

    end if

    deallocate(ium)

    call p_tuple_deallocate(merged_p_tuple(1))

    call QcMatDst(A)
    call QcMatDst(B)
    call QcMatDst(C)
    call QcMatDst(T)


  end subroutine

  ! Get correct index for perturbed S, D, F matrix retrieval call
  function get_fds_data_index(pert_tuple, total_num_perturbations, which_index_is_pid, &
                              indices_len, indices)

    implicit none

    type(p_tuple) :: pert_tuple
    integer :: i, total_num_perturbations, indices_len
    integer, allocatable, dimension(:) :: get_fds_data_index
    integer, dimension(total_num_perturbations) :: which_index_is_pid
    integer, dimension(indices_len) :: indices

    if (pert_tuple%npert == 0) then

       allocate(get_fds_data_index(1))
       get_fds_data_index(1) = 1

    else

       allocate(get_fds_data_index(pert_tuple%npert))

       do i = 1, pert_tuple%npert
! FIXME: UNSURE ABOUT NEXT LINE                                                                              
! Will use just the regular index, but change back if this causes problems                         
          get_fds_data_index(i) = indices(i)

! NEXT LINE WAS THE ORIGINAL LINE, IT HAS GONE OUT OF BOUNDS AT LEAST ONCE
!          get_fds_data_index(i) = indices(which_index_is_pid(pert_tuple%pid(i)))          





       end do

    end if

  end function

  ! Return the sum of frequencies in perturbation tuple (or zero if unperturbed)
  function frequency_zero_or_sum(pert_tuple)

    implicit none

    type(p_tuple) :: pert_tuple
    complex(8) :: frequency_zero_or_sum
    integer :: i

    frequency_zero_or_sum = 0.0

    if (pert_tuple%npert > 0) then

       do i = 1, pert_tuple%npert

          frequency_zero_or_sum = frequency_zero_or_sum + pert_tuple%freq(i)

       end do

    end if

  end function

end module
