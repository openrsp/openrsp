! Copyright 2012 Magnus Ringholm
!
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

! Contains routines for the calculation of perturbed
! overlap, density and Fock matrices used throughout the
! rsp_general calculation.

module rsp_perturbed_sdf

  use matrix_defop
  use rsp_contribs
  use rsp_field_tuple
  use rsp_indices_and_addressing
  use rsp_perturbed_matrices
  use rsp_sdf_caching
  use rsp_lof_caching

  implicit none

  public rsp_fds
  public get_fds
  public rsp_fock_lowerorder
  public get_fock_lowerorder

!   type rsp_cfg
! 
!      type(matrix) :: zeromat
! 
!   end type

  contains

  recursive subroutine rsp_fds(zeromat, pert, kn, F, D, S)

    implicit none

    
    type(p_tuple) :: pert
    type(p_tuple), dimension(pert%n_perturbations) :: psub
    integer, dimension(2) :: kn
    integer :: i, j, k
    type(SDF) :: F, D, S
    type(matrix) :: zeromat

    ! Unless at final recursion level, recurse further
    ! Make all size (n - 1) subsets of the perturbations and recurse
    ! Then (at final recursion level) get perturbed F, D, S 

    if (pert%n_perturbations > 1) then

       call make_p_tuple_subset(pert, psub)

       do i = 1, size(psub)

          call rsp_fds(zeromat, psub(i), kn, F, D, S)

       end do       

    end if

    if (sdf_already(D, pert) .eqv. .FALSE.) then
         
       if (kn_skip(pert%n_perturbations, pert%pid, kn) .eqv. .FALSE.) then

          write(*,*) 'Calling ovlint/fock/density with labels ', pert%plab, &
                     ' and perturbation id ', pert%pid
          write(*,*) ' '
                 
          ! FIXME (MaR) Quick fix: Reenumerate pids from 1 and up so that 
          ! get_fds doesn't stumble. It seems to work, but consider rewrite.
         
          k = 1

          do j = 1, pert%n_perturbations

             pert%pid(j) = k
             k = k + 1

          end do

          call get_fds(zeromat, p_tuple_standardorder(pert), F, D, S)

       else

          write(*,*) 'Would have called ovlint/fock/density with labels ', &
                     pert%plab, ' and perturbation id ', pert%pid, &
                     ' but it was k-n forbidden'
          write(*,*) ' '

       end if

    else

       write(*,*) 'FDS for labels ', pert%plab, &
                  'and perturbation id ', pert%pid, ' was found in cache'
       write(*,*) ' '

    end if

  end subroutine



  ! ASSUMES THAT PERTURBATION TUPLE IS IN STANDARD ORDER
  subroutine get_fds(zeromat, pert, F, D, S)
#ifndef VAR_LSDALTON
    !host program specific solver
    use interface_rsp_solver, only: rsp_mosolver_exec
#endif
    implicit none

    
    integer :: sstr_incr, i, j, superstructure_size, nblks, perturbed_matrix_size
    integer, allocatable, dimension(:) :: ind, blk_sizes
    integer, allocatable, dimension(:,:) :: blk_info, indices
    integer, dimension(0) :: noc
    character(4), dimension(0) :: nof
    type(p_tuple) :: pert
    type(p_tuple), allocatable, dimension(:,:) :: derivative_structure
    type(SDF) :: F, D, S
    type(matrix) :: X(1), RHS(1), A, B, zeromat
    type(matrix), allocatable, dimension(:) :: Fp, Dp, Sp, Dh
    type(f_l_cache), pointer :: fock_lowerorder_cache


! ASSUME CLOSED SHELL
call mat_init(A, zeromat%nrow, zeromat%ncol, .true.)

! ASSUME CLOSED SHELL
call mat_init(B, zeromat%nrow, zeromat%ncol, .true.)

!     A = mat_alloc_like(zeromat)
!     A = mat_zero_like(zeromat)
!     call mat_ensure_alloc(A)
! 
!     B = mat_alloc_like(zeromat)
!     B = mat_zero_like(zeromat)
!     call mat_ensure_alloc(B)

    call sdf_getdata_s(D, get_emptypert(), (/1/), A)
    call sdf_getdata_s(S, get_emptypert(), (/1/), B)


    nblks = get_num_blks(pert)
    allocate(blk_info(nblks, 3))
allocate(blk_sizes(pert%n_perturbations))


    blk_info = get_blk_info(nblks, pert)
    
    perturbed_matrix_size = get_triangulated_size(nblks, blk_info)


blk_sizes = get_triangular_sizes(nblks, blk_info(:,2), blk_info(:,3))

! write(*,*) 'perturbed matrix size', perturbed_matrix_size


    allocate(Fp(perturbed_matrix_size))
    allocate(Dp(perturbed_matrix_size))
    allocate(Sp(perturbed_matrix_size))
    allocate(Dh(perturbed_matrix_size))

    ! Get the appropriate Fock/density/overlap matrices

    ! 1. Call ovlint and store perturbed overlap matrix

    call rsp_ovlint_tr(zeromat%nrow, pert%n_perturbations, pert%plab, &
                       (/ (1, j = 1, pert%n_perturbations) /), pert%pdim, &
                       nblks, blk_info, blk_sizes, &
                       perturbed_matrix_size, Sp)
    call sdf_add(S, pert, perturbed_matrix_size, Sp)

! write(*,*) 'Got Sp'
! do i = 1, perturbed_matrix_size
! 
! write(*,*) 'Sp at', i
!        write(*,*) Sp(i)%elms_alpha
! 
! end do

deallocate(blk_sizes)

    ! INITIALIZE AND STORE D INSTANCE WITH ZEROES
    ! THE ZEROES WILL ENSURE THAT TERMS INVOLVING THE HIGHEST ORDER DENSITY MATRICES
    ! WILL BE ZERO IN THE CONSTRUCTION OF Dp


!  write(*,*) 'zeroing'


    do i = 1, perturbed_matrix_size

!  write(*,*) '1'

! ASSUME CLOSED SHELL
call mat_init(Dp(i), zeromat%nrow, zeromat%ncol, .true.)

!        Dp(i) = mat_alloc_like(zeromat)
!        Dp(i) = mat_zero_like(zeromat)
!        call mat_ensure_alloc(Dp(i))
! write(*,*) 'Dp tag 1', Dp(i)%magic_tag
! write(*,*) '2'

! ASSUME CLOSED SHELL
call mat_init(Dh(i), zeromat%nrow, zeromat%ncol, .true.)

!        Dh(i) = mat_alloc_like(zeromat)
!        Dh(i) = mat_zero_like(zeromat)
!        call mat_ensure_alloc(Dh(i))
call mat_init(Fp(i), zeromat%nrow, zeromat%ncol, .true.)

! ! write(*,*) '3'
!        Fp(i) = mat_alloc_like(zeromat)
! ! write(*,*) '3a'
!        Fp(i) = mat_zero_like(zeromat)
! ! write(*,*) '3b'
!        call mat_ensure_alloc(Fp(i))
! ! write(*,*) '3c'

    end do

!  write(*,*) 'zeroed'
! write(*,*) 'Dp tag 2', Dp(1)%magic_tag
    call sdf_add(D, pert, perturbed_matrix_size, Dp)
! write(*,*) 'Dp tag 3', Dp(1)%magic_tag
!  write(*,*) 'zeroed Dp'

    ! 2. Construct Dp and the initial part of Fp
    ! a) For the initial part of Fp: Make the initial recursive (lower order) 
    ! oneint, twoint, and xcint calls as needed

    call f_l_cache_allocate(fock_lowerorder_cache)
! write(*,*) 'Dp tag 3b', Dp(1)%magic_tag

!  write(*,*) 'allocated f l cache'

    call rsp_fock_lowerorder(zeromat, pert, pert%n_perturbations, 1, (/get_emptypert()/), &
                         0, D, perturbed_matrix_size, Fp, fock_lowerorder_cache)


! write(*,*) 'Dp tag 4', Dp(1)%magic_tag

!  write(*,*) 'got fock lowerorder'
! 
! do i = 1, perturbed_matrix_size
! 
! write(*,*) 'Fp', i
! write(*,*) Fp(i)%elms_alpha
! 
! end do

    deallocate(fock_lowerorder_cache)


! do i = 1, perturbed_matrix_size
! 
! write(*,*) 'Fp again', i
! write(*,*) Fp(i)%elms_alpha
! 
! end do

    call sdf_add(F, pert, perturbed_matrix_size, Fp)

!  write(*,*) 'added to Fp'
! write(*,*) 'Dp tag 5', Dp(1)%magic_tag

    ! b) For Dp: Create differentiation superstructure: First dryrun for size, and
    ! then the actual superstructure call

    superstructure_size = derivative_superstructure_getsize(pert, &
                          (/pert%n_perturbations, pert%n_perturbations/), .FALSE., &
                          (/get_emptypert(), get_emptypert(), get_emptypert()/))

    sstr_incr = 0

! write(*,*) 'sstr size'

    allocate(derivative_structure(superstructure_size, 3))
    allocate(indices(perturbed_matrix_size, pert%n_perturbations))
    allocate(ind(pert%n_perturbations))

    call derivative_superstructure(pert, (/pert%n_perturbations, &
         pert%n_perturbations/), .FALSE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         superstructure_size, sstr_incr, derivative_structure)

! write(*,*) 'Dp tag 6', Dp(1)%magic_tag
! write(*,*) 'sstr'

    call make_triangulated_indices(nblks, blk_info, perturbed_matrix_size, indices)

! write(*,*) 'tri ind'

    do i = 1, size(indices, 1)

! write(*,*) 'i is', i

       ind = indices(i, :)

! write(*,*) 'index is', ind
!         write(*,*) 'Dp before z', Dp(i)%elms_alpha
! write(*,*) 'Dp tag', Dp(i)%magic_tag

       call rsp_get_matrix_z(zeromat, superstructure_size, derivative_structure, &
               (/pert%n_perturbations,pert%n_perturbations/), pert%n_perturbations, &
               (/ (j, j = 1, pert%n_perturbations) /), pert%n_perturbations, &
               ind, F, D, S, Dp(i))

! write(*,*) 'got z', Dp(i)%elms_alpha

       Dp(i) = Dp(i) - A * B * Dp(i) - Dp(i) * B * A

! write(*,*) 'projected dp'
!        write(*,*) 'Dp at projection', Dp(i)%elms_alpha

       call sdf_add(D, pert, perturbed_matrix_size, Dp)

! write(*,*) 'added sdf'

       ! 3. Complete the particular contribution to Fp
       ! NOTE (MaR): THERE SHOULD BE A CALL TO rsp_xcint HERE TOO
       ! THE IF CRITERION HERE NEEDS ANOTHER LOOK


! write(*,*) 'F particular before', Fp(i)%elms_alpha

!        if (pert%n_perturbations <=2) then
! 
          call rsp_twoint_tr(zeromat%nrow, 0, nof, noc, pert%pdim, Dp(i), &
                             1, Fp(i:i))

          call rsp_xcint_tr_adapt(zeromat%nrow, 0, nof, noc, pert%pdim, &
               (/ sdf_getdata(D, get_emptypert(), (/1/)), Dp(i) /) , &
                             1, Fp(i:i))
! 
!        end if




! write(*,*) 'did twoint_tr', Fp(i)%elms_alpha

       call sdf_add(F, pert, perturbed_matrix_size, Fp)


! write(*,*) 'F particular after', Fp(i)%elms_alpha

       ! 4. Make right-hand side using Dp

! ASSUME CLOSED SHELL
call mat_init(RHS(1), zeromat%nrow, zeromat%ncol, .true.)

!        RHS(1) = mat_alloc_like(zeromat)
!        RHS(1) = mat_zero_like(zeromat)
!        call mat_ensure_alloc(RHS(1))

       call rsp_get_matrix_y(zeromat, superstructure_size, derivative_structure, &
                pert%n_perturbations, (/ (j, j = 1, pert%n_perturbations) /), &
                pert%n_perturbations, ind, F, D, S, RHS(1))


! ASSUME CLOSED SHELL
call mat_init(X(1), zeromat%nrow, zeromat%ncol, .true.)
     
!        X(1) = mat_alloc_like(zeromat)
!        X(1) = mat_zero_like(zeromat)
!        call mat_ensure_alloc(X(1))

! write(*,*) 'made rhs', RHS(1)%elms_alpha

       ! Note (MaR): What does the second argument in rsp_mosolver_exec mean?
#ifndef VAR_LSDALTON
    !host program specific solver
       call rsp_mosolver_exec(RHS(1), (/0d0/), X)
#endif
       ! Note (MaR): Why multiply by -2 like below?
       X(1) = -2d0*X(1)
       RHS(1) = 0

! write(*,*) 'solved rsp eq'

       ! 5. Get Dh using the rsp equation solution X

       Dh(i) = X(1) * B * A - A * B * X(1)

       ! 6. Make homogeneous contribution to Fock matrix

       ! THE IF CRITERION HERE NEEDS ANOTHER LOOK

!        if (pert%n_perturbations <=2) then

          call rsp_twoint_tr(zeromat%nrow, 0, nof, noc, pert%pdim, Dh(i), &
                          1, Fp(i:i))

          call rsp_xcint_tr_adapt(zeromat%nrow, 0, nof, noc, pert%pdim, &
               (/ sdf_getdata(D, get_emptypert(), (/1/)), Dh(i) /) , &
                             1, Fp(i:i))

!        end if

       ! 7. Complete perturbed D with homogeneous part

       Dp(i) = Dp(i) + Dh(i)


write(*,*) 'Finished component', i

! write(*,*) 'Finally, Dp is', i, ' at indices', ind 
! write(*,*) Dp(i)%elms_alpha
! 
! write(*,*) 'Finally, Fp is', i
! write(*,*) Fp(i)%elms_alpha

    end do




    ! Add the final values to cache

    call sdf_add(F, pert, perturbed_matrix_size, Fp)
    call sdf_add(D, pert, perturbed_matrix_size, Dp)

!     do i = 1, size(indices, 1)
! 
! write(*,*) 'i is', i


! call sdf_getdata_s(D, pert, indices(i,:), A)
! call sdf_getdata_s(F, pert, indices(i,:), B)
! 
! write(*,*) 'Finally, in cache, Dp is', i
! write(*,*) A%elms_alpha
! 
! write(*,*) 'Finally, in cache, Fp is', i
! write(*,*) B%elms_alpha
! 
! end do

    do i = 1, size(indices, 1)
! write(*,*) 'deallocation 1'
Dh(i) = 0
Dp(i) = 0
Fp(i) = 0
Sp(i) = 0
! write(*,*) 'deallocation 2', i

end do


    deallocate(derivative_structure)
    deallocate(ind)

! write(*,*) 'deallocation 3'


    deallocate(Fp)
    deallocate(Dp)
    deallocate(Sp)
    deallocate(Dh)

    deallocate(blk_info)

! write(*,*) 'deallocation 4'

  end subroutine


  recursive subroutine rsp_fock_lowerorder(zeromat, pert, total_num_perturbations, &
                       num_p_tuples, p_tuples, density_order, D, property_size, Fp, &
                       fock_lowerorder_cache)

    implicit none

    logical :: density_order_skip
    
    type(p_tuple) :: pert
    integer :: num_p_tuples, density_order, i, j, total_num_perturbations, property_size
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, t_new
    type(SDF) :: D
    type(matrix) :: zeromat
    type(matrix), dimension(property_size) :: Fp
    type(f_l_cache) :: fock_lowerorder_cache


    if (pert%n_perturbations >= 1) then

       ! The differentiation can do three things:
       ! 1. Differentiate the expression 'directly'

       if (p_tuples(1)%n_perturbations == 0) then

          call rsp_fock_lowerorder(zeromat, p_tuple_remove_first(pert), & 
               total_num_perturbations, num_p_tuples, &
               (/p_tuple_getone(pert,1), p_tuples(2:size(p_tuples))/), &
               density_order, D, property_size, Fp, fock_lowerorder_cache)

       else

          call rsp_fock_lowerorder(zeromat, p_tuple_remove_first(pert), &
               total_num_perturbations, num_p_tuples, &
               (/p_tuple_extend(p_tuples(1), p_tuple_getone(pert,1)), &
               p_tuples(2:size(p_tuples))/), &
               density_order, D, property_size, Fp, fock_lowerorder_cache)

       end if
    
       ! 2. Differentiate all of the contraction densities in turn

       do i = 2, num_p_tuples

          t_new = p_tuples

          if (p_tuples(i)%n_perturbations == 0) then

             t_new(i) = p_tuple_getone(pert, 1)

          else

             t_new(i) = p_tuple_extend(t_new(i), p_tuple_getone(pert, 1))

          end if

          call rsp_fock_lowerorder(zeromat, p_tuple_remove_first(pert), &
               total_num_perturbations, num_p_tuples, &
               t_new, density_order + 1, D, property_size, Fp, fock_lowerorder_cache)

       end do

       ! 3. Chain rule differentiate w.r.t. the density (giving 
       ! a(nother) pert D contraction)

       call rsp_fock_lowerorder(zeromat, p_tuple_remove_first(pert), &
            total_num_perturbations, num_p_tuples + 1, &
            (/p_tuples(:), p_tuple_getone(pert, 1)/), &
            density_order + 1, D, property_size, Fp, fock_lowerorder_cache)

    else

       write(*,*) 'Getting perturbed Fock matrix lower order contribution'

       do i = 1, num_p_tuples
 
          if (i == 1) then

             write(*,*) 'F', p_tuples(i)%pid

          else

             write(*,*) 'D', p_tuples(i)%pid

          end if

       end do

       density_order_skip = .FALSE.

       do i = 2, num_p_tuples

          if (p_tuples(i)%n_perturbations >= total_num_perturbations) then

             density_order_skip = .TRUE.

          end if

       end do
      
       if (density_order_skip .EQV. .FALSE.) then

          if (f_l_cache_already(fock_lowerorder_cache, &
          num_p_tuples, p_tuples) .EQV. .FALSE.) then

             call get_fock_lowerorder(zeromat, num_p_tuples, total_num_perturbations, &
                                      p_tuples_standardorder(num_p_tuples, p_tuples), &
                                      density_order, D, property_size, Fp, &
                                      fock_lowerorder_cache)

             write(*,*) 'Calculated perturbed Fock matrix lower order contribution'
             write(*,*) ' '

! do i = 1, property_size
! write(*,*) 'LOF contrib. at i = ', i
! write(*,*) Fp(i)%elms_alpha
! 
! 
! end do

          else

             call f_l_cache_getdata(fock_lowerorder_cache, num_p_tuples, &
                                    p_tuples_standardorder(num_p_tuples, p_tuples), &
                                    property_size, Fp)

             write(*,*) ' '

! do i = 1, property_size
! write(*,*) 'LOF contrib. at i = ', i
! write(*,*) Fp(i)%elms_alpha
! 
! 
! end do

          end if

       else

          write(*,*) 'Skipping contribution: At least one contraction D perturbed' 
          write(*,*) 'at order for which perturbed D is to be found '
          write(*,*) ' '

       end if


    end if

  end subroutine




  subroutine get_fock_lowerorder(zeromat, num_p_tuples, total_num_perturbations, p_tuples, &
                                 density_order, D, property_size, Fp, &
                                 fock_lowerorder_cache)

    implicit none

    
    type(p_tuple) :: merged_p_tuple
    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    type(SDF) :: D
    type(matrix), allocatable, dimension(:) :: dens_tuple
    integer :: i, j, k, m, num_p_tuples, total_num_perturbations, merged_nblks, &
               density_order, property_size, fp_offset, lo_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, offset
    integer, dimension(0) :: noc
    integer, dimension(total_num_perturbations) :: ncarray, ncouter, ncinner, pidouter
    integer, allocatable, dimension(:) :: o_whichpert, o_whichpertbig, o_wh_forave
    integer, allocatable, dimension(:) :: ncoutersmall, pidoutersmall, ncinnersmall
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:) :: blk_sizes_merged
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    integer, allocatable, dimension(:,:) :: triang_indices_fp, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    type(matrix) :: zeromat
    type(matrix), allocatable, dimension(:) :: tmp, lower_order_contribution
    type(matrix), dimension(property_size) :: Fp
    type(f_l_cache) :: fock_lowerorder_cache
    
! write(*,*) 'start'
! 
! do i = 1, num_p_tuples
! 
! write(*,*) 'tuple', i, ':', p_tuples(i)%plab
! 
! 
! end do


    ncarray = get_ncarray(total_num_perturbations, num_p_tuples, p_tuples)
    ncouter = nc_only(total_num_perturbations, total_num_perturbations - & 
                      p_tuples(1)%n_perturbations, num_p_tuples - 1, &
                      p_tuples(2:num_p_tuples), ncarray)
    ncinner = nc_only(total_num_perturbations, p_tuples(1)%n_perturbations, 1, &
                      p_tuples(1), ncarray)

    allocate(ncoutersmall(total_num_perturbations - p_tuples(1)%n_perturbations))
    allocate(ncinnersmall(p_tuples(1)%n_perturbations))
    allocate(pidoutersmall(total_num_perturbations - p_tuples(1)%n_perturbations))

    ncoutersmall = nc_onlysmall(total_num_perturbations, total_num_perturbations - &
                                p_tuples(1)%n_perturbations, num_p_tuples - 1, &
                                p_tuples(2:num_p_tuples), ncarray)
    ncinnersmall = nc_onlysmall(total_num_perturbations, p_tuples(1)%n_perturbations, &
                   1, p_tuples(1), ncarray)
    pidoutersmall = get_pidoutersmall(total_num_perturbations - &
                    p_tuples(1)%n_perturbations, num_p_tuples - 1, &
                    p_tuples(2:num_p_tuples))

    allocate(o_whichpert(total_num_perturbations))
    allocate(o_wh_forave(total_num_perturbations))
    allocate(dens_tuple(num_p_tuples))

allocate(nfields(num_p_tuples))
allocate(nblks_tuple(num_p_tuples))




! Note: Second way of blks_tuple_info can in the general case be larger than
! needed, but is allocated this way to get a prismic data structure



! write(*,*) 'allocations: num_p_tuples:', num_p_tuples

do i = 1, num_p_tuples

! write(*,*) 'i is', i

nfields(i) = p_tuples(i)%n_perturbations
! write(*,*) 'nfields', nfields(i)

nblks_tuple(i) = get_num_blks(p_tuples(i))
! write(*,*) 'nblks tuple', nblks_tuple(i)


end do

allocate(blks_tuple_info(num_p_tuples, total_num_perturbations, 3))
allocate(blks_tuple_triang_size(num_p_tuples))


allocate(blk_sizes(num_p_tuples, total_num_perturbations))
allocate(blk_sizes_merged(total_num_perturbations))

do i = 1, num_p_tuples

blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p_tuples(i))

! write(*,*) 'blks_tuple_info', blks_tuple_info(i, :, :)
blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                            blks_tuple_info(i, 1:nblks_tuple(i), :))

blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))


! write(*,*) 'blks_tuple_triang_size(i)', blks_tuple_triang_size(i)
end do

! write(*,*) 'triang'

    outer_indices_size = product(blks_tuple_triang_size(2:num_p_tuples))


if (p_tuples(1)%n_perturbations == 0) then

    inner_indices_size = 1

else

    inner_indices_size = blks_tuple_triang_size(1)

end if

! write(*,*) 'index sizes'


    allocate(tmp(inner_indices_size))
    allocate(lower_order_contribution(inner_indices_size * outer_indices_size))

    o_whichpert = make_outerwhichpert(total_num_perturbations, num_p_tuples, p_tuples)



    call sortdimbypid(total_num_perturbations, total_num_perturbations - &
                      p_tuples(1)%n_perturbations, pidoutersmall, &
                      ncarray, ncoutersmall, o_whichpert)

! write(*,*) 'o_whichpert', o_whichpert

! write(*,*) 'near index loop'

    if (total_num_perturbations > p_tuples(1)%n_perturbations) then

k = 1

do i = 2, num_p_tuples

do j = 1, p_tuples(i)%n_perturbations


o_wh_forave(p_tuples(i)%pid(j)) = k

k = k + 1


end do



end do




!        do i = 1, size(o_whichpert)
! 
!           if (.NOT.(o_whichpert(i) == 0)) then
! 
!              o_wh_forave(o_whichpert(i)) = i
! 
!           end if
!   
!        end do


!  write(*,*) 'o_wh_forave', o_wh_forave
      do j = 1, size(lower_order_contribution)

! ASSUME CLOSED SHELL
call mat_init(lower_order_contribution(j), zeromat%nrow, zeromat%ncol, .true.)

!           lower_order_contribution(j) = mat_alloc_like(zeromat)
!           lower_order_contribution(j) = mat_zero_like(zeromat)
!           call mat_ensure_alloc(lower_order_contribution(j))

       end do


! write(*,*) '1b'
      do j = 1, size(tmp)

! ASSUME CLOSED SHELL
call mat_init(tmp(j), zeromat%nrow, zeromat%ncol, .true.)


!           tmp(j) = mat_alloc_like(zeromat)
!           tmp(j) = mat_zero_like(zeromat)
!           call mat_ensure_alloc(tmp(j))

       end do


! write(*,*) '1c'

       do i = 2, num_p_tuples
! write(*,*) 'i is', i
! write(*,*) 'size of dens tuple', size(dens_tuple)

! ASSUME CLOSED SHELL
call mat_init(dens_tuple(i), zeromat%nrow, zeromat%ncol, .true.)

!           dens_tuple(i) = mat_alloc_like(zeromat)
!           dens_tuple(i) = mat_zero_like(zeromat)
!           call mat_ensure_alloc(dens_tuple(i))


       end do

 
! write(*,*) '1d'
       allocate(outer_indices(outer_indices_size,size(ncoutersmall)))
       allocate(inner_indices(inner_indices_size,size(ncinnersmall)))
! write(*,*) '1'
! write(*,*) 'blks_tuple_info 2 a', blks_tuple_info(2, :, :)
! write(*,*) 'blks_tuple_info 3 a', blks_tuple_info(2, &
!                1:nblks_tuple(2), :)
       call make_triangulated_tuples_indices(num_p_tuples - 1, total_num_perturbations, & 
            nblks_tuple(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, &
            :, :), blks_tuple_triang_size(2:num_p_tuples), outer_indices)
! write(*,*) '2 2 2'
       if (p_tuples(1)%n_perturbations > 0) then

          call make_triangulated_indices(nblks_tuple(1), blks_tuple_info(1, &
               1:nblks_tuple(1), :), blks_tuple_triang_size(1), inner_indices)
! write(*,*) '3'
       end if

!        allocate(inner_offsets(product(ncinner)))

       do i = 1, size(outer_indices, 1)

! write(*,*) 'current outer index', outer_indices(i, :)

          do j = 2, num_p_tuples

             call sdf_getdata_s(D, p_tuples(j), (/ &
                             (outer_indices(i,o_wh_forave(p_tuples(j)%pid(k))), &
                             k = 1, p_tuples(j)%n_perturbations) /), dens_tuple(j))
! write(*,*) '4'
          end do

! NOTE: IS THIS ZEROING OF tmp REDUNDANT?

          do j = 1, size(tmp)

! ASSUME CLOSED SHELL
call mat_init(tmp(j), zeromat%nrow, zeromat%ncol, .true.)

!              tmp(j) = mat_alloc_like(zeromat)
!              tmp(j)%elms_alpha = 0.0
!              call mat_ensure_alloc(tmp(j))

          end do

! write(*,*) '5'
          if (num_p_tuples <= 1) then

             call rsp_oneint_tr(zeromat%nrow, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                             (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                             p_tuples(1)%pdim, nblks_tuple(1), blks_tuple_info(1, &
                   1:nblks_tuple(1), :), blk_sizes(1, 1:nblks_tuple(1)), size(tmp), tmp)

          end if


          if (num_p_tuples <= 2) then

             call rsp_twoint_tr(zeromat%nrow, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                             (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                             p_tuples(1)%pdim, dens_tuple(2), size(tmp), tmp)

          end if

! write(*,*) '6'

!              write(*,*) 'Calling xcint (NOTE: XCINT CALL SKIPPED FOR NOW)'
!
             call rsp_xcint_tr_adapt(zeromat%nrow, p_tuples(1)%n_perturbations, &
                           p_tuples(1)%plab, (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                           p_tuples(1)%pdim, &
                           (/ sdf_getdata(D, get_emptypert(), (/1/)), &
                           (dens_tuple(k), k = 2, num_p_tuples) /), &
                           property_size, tmp)


! Make tmp the correct ("inner-triangulated") size
! Make some storage array the size of the full amount of data to be cached
! Make sure that the data from the integral is put into the 
! appropriate positions in tmp
! After the accumulation loop: Add the appropriate elements of tmp to
! Fp and then add tmp to fock_lowerorder_cache



if (p_tuples(1)%n_perturbations > 0) then

do j = 1, size(inner_indices,1)

offset = get_triang_blks_tuple_offset(num_p_tuples, total_num_perturbations, nblks_tuple, &
         (/ (p_tuples(k)%n_perturbations, k = 1, num_p_tuples) /), &
         blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
         (/inner_indices(j, :), outer_indices(i, :) /)) 


! write(*,*) 'full index tuple is', (/inner_indices(j, :), outer_indices(i, :) /)
! write(*,*) 'offset is', offset

! ASSUME CLOSED SHELL
call mat_init(lower_order_contribution(offset), zeromat%nrow, zeromat%ncol, .true.)


! lower_order_contribution(offset) = mat_alloc_like(zeromat)
! lower_order_contribution(offset)%elms_alpha = 0.0
! call mat_ensure_alloc(lower_order_contribution(offset))


lower_order_contribution(offset) = tmp(j)

end do

else

! Note: There might be problems with this call (since the first p_tuple is empty)

offset = get_triang_blks_tuple_offset(num_p_tuples - 1, total_num_perturbations,  &
         nblks_tuple(2:num_p_tuples), &
         (/ (p_tuples(k)%n_perturbations, k = 2, num_p_tuples) /), &
         blks_tuple_info(2:num_p_tuples, :, :), blk_sizes(2:num_p_tuples,:), & 
         blks_tuple_triang_size(2:num_p_tuples), &
         (/outer_indices(i, :) /)) 

! ASSUME CLOSED SHELL
call mat_init(lower_order_contribution(offset), zeromat%nrow, zeromat%ncol, .true.)


! lower_order_contribution(offset) = mat_alloc_like(zeromat)
! lower_order_contribution(offset)%elms_alpha = 0.0
! call mat_ensure_alloc(lower_order_contribution(offset))

lower_order_contribution(offset) = tmp(1)

end if
! write(*,*) '6'
       end do

! get triangulated indices for Fp

if (p_tuples(1)%n_perturbations > 0) then

call p_tuple_p1_cloneto_p2(p_tuples(1), merged_p_tuple)

do i = 2, num_p_tuples

! This can be problematic - consider rewriting merge_p_tuple as subroutine

merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))

end do
! write(*,*) '7'
else
! write(*,*) '7a'
call p_tuple_p1_cloneto_p2(p_tuples(2), merged_p_tuple)
! write(*,*) '7b'
do i = 3, num_p_tuples

! This can be problematic - consider rewriting merge_p_tuple as subroutine

merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))


end do

end if
! write(*,*) '8'
merged_nblks = get_num_blks(merged_p_tuple)
allocate(merged_blk_info(1, merged_nblks, 3))
merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)

blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))

! do i = 1, merged_nblks
! write(*,*) 'merged block info', i, ':', merged_blk_info(1,i,:)
! end do

merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)

! write(*,*) 'about to allocate'
! write(*,*) 'dimensions:', merged_triang_size, sum(merged_blk_info(1, :,2))

allocate(triang_indices_fp(merged_triang_size, sum(merged_blk_info(1, :,2))))

! write(*,*) 'made allocation'

call make_triangulated_indices(merged_nblks, merged_blk_info, & 
     merged_triang_size, triang_indices_fp)

! write(*,*) '8b'
do i = 1, size(triang_indices_fp, 1)

fp_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
            (/sum(nfields)/), &
            (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
            (/triang_indices_fp(i, :) /))
! write(*,*) 'indices are', triang_indices_fp(i, :)

! write(*,*) 'fp offset is', fp_offset

if (p_tuples(1)%n_perturbations > 0) then

lo_offset = get_triang_blks_tuple_offset(num_p_tuples, &
            total_num_perturbations, nblks_tuple, &
            nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
            (/triang_indices_fp(i, :) /))


else

lo_offset = get_triang_blks_tuple_offset(num_p_tuples - 1, &
            total_num_perturbations, nblks_tuple(2:num_p_tuples), &
            nfields(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, :, :), &
            blk_sizes(2:num_p_tuples,:), blks_tuple_triang_size(2:num_p_tuples), & 
            (/triang_indices_fp(i, :) /))



end if


! write(*,*) 'lo offset is', lo_offset 
! write(*,*) '8d'
Fp(fp_offset) = Fp(fp_offset) + lower_order_contribution(lo_offset)
! write(*,*) '8e'
end do
! write(*,*) '9'
! write(*,*) 'prop size', property_size
! write(*,*) 'lo size', size(lower_order_contribution)
call f_l_cache_add_element(fock_lowerorder_cache, num_p_tuples, p_tuples, &
     inner_indices_size * outer_indices_size, lower_order_contribution)

! ! write(*,*) '10'
deallocate(merged_blk_info)
deallocate(triang_indices_fp)
    deallocate(outer_indices)
    deallocate(inner_indices)

    else

       if (num_p_tuples <= 1) then
! write(*,*) '10'
          call rsp_oneint_tr(zeromat%nrow, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                          (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                          p_tuples(1)%pdim, nblks_tuple(1), blks_tuple_info(1, &
                   1:nblks_tuple(1), :), blk_sizes(1, 1:nblks_tuple(1)), property_size, Fp)
! Waiting for developments in rsp_contribs
! NOTE: Find out if necessary ovlint/oneint in "outer indices case" above
! NOTE: Add (a) corresponding call(s) for the energy term contributions (see e.g. eqn.
! 209 in ajt_rsp)
!           call rsp_ovlint(zeromat%nrow, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
!                           (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
!                           p_tuples(1)%pdim, w = p_tuples(1)%freq, fock = Fp)

       end if
! write(*,*) '11'
       if (num_p_tuples <= 2) then
! write(*,*) '12'
          call rsp_twoint_tr(zeromat%nrow, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
               (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
               p_tuples(1)%pdim, sdf_getdata(D, get_emptypert(), (/1/) ), &
               property_size, Fp)

       end if
! write(*,*) '13'
! do i = 1, p_tuples(1)%pdim(1)
! 
! write(*,*) 'for element ', i
! write(*,*) Fp(i)%elms_alpha
! 
! end do

!        write(*,*) 'Calling xcint (NOTE: XCINT CALL SKIPPED FOR NOW)'

       call rsp_xcint_tr_adapt(zeromat%nrow, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                      (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                      p_tuples(1)%pdim, &
                      (/ sdf_getdata(D, get_emptypert(), (/1/)) /), &
                      property_size, Fp)


       ! NOTE: THERE IS NO NEED TO CACHE THE "ALL INNER" CONTRIBUTION
       ! It should be possible to just add it to Fp like already done above
       ! even with the extra complexity from the triangularization 

    end if

deallocate(nfields)
deallocate(nblks_tuple)
deallocate(blks_tuple_info)
deallocate(blks_tuple_triang_size)


deallocate(blk_sizes)
deallocate(blk_sizes_merged)


    deallocate(ncoutersmall)
    deallocate(ncinnersmall)
    deallocate(pidoutersmall)
    deallocate(o_whichpert)
    deallocate(o_wh_forave)
!     deallocate(dens_tuple)


    deallocate(tmp)
    deallocate(lower_order_contribution)

  end subroutine

end module
