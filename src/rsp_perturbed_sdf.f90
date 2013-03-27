! Copyright 2012 Magnus Ringholm
!
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

! Contains routines for the calculation of perturbed
! overlap, density and Fock matrices used throughout the
! rsp_general calculation.

module rsp_perturbed_sdf

  use matrix_defop
  use matrix_lowlevel, only: mat_init, mat_zero_like
  use rsp_contribs
  use rsp_field_tuple
  use rsp_indices_and_addressing
  use rsp_perturbed_matrices
  use rsp_sdf_caching
  use rsp_lof_caching
  use interface_2el

  implicit none

  public rsp_fds
  public get_fds
  public rsp_fock_lowerorder
  public get_fock_lowerorder

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

          if (sdf_already(D, psub(i)) .eqv. .FALSE.) then

             call rsp_fds(zeromat, psub(i), kn, F, D, S)

          end if

       end do       

    end if

    if (sdf_already(D, pert) .eqv. .FALSE.) then
         
       if (kn_skip(pert%n_perturbations, pert%pid, kn) .eqv. .FALSE.) then

          write(*,*) 'Calling ovlint/fock/density with labels ', pert%plab, &
                     ' and perturbation id ', pert%pid, ' with frequencies (real part)', &
                     real(pert%freq)
          write(*,*) ' '
                 
          ! MaR: Quick fix: Reenumerate pids from 1 and up so that 
          ! get_fds doesn't stumble. It seems to work, but consider rewrite.
         
          k = 1

          do j = 1, pert%n_perturbations

             pert%pid(j) = k
             k = k + 1

          end do

          call get_fds(zeromat, p_tuple_standardorder(pert), F, D, S)

       else

!           write(*,*) 'Would have called ovlint/fock/density with labels ', &
!                      pert%plab, ' and perturbation id ', pert%pid, &
!                      ' but it was k-n forbidden'
!           write(*,*) ' '

       end if

    else

!        write(*,*) 'FDS for labels ', pert%plab, &
!                   'and perturbation id ', pert%pid, ' was found in cache'
!        write(*,*) ' '

    end if

  end subroutine



  ! ASSUMES THAT PERTURBATION TUPLE IS IN STANDARD ORDER
  subroutine get_fds(zeromat, pert, F, D, S)

    use interface_rsp_solver, only: rsp_solver_exec
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
!     call mat_init(A, zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call mat_zero_like(zeromat, A)

    ! ASSUME CLOSED SHELL
!     call mat_init(B, zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call mat_zero_like(zeromat, B)

    call sdf_getdata_s(D, get_emptypert(), (/1/), A)
    call sdf_getdata_s(S, get_emptypert(), (/1/), B)

    nblks = get_num_blks(pert)

    allocate(blk_info(nblks, 3))
    allocate(blk_sizes(pert%n_perturbations))

    blk_info = get_blk_info(nblks, pert)
    perturbed_matrix_size = get_triangulated_size(nblks, blk_info)
    blk_sizes = get_triangular_sizes(nblks, blk_info(:,2), blk_info(:,3))

    allocate(Fp(perturbed_matrix_size))
    allocate(Dp(perturbed_matrix_size))
    allocate(Sp(perturbed_matrix_size))
    allocate(Dh(perturbed_matrix_size))

    ! Get the appropriate Fock/density/overlap matrices

    ! 1. Call ovlint and store perturbed overlap matrix


    do i = 1, perturbed_matrix_size

       ! ASSUME CLOSED SHELL
!        call mat_init(Sp(i), zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call mat_zero_like(zeromat, Sp(i))

    end do
! write(*,*) 'Sp a', Sp(1)%elms

    call rsp_ovlint(zeromat%nrow, pert%n_perturbations, pert%plab, &
                       (/ (1, j = 1, pert%n_perturbations) /), pert%pdim, &
                       nblks, blk_info, blk_sizes, &
                       perturbed_matrix_size, Sp)

! write(*,*) 'Sp b', Sp(1)%elms

    call sdf_add(S, pert, perturbed_matrix_size, Sp)

    deallocate(blk_sizes)

    ! INITIALIZE AND STORE D INSTANCE WITH ZEROES
    ! THE ZEROES WILL ENSURE THAT TERMS INVOLVING THE HIGHEST ORDER DENSITY MATRICES
    ! WILL BE ZERO IN THE CONSTRUCTION OF Dp

    do i = 1, perturbed_matrix_size

       ! ASSUME CLOSED SHELL
!        call mat_init(Dp(i), zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call mat_zero_like(zeromat, Dp(i))

!        call mat_init(Dh(i), zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call mat_zero_like(zeromat, Dh(i))

!        call mat_init(Fp(i), zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call mat_zero_like(zeromat, Fp(i))

    end do

    call sdf_add(D, pert, perturbed_matrix_size, Dp)

    ! 2. Construct Dp and the initial part of Fp
    ! a) For the initial part of Fp: Make the initial recursive (lower order) 
    ! oneint, twoint, and xcint calls as needed

! write(*,*) 'Fp a', Fp(1)%elms

    call f_l_cache_allocate(fock_lowerorder_cache)
    call rsp_fock_lowerorder(zeromat, pert, pert%n_perturbations, 1, (/get_emptypert()/), &
                         0, D, perturbed_matrix_size, Fp, fock_lowerorder_cache)

! write(*,*) 'Fp b', Fp(1)%elms

    deallocate(fock_lowerorder_cache)

    call sdf_add(F, pert, perturbed_matrix_size, Fp)

    ! b) For Dp: Create differentiation superstructure: First dryrun for size, and
    ! then the actual superstructure call

    superstructure_size = derivative_superstructure_getsize(pert, &
                          (/pert%n_perturbations, pert%n_perturbations/), .FALSE., &
                          (/get_emptypert(), get_emptypert(), get_emptypert()/))

    sstr_incr = 0

    allocate(derivative_structure(superstructure_size, 3))
    allocate(indices(perturbed_matrix_size, pert%n_perturbations))
    allocate(ind(pert%n_perturbations))

    call derivative_superstructure(pert, (/pert%n_perturbations, &
         pert%n_perturbations/), .FALSE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         superstructure_size, sstr_incr, derivative_structure)
    call make_triangulated_indices(nblks, blk_info, perturbed_matrix_size, indices)

    do i = 1, size(indices, 1)

       ind = indices(i, :)

! write(*,*) 'Dp 0', Dp(1)%elms

       call rsp_get_matrix_z(zeromat, superstructure_size, derivative_structure, &
               (/pert%n_perturbations,pert%n_perturbations/), pert%n_perturbations, &
               (/ (j, j = 1, pert%n_perturbations) /), pert%n_perturbations, &
               ind, F, D, S, Dp(i))

! write(*,*) 'Dp 1', Dp(1)%elms


       Dp(i) = Dp(i) - A * B * Dp(i) - Dp(i) * B * A
! write(*,*) 'Dp 2', Dp(1)%elms


       call sdf_add(D, pert, perturbed_matrix_size, Dp)

       ! 3. Complete the particular contribution to Fp
! write(*,*) 'Fp b2', Fp(1)%elms
       call rsp_twoint(zeromat%nrow, 0, nof, noc, pert%pdim, Dp(i), &
                          1, Fp(i:i))
! write(*,*) 'Fp b3', Fp(1)%elms
       call rsp_xcint_adapt(zeromat%nrow, 0, nof, noc, pert%pdim, &
            (/ A, Dp(i) /) , 1, Fp(i:i))
! write(*,*) 'Fp b4', Fp(1)%elms

       call sdf_add(F, pert, perturbed_matrix_size, Fp)
! write(*,*) 'Fp c', Fp(1)%elms


       ! 4. Make right-hand side using Dp

       ! ASSUME CLOSED SHELL
!        call mat_init(RHS(1), zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call mat_zero_like(zeromat, RHS(1))

!        call mat_init(X(1), zeromat%nrow, zeromat%ncol, is_zero=.true.)
    call mat_zero_like(zeromat, X(1))

! write(*,*) 'RHS a', RHS(1)%elms

       call rsp_get_matrix_y(zeromat, superstructure_size, derivative_structure, &
                pert%n_perturbations, (/ (j, j = 1, pert%n_perturbations) /), &
                pert%n_perturbations, ind, F, D, S, RHS(1))

! write(*,*) 'RHS b', RHS(1)%elms

       ! Note (MaR): Passing only real part of freq. Is this OK?
       call rsp_solver_exec(RHS(1), (/sum(real(pert%freq(:)))/), X)
       RHS(1) = 0


       ! 5. Get Dh using the rsp equation solution X

       Dh(i) = A*B*X(1) - X(1)*B*A

       ! 6. Make homogeneous contribution to Fock matrix

       call rsp_twoint(zeromat%nrow, 0, nof, noc, pert%pdim, Dh(i), &
                          1, Fp(i:i))

       call rsp_xcint_adapt(zeromat%nrow, 0, nof, noc, pert%pdim, &
            (/ A, Dh(i) /) , 1, Fp(i:i))

       ! 7. Complete perturbed D with homogeneous part

       Dp(i) = Dp(i) + Dh(i)



if (perturbed_matrix_size < 10) then

write(*,*) 'Finished component', i, ':', (float(i)/float(perturbed_matrix_size))*100, '% done'
  

else

if (mod(i, perturbed_matrix_size/10) == 1) then

write(*,*) 'Finished component', i, ':', (float(i)/float(perturbed_matrix_size))*100, '% done'

end if

end if

! write(*,*) perturbed_matrix_size/10
! 
! if (mod(i, ) == 1) then
! 
!
! 
! end if


!        write(*,*) ' '
!        write(*,*) 'Finally, Dp is:'
!        write(*,*) Dp(i)%elms
!        write(*,*) ' '
!        write(*,*) 'Finally, Fp is:'
!        write(*,*) Fp(i)%elms
!        write(*,*) ' '
!        write(*,*) 'Finally, Sp is:'
!        write(*,*) Sp(i)%elms_alpha
!        write(*,*) ' '

    end do

    ! Add the final values to cache

    call sdf_add(F, pert, perturbed_matrix_size, Fp)
    call sdf_add(D, pert, perturbed_matrix_size, Dp)

    do i = 1, size(indices, 1)

       Dh(i) = 0
       Dp(i) = 0
       Fp(i) = 0
       Sp(i) = 0

    end do

    A = 0
    B = 0

    deallocate(derivative_structure)
    deallocate(ind)
    deallocate(Fp)
    deallocate(Dp)
    deallocate(Sp)
    deallocate(Dh)
    deallocate(blk_info)

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

!        p_tuples = p_tuples_standardorder(num_p_tuples, p_tuples)

       density_order_skip = .FALSE.

       do i = 2, num_p_tuples

          if (p_tuples(i)%n_perturbations >= total_num_perturbations) then

             density_order_skip = .TRUE.

          end if

       end do
      
       if (density_order_skip .EQV. .FALSE.) then

          if (f_l_cache_already(fock_lowerorder_cache, &
          num_p_tuples, p_tuples_standardorder(num_p_tuples, p_tuples)) .EQV. .FALSE.) then

       write(*,*) 'Calculating perturbed Fock matrix lower order contribution'

       do i = 1, num_p_tuples
 
          if (i == 1) then

             write(*,*) 'F', p_tuples(i)%pid

          else

             write(*,*) 'D', p_tuples(i)%pid

          end if

       end do

             call get_fock_lowerorder(zeromat, num_p_tuples, total_num_perturbations, &
                                      p_tuples_standardorder(num_p_tuples, p_tuples), &
                                      density_order, D, property_size, Fp, &
                                      fock_lowerorder_cache)

             write(*,*) 'Calculated perturbed Fock matrix lower order contribution'
             write(*,*) ' '

          else

             call f_l_cache_getdata(fock_lowerorder_cache, num_p_tuples, &
                                    p_tuples_standardorder(num_p_tuples, p_tuples), &
                                    property_size, Fp)

!              write(*,*) ' '

          end if

       else

!           write(*,*) 'Skipping contribution: At least one contraction D perturbed' 
!           write(*,*) 'at order for which perturbed D is to be found '
!           write(*,*) ' '

       end if

    end if

  end subroutine




  subroutine get_fock_lowerorder(zeromat, num_p_tuples, total_num_perturbations, p_tuples, &
                                 density_order, D, property_size, Fp, &
                                 fock_lowerorder_cache)

    implicit none
    
    type(p_tuple) :: merged_p_tuple, t_matrix_bra, t_matrix_ket
    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    type(SDF) :: D
    type(matrix), allocatable, dimension(:) :: dens_tuple
    integer :: i, j, k, m, num_p_tuples, total_num_perturbations, merged_nblks, &
               density_order, property_size, fp_offset, lo_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, offset
    integer, dimension(0) :: noc
    integer, dimension(total_num_perturbations) :: ncarray, ncouter, ncinner, pidouter, &
                                                pids_current_contribution, translated_index
    integer, allocatable, dimension(:) :: o_whichpert, o_whichpertbig, o_wh_forave
    integer, allocatable, dimension(:) :: ncoutersmall, pidoutersmall, ncinnersmall
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:) :: blk_sizes_merged
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    integer, allocatable, dimension(:,:) :: triang_indices_fp, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    type(matrix) :: zeromat, D_unp
    type(matrix), allocatable, dimension(:) :: tmp, lower_order_contribution
    type(matrix), dimension(property_size) :: Fp
    type(f_l_cache) :: fock_lowerorder_cache

!    ncarray = get_ncarray(total_num_perturbations, num_p_tuples, p_tuples)
!    ncouter = nc_only(total_num_perturbations, total_num_perturbations - & 
!                      p_tuples(1)%n_perturbations, num_p_tuples - 1, &
!                      p_tuples(2:num_p_tuples), ncarray)
!    ncinner = nc_only(total_num_perturbations, p_tuples(1)%n_perturbations, 1, &
!                      p_tuples(1), ncarray)

    allocate(ncoutersmall(total_num_perturbations - p_tuples(1)%n_perturbations))
    allocate(ncinnersmall(p_tuples(1)%n_perturbations))
    allocate(pidoutersmall(total_num_perturbations - p_tuples(1)%n_perturbations))

!    ncoutersmall = nc_onlysmall(total_num_perturbations, total_num_perturbations - &
!                                p_tuples(1)%n_perturbations, num_p_tuples - 1, &
!                                p_tuples(2:num_p_tuples), ncarray)
!    ncinnersmall = nc_onlysmall(total_num_perturbations, p_tuples(1)%n_perturbations, &
!                   1, p_tuples(1), ncarray)
!    pidoutersmall = get_pidoutersmall(total_num_perturbations - &
!                    p_tuples(1)%n_perturbations, num_p_tuples - 1, &
 !                   p_tuples(2:num_p_tuples))

    ! MaR: Second way of blks_tuple_info can in the general case be larger than
    ! needed, but is allocated this way to get a prismic data structure
    allocate(blks_tuple_info(num_p_tuples, total_num_perturbations, 3))
    allocate(blks_tuple_triang_size(num_p_tuples))
    allocate(blk_sizes(num_p_tuples, total_num_perturbations))
    allocate(blk_sizes_merged(total_num_perturbations))
    allocate(o_whichpert(total_num_perturbations))
    allocate(o_wh_forave(total_num_perturbations))
    allocate(dens_tuple(num_p_tuples))
    allocate(nfields(num_p_tuples))
    allocate(nblks_tuple(num_p_tuples))

    do i = 1, num_p_tuples

       nfields(i) = p_tuples(i)%n_perturbations
       nblks_tuple(i) = get_num_blks(p_tuples(i))

    end do

    do i = 1, num_p_tuples

       call get_blk_info_s(nblks_tuple(i), p_tuples(i), blks_tuple_info(i, 1:nblks_tuple(i), :))

! write(*,*) blks_tuple_info(i, :, :)
! write(*,*) 'sanitized'
! write(*,*) blks_tuple_info(i, 1:nblks_tuple(i), :)

       blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                                   blks_tuple_info(i, 1:nblks_tuple(i), :))


! write(*,*)  blks_tuple_triang_size(i) 

       blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
       blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))

    end do

    outer_indices_size = product(blks_tuple_triang_size(2:num_p_tuples))

    if (p_tuples(1)%n_perturbations == 0) then

       inner_indices_size = 1

    else

       inner_indices_size = blks_tuple_triang_size(1)

    end if

    allocate(tmp(inner_indices_size))
    allocate(lower_order_contribution(inner_indices_size * outer_indices_size))

    o_whichpert = make_outerwhichpert(total_num_perturbations, num_p_tuples, p_tuples)
!    call sortdimbypid(total_num_perturbations, total_num_perturbations - &
!                      p_tuples(1)%n_perturbations, pidoutersmall, &
!                      ncarray, ncoutersmall, o_whichpert)

    call sdf_getdata_s(D, get_emptypert(), (/1/), D_unp)

    if (total_num_perturbations > p_tuples(1)%n_perturbations) then

       k = 1

       do i = 2, num_p_tuples
          do j = 1, p_tuples(i)%n_perturbations

             o_wh_forave(p_tuples(i)%pid(j)) = k
             k = k + 1

          end do
       end do

       do j = 1, size(lower_order_contribution)

          ! ASSUME CLOSED SHELL
          call mat_init(lower_order_contribution(j), zeromat%nrow, zeromat%ncol)
          call mat_init_like_and_zero(zeromat, lower_order_contribution(j))

       end do

       do j = 1, size(tmp)

          ! ASSUME CLOSED SHELL
          call mat_init(tmp(j), zeromat%nrow, zeromat%ncol)
          call mat_init_like_and_zero(zeromat, tmp(j))

       end do

       do i = 2, num_p_tuples

          ! ASSUME CLOSED SHELL
          call mat_init(dens_tuple(i), zeromat%nrow, zeromat%ncol)
          call mat_init_like_and_zero(zeromat, dens_tuple(i))

       end do

       allocate(outer_indices(outer_indices_size,total_num_perturbations - &
                p_tuples(1)%n_perturbations))
       allocate(inner_indices(inner_indices_size,p_tuples(1)%n_perturbations))

       call make_triangulated_tuples_indices(num_p_tuples - 1, total_num_perturbations, & 
            nblks_tuple(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, &
            :, :), blks_tuple_triang_size(2:num_p_tuples), outer_indices)

       if (p_tuples(1)%n_perturbations > 0) then

          call make_triangulated_indices(nblks_tuple(1), blks_tuple_info(1, &
               1:nblks_tuple(1), :), blks_tuple_triang_size(1), inner_indices)

       end if

       do i = 1, size(outer_indices, 1)

          do j = 2, num_p_tuples

             call sdf_getdata_s(D, p_tuples(j), (/ &
                             (outer_indices(i,o_wh_forave(p_tuples(j)%pid(k))), &
                             k = 1, p_tuples(j)%n_perturbations) /), dens_tuple(j))

          end do

          ! MaR: IS THIS ZEROING OF tmp REDUNDANT?

          do j = 1, size(tmp)

             ! ASSUME CLOSED SHELL
             call mat_init(tmp(j), zeromat%nrow, zeromat%ncol)
             call mat_init_like_and_zero(zeromat, tmp(j))

          end do

          if (num_p_tuples <= 1) then

             call rsp_oneint(zeromat%nrow, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                             (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                             p_tuples(1)%pdim, nblks_tuple(1), blks_tuple_info(1, &
                   1:nblks_tuple(1), :), blk_sizes(1, 1:nblks_tuple(1)), size(tmp), tmp)

          end if

          if (num_p_tuples <= 2) then

             call rsp_twoint(zeromat%nrow, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                             (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                             p_tuples(1)%pdim, dens_tuple(2), size(tmp), tmp)

          end if

          call rsp_xcint_adapt(zeromat%nrow, p_tuples(1)%n_perturbations, &
               p_tuples(1)%plab, (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
               p_tuples(1)%pdim, (/ D_unp, &
               (dens_tuple(k), k = 2, num_p_tuples) /), property_size, tmp)

          if (p_tuples(1)%n_perturbations > 0) then

             do j = 1, size(inner_indices,1)

                offset = get_triang_blks_tuple_offset(num_p_tuples, total_num_perturbations, &
                nblks_tuple, (/ (p_tuples(k)%n_perturbations, k = 1, num_p_tuples) /), &
                blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                (/inner_indices(j, :), outer_indices(i, :) /)) 

                ! ASSUME CLOSED SHELL
                call mat_init(lower_order_contribution(offset), zeromat%nrow, zeromat%ncol)
                call mat_init_like_and_zero(zeromat, lower_order_contribution(offset))

                lower_order_contribution(offset) = tmp(j)

             end do

          else

             ! MaR: There might be problems with this call (since the first p_tuple is empty)

             offset = get_triang_blks_tuple_offset(num_p_tuples - 1, total_num_perturbations, &
             nblks_tuple(2:num_p_tuples), &
             (/ (p_tuples(k)%n_perturbations, k = 2, num_p_tuples) /), &
             blks_tuple_info(2:num_p_tuples, :, :), blk_sizes(2:num_p_tuples,:), & 
             blks_tuple_triang_size(2:num_p_tuples), (/outer_indices(i, :) /)) 

             ! ASSUME CLOSED SHELL
             call mat_init(lower_order_contribution(offset), zeromat%nrow, zeromat%ncol)
             call mat_init_like_and_zero(zeromat, lower_order_contribution(offset))

             lower_order_contribution(offset) = tmp(1)

          end if

       end do

       if (p_tuples(1)%n_perturbations > 0) then

          call p_tuple_p1_cloneto_p2(p_tuples(1), merged_p_tuple)

          do i = 2, num_p_tuples

             ! MaR: This can be problematic - consider rewriting merge_p_tuple as subroutine
             merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))

          end do

       else

          call p_tuple_p1_cloneto_p2(p_tuples(2), merged_p_tuple)

          do i = 3, num_p_tuples

             ! MaR: This can be problematic - consider rewriting merge_p_tuple as subroutine
             merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))

          end do

       end if

       merged_p_tuple = p_tuple_standardorder(merged_p_tuple)

       k = 1
       do i = 1, num_p_tuples
          do j = 1, p_tuples(i)%n_perturbations
             pids_current_contribution(k) = p_tuples(i)%pid(j)
             k = k + 1
          end do
       end do

! write(*,*) 'merged plab', merged_p_tuple%plab

       merged_nblks = get_num_blks(merged_p_tuple)

! write(*,*) 'merged plab 2', merged_p_tuple%plab

       allocate(merged_blk_info(1, merged_nblks, 3))

! write(*,*) 'allocate OK', merged_p_tuple%plab

       call get_blk_info_s(merged_nblks, merged_p_tuple, merged_blk_info(1, :, :))

! write(*,*) 'merged plab 3', merged_p_tuple%plab
! 
! do i = 1, merged_nblks
! 
! write(*,*) 'merged block info', merged_blk_info(1,i,:)
! 
! end do

       blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
       merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
       merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)

       allocate(triang_indices_fp(merged_triang_size, sum(merged_blk_info(1, :,2))))

       call make_triangulated_indices(merged_nblks, merged_blk_info, & 
            merged_triang_size, triang_indices_fp)

! do i = 1, size(triang_indices_fp,1)
! 
! write(*,*) 'triang indices', triang_indices_fp(i,:)
! 
! end do
! 
! write(*,*) 'size Fp', size(Fp)
! write(*,*) 'size loc', size(lower_order_contribution)


       do i = 1, size(triang_indices_fp, 1)

          fp_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
                      (/sum(nfields)/), &
                      (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
                      (/triang_indices_fp(i, :) /))

          do j = 1, total_num_perturbations
    
             translated_index(j) = triang_indices_fp(i,pids_current_contribution(j))
    
          end do

          if (p_tuples(1)%n_perturbations > 0) then

             lo_offset = get_triang_blks_tuple_offset(num_p_tuples, &
                         total_num_perturbations, nblks_tuple, &
                         nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                         (/translated_index(:)/))

          else

             lo_offset = get_triang_blks_tuple_offset(num_p_tuples - 1, &
                         total_num_perturbations, nblks_tuple(2:num_p_tuples), &
                         nfields(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, :, :), &
                         blk_sizes(2:num_p_tuples,:), &
                         blks_tuple_triang_size(2:num_p_tuples), & 
                         (/translated_index(:)/))

          end if

!           write(*,*) 'ind', triang_indices_fp(i, :)
! write(*,*) 'fp_offset', fp_offset
! write(*,*) 'lo_offset', lo_offset

          Fp(fp_offset) = Fp(fp_offset) + lower_order_contribution(lo_offset)

       end do

       call f_l_cache_add_element(fock_lowerorder_cache, num_p_tuples, p_tuples, &
            inner_indices_size * outer_indices_size, lower_order_contribution)

       deallocate(merged_blk_info)
       deallocate(triang_indices_fp)
       deallocate(outer_indices)
       deallocate(inner_indices)

    else

       if (num_p_tuples <= 1) then

          call rsp_oneint(zeromat%nrow, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                          (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                          p_tuples(1)%pdim, nblks_tuple(1), blks_tuple_info(1, &
                   1:nblks_tuple(1), :), blk_sizes(1, 1:nblks_tuple(1)), property_size, Fp)

! NOTE: Find out if necessary ovlint/oneint in "outer indices case" above
! NOTE (Oct 12): Probably not unless some hidden density matrix dependence

          t_matrix_bra = get_emptypert()
          t_matrix_ket = get_emptypert()

          call rsp_ovlint_t_matrix(zeromat%nrow, p_tuples(1)%n_perturbations, p_tuples(1), &
                                   t_matrix_bra, t_matrix_ket, property_size, Fp)

       end if

       if (num_p_tuples <= 2) then

          call rsp_twoint(zeromat%nrow, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
               (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
               p_tuples(1)%pdim, D_unp, &
               property_size, Fp)

       end if

       call rsp_xcint_adapt(zeromat%nrow, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                      (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                      p_tuples(1)%pdim, &
                      (/ D_unp /), &
                      property_size, Fp)

       ! MaR: THERE IS NO NEED TO CACHE THE "ALL INNER" CONTRIBUTION
       ! It should be possible to just add it to Fp like already done above
       ! even with the extra complexity from the triangularization 

    end if

    D_unp = 0

    do i = 1, num_p_tuples
   
       dens_tuple(i) = 0
   
    end do

    do i = 1, size(tmp)

       tmp(i) = 0

    end do

    do i = 1, size(lower_order_contribution)

       lower_order_contribution(i) = 0

    end do

    deallocate(dens_tuple)


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
    deallocate(tmp)
    deallocate(lower_order_contribution)

    ! MaR: Why is the next line commented? Find out
!     deallocate(dens_tuple)

  end subroutine

end module
