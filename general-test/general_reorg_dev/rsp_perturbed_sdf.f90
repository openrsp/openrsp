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

  implicit none

  public rsp_fds
  public get_fds
  public rsp_fock_lowerorder
  public get_fock_lowerorder

  contains

  recursive subroutine rsp_fds(mol, pert, kn, F, D, S)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    type(p_tuple), dimension(pert%n_perturbations) :: psub
    integer, dimension(2) :: kn
    integer :: i, j, k
    type(SDF) :: F, D, S

    ! Unless at final recursion level, recurse further
    ! Make all size (n - 1) subsets of the perturbations and recurse
    ! Then (at final recursion level) get perturbed F, D, S 

    if (pert%n_perturbations > 1) then

       call make_p_tuple_subset(pert, psub)

       do i = 1, size(psub)

          call rsp_fds(mol, psub(i), kn, F, D, S)

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

          call get_fds(mol, p_tuple_standardorder(pert), F, D, S)

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
  subroutine get_fds(mol, pert, F, D, S)

    use interface_rsp_solver, only: rsp_mosolver_exec

    implicit none

    type(rsp_cfg) :: mol
    integer :: sstr_incr, i, j, superstructure_size
    integer, allocatable, dimension(:) :: ind
    integer, dimension(0) :: noc
    character, dimension(0) :: nof
    type(p_tuple) :: pert
    type(p_tuple), allocatable, dimension(:,:) :: derivative_structure
    type(SDF) :: F, D, S
    type(matrix) :: X(1), RHS(1), A, B
    type(matrix), dimension(product(pert%pdim)) :: Fp, Dp, Sp, Dh

    A = mat_alloc_like(mol%zeromat)
    A = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(A)

    B = mat_alloc_like(mol%zeromat)
    B = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(B)

    call sdf_getdata_s(D, get_emptypert(), (/1/), A)
    call sdf_getdata_s(S, get_emptypert(), (/1/), B)

    ! Get the appropriate Fock/density/overlap matrices

    ! 1. Call ovlint and store perturbed overlap matrix

    call rsp_ovlint(mol, pert%n_perturbations, pert%plab, &
                    (/ (1, j = 1, pert%n_perturbations) /), pert%pdim, Sp)
    call sdf_add(S, pert, Sp)

    ! INITIALIZE AND STORE D INSTANCE WITH ZEROES
    ! THE ZEROES WILL ENSURE THAT TERMS INVOLVING THE HIGHEST ORDER DENSITY MATRICES
    ! WILL BE ZERO IN THE CONSTRUCTION OF Dp

    do i = 1, product(pert%pdim)

       Dp(i) = mat_alloc_like(mol%zeromat)
       Dp(i) = mat_zero_like(mol%zeromat)
       call mat_ensure_alloc(Dp(i))

       Dh(i) = mat_alloc_like(mol%zeromat)
       Dh(i) = mat_zero_like(mol%zeromat)
       call mat_ensure_alloc(Dh(i))

       Fp(i) = mat_alloc_like(mol%zeromat)
       Fp(i) = mat_zero_like(mol%zeromat)
       call mat_ensure_alloc(Fp(i))

    end do

    call sdf_add(D, pert, Dp)

    ! 2. Construct Dp and the initial part of Fp
    ! a) For the initial part of Fp: Make the initial recursive (lower order) 
    ! oneint, twoint, and xcint calls as needed

    call fock_lowerorder_cache_allocate(fock_lowerorder_cache)

    call rsp_fock_lowerorder(mol, pert, pert%n_perturbations, 1, (/get_emptypert()/), &
                         0, D, product(pert%pdim), Fp, fock_lowerorder_cache)

    deallocate(fock_lowerorder_cache)    

    call sdf_add(F, pert, Fp)

    ! b) For Dp: Create differentiation superstructure: First dryrun for size, and
    ! then the actual superstructure call

    superstructure_size = derivative_superstructure_getsize(mol, pert, &
                          (/pert%n_perturbations, pert%n_perturbations/), .FALSE., &
                          (/get_emptypert(), get_emptypert(), get_emptypert()/))

    sstr_incr = 0

    allocate(derivative_structure(superstructure_size,3))
    allocate(ind(pert%n_perturbations))

    call derivative_superstructure(mol, pert, (/pert%n_perturbations, &
         pert%n_perturbations/), .FALSE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         superstructure_size, sstr_incr, derivative_structure)

    call make_triangulated_indices(..., indices)

    do i = 1, size(indices, 1)

       ind = indices(i, :)

       Dp(i) = 1.0d0 * rsp_get_matrix_z(mol, superstructure_size, derivative_structure, &
               (/pert%n_perturbations,pert%n_perturbations/), pert%n_perturbations, &
               (/ (j, j = 1, pert%n_perturbations) /), pert%n_perturbations, &
               ind, F, D, S)

       Dp(i) = Dp(i) - A * B * Dp(i) - Dp(i) * B * A

       call sdf_add(D, pert, Dp)

       ! 3. Complete the particular contribution to Fp
       ! NOTE (MaR): THERE SHOULD BE A CALL TO rsp_xcint HERE TOO
       ! THE IF CRITERION HERE NEEDS ANOTHER LOOK

!        if (pert%n_perturbations <=2) then
! 
          call rsp_twoint(mol, 0, nof, noc, pert%pdim, Dp(i), Fp(i:i))
! 
!        end if


       call sdf_add(F, pert, Fp)

       ! 4. Make right-hand side using Dp

       RHS(1) = mat_alloc_like(mol%zeromat)
       RHS(1) = mat_zero_like(mol%zeromat)
       call mat_ensure_alloc(RHS(1))

       RHS(1) = rsp_get_matrix_y(mol, superstructure_size, derivative_structure, &
                pert%n_perturbations, (/ (j, j = 1, pert%n_perturbations) /), &
                pert%n_perturbations, ind, F, D, S)
     
       X(1) = mat_alloc_like(mol%zeromat)
       X(1) = mat_zero_like(mol%zeromat)
       call mat_ensure_alloc(X(1))

       ! Note (MaR): What does the second argument in rsp_mosolver_exec mean?
       call rsp_mosolver_exec(RHS(1), (/0d0/), X)
       ! Note (MaR): Why multiply by -2 like below?
       X(1) = -2d0*X(1)
       RHS(1) = 0

       ! 5. Get Dh using the rsp equation solution X

       Dh(i) = X(1) * B * A - A * B * X(1)

       ! 6. Make homogeneous contribution to Fock matrix

       ! THE IF CRITERION HERE NEEDS ANOTHER LOOK

!        if (pert%n_perturbations <=2) then

          call rsp_twoint(mol, 0, nof, noc, pert%pdim, Dh(i), Fp(i:i))

!        end if


       ! 'NOTE (MaR): XCINT CALL SKIPPED FOR NOW'

       ! call rsp_xcint(mol, 0, nof, noc, pert%pdim, 2, &
       ! (/sdf_getdata(D, get_emptypert(), (/1/)), Dh(i)/), Fp(i:i))

       ! 7. Complete perturbed D with homogeneous part

       Dp(i) = Dp(i) + Dh(i)


    end do

    deallocate(derivative_structure)
    deallocate(ind)

    ! Add the final values to cache

    call sdf_add(F, pert, Fp)
    call sdf_add(D, pert, Dp)

  end subroutine


  recursive subroutine rsp_fock_lowerorder(mol, pert, total_num_perturbations, &
                       num_p_tuples, p_tuples, density_order, D, property_size, Fp, &
                       fock_lowerorder_cache)

    implicit none

    logical :: density_order_skip
    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    integer :: num_p_tuples, density_order, i, j, total_num_perturbations, property_size
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, t_new
    type(SDF) :: D
    type(property_cache) :: energy_cache
    type(matrix), dimension(property_size) :: Fp
    type(f_l_cache) :: fock_lowerorder_cache


    if (pert%n_perturbations >= 1) then

       ! The differentiation can do three things:
       ! 1. Differentiate the expression 'directly'

       if (p_tuples(1)%n_perturbations == 0) then

          call rsp_fock_lowerorder(mol, p_tuple_remove_first(pert), & 
               total_num_perturbations, num_p_tuples, &
               (/p_tuple_getone(pert,1), p_tuples(2:size(p_tuples))/), &
               density_order, D, property_size, Fp, fock_lowerorder_cache)

       else

          call rsp_fock_lowerorder(mol, p_tuple_remove_first(pert), &
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

          call rsp_fock_lowerorder(mol, p_tuple_remove_first(pert), &
               total_num_perturbations, num_p_tuples, &
               t_new, density_order + 1, D, property_size, Fp, fock_lowerorder_cache)

       end do

       ! 3. Chain rule differentiate w.r.t. the density (giving 
       ! a(nother) pert D contraction)

       call rsp_fock_lowerorder(mol, p_tuple_remove_first(pert), &
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

          call get_fock_lowerorder(mol, num_p_tuples, total_num_perturbations, &
                                   p_tuples, density_order, D, property_size, Fp, &
                                   fock_lowerorder_cache)

          write(*,*) 'Calculated perturbed Fock matrix lower order contribution'
          write(*,*) ' '

       else

          write(*,*) 'Skipping contribution: At least one contraction D perturbed' 
          write(*,*) 'at order for which perturbed D is to be found '
          write(*,*) ' '

       end if


    end if

  end subroutine




  subroutine get_fock_lowerorder(mol, num_p_tuples, total_num_perturbations, p_tuples, &
                                 density_order, D, property_size, Fp, &
                                 fock_lowerorder_cache)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    type(SDF) :: D
    type(property_cache) :: energy_cache
    type(matrix), allocatable, dimension(:) :: dens_tuple
    integer :: i, j, k, m, num_p_tuples, total_num_perturbations, &
               density_order, property_size, offset
    integer, dimension(0) :: noc
    integer, dimension(total_num_perturbations) :: ncarray, ncouter, ncinner, pidouter
    integer, allocatable, dimension(:) :: o_whichpert, o_whichpertbig, o_wh_forave
    integer, allocatable, dimension(:) :: inner_offsets, ncoutersmall, &
                                          pidoutersmall, ncinnersmall
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    type(matrix), allocatable, dimension(:) :: tmp
    type(matrix), dimension(property_size) :: Fp
    type(f_l_cache) :: fock_lowerorder_cache
    
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
    allocate(tmp(product(ncinnersmall)))

    o_whichpert = make_outerwhichpert(total_num_perturbations, num_p_tuples, p_tuples)

    call sortdimbypid(total_num_perturbations, total_num_perturbations - &
                      p_tuples(1)%n_perturbations, pidoutersmall, &
                      ncarray, ncoutersmall, o_whichpert)

    if (total_num_perturbations > p_tuples(1)%n_perturbations) then

       do i = 1, size(o_whichpert)

          if (.NOT.(o_whichpert(i) == 0)) then

             o_wh_forave(o_whichpert(i)) = i

          end if
  
       end do


       do i = 1, num_p_tuples

          dens_tuple(i) = mat_alloc_like(mol%zeromat)
          dens_tuple(i) = mat_zero_like(mol%zeromat)
          call mat_ensure_alloc(dens_tuple(i))

       end do

       do j = 1, size(tmp)

          tmp(j) = mat_alloc_like(mol%zeromat)
          tmp(j) = mat_zero_like(mol%zeromat)
          call mat_ensure_alloc(tmp(j))

       end do

       outer_indices_size = get_triang_tuples_size(outer tuples)
       inner_indices_size = get_triang_tuples_size(inner tuples)
    
       allocate(outer_indices(outer_indices_size,size(ncoutersmall)))
       allocate(inner_indices(inner_indices_size,size(ncinnersmall)))

       call make_triangulated_tuples_indices(num_p_tuples - 1, &
            p_tuples(2:num_p_tuples), outer_indices_size, outer_indices)

       if (p_tuples(1)%n_perturbations > 0) then

          call make_triangulated_tuples_indices(1, p_tuples(1), inner_indices_size, &
          outer_indices)

       end if

       allocate(inner_offsets(product(ncinner)))

       do i = 1, size(outer_indices, 1)

          do j = 2, num_p_tuples

             call sdf_getdata_s(D, p_tuples(j), (/ &
                             (outer_indices(i,o_wh_forave(p_tuples(j)%pid(k))), &
                             k = 1, p_tuples(j)%n_perturbations) /), dens_tuple(j))

          end do

! NOTE: IS THIS ZEROING OF tmp REDUNDANT?

          do j = 1, size(tmp)

             tmp(j) = mat_alloc_like(mol%zeromat)
             tmp(j) = mat_zero_like(mol%zeromat)
             call mat_ensure_alloc(tmp(j))

          end do


          if (num_p_tuples <= 1) then

             call rsp_oneint(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                             (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                             p_tuples(1)%pdim, tmp)

          end if


          if (num_p_tuples <= 2) then

             call rsp_twoint(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                             (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                             p_tuples(1)%pdim, dens_tuple(2), tmp)

          end if

!              write(*,*) 'Calling xcint (NOTE: XCINT CALL SKIPPED FOR NOW)'
!
!              call rsp_xcint(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
!                            (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
!                            p_tuples(1)%pdim, density_order, &
!                            (/ sdf_getdata(D, get_emptypert(), (/1/)), &
!                            (dens_tuple(k), k = 2, num_p_tuples) /), &
!                            tmp)


! Make tmp the correct ("inner-triangulated") size
! Make some storage array the size of the full amount of data to be cached
! Make sure that the data from the integral is put into the 
! appropriate positions in tmp
! After the accumulation loop: Add the appropriate elements of tmp to
! Fp and then add tmp to fock_lowerorder_cache



if (p_tuples(1)%n_perturbations > 0) then

do j = 1, size(inner_indices,1)

offset = get_one_triangular_offset((/inner_indices(j, :), outer_indices(i, :) /), &
         blk_info, ncarray)

lower_order_contribution(offset) = tmp(j)

end do

else

offset = get_one_triangular_offset((/outer_indices(i, :) /), &
         blk_info, ncarray)


lower_order_contribution(offset) = tmp(1)

end if







! This existing handling is commented out and is 
! to be replaced according to the plan above
!  
!           if (p_tuples(1)%n_perturbations > 0) then
! 
!              do j = 1, size(inner_indices, 1)
! 
!                 offset = get_one_tensor_offset( sum( (/ (p_tuples(k)%n_perturbations, & 
!                          k=1, num_p_tuples ) /) ), &
!                          (/ inner_indices(j,:), outer_indices(i,:) /), &
!                          (/ (p_tuples(k)%pid, k=1, num_p_tuples ) /), ncarray)
! 
!                 Fp(offset) = Fp(offset) + tmp(j)
! 
!              end do
! 
!           else
! 
!              offset = get_one_tensor_offset( sum( (/ (p_tuples(k)%n_perturbations, &
!                       k=2, num_p_tuples ) /) ), &
!                       (/ outer_indices(i,:) /), &
!                       (/ (p_tuples(k)%pid, k=2, num_p_tuples ) /), ncarray)
! 
!              Fp(offset) = Fp(offset) + tmp(1)
! 
!           end if

       end do

! get triangulated indices for Fp

do i = 1, size(triang_indices_fp, 1)

fp_offset = get_one_triangular_offset(unsplit p_tuples indices)
lo_offset = get_one_triangular_offset(split p_tuples indices)

Fp(fp_offset) = Fp(fp_offset) + lower_order_contribution(lo_offset)

end do

fock_lowerorder_cache_add(lower_order_contribution)


    deallocate(inner_offsets)

    else

       if (num_p_tuples <= 1) then

          call rsp_oneint(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                          (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                          p_tuples(1)%pdim, &
                          Fp)
! Waiting for developments in rsp_contribs
! NOTE: Find out if necessary ovlint/oneint in "outer indices case" above
! NOTE: Add (a) corresponding call(s) for the energy term contributions (see e.g. eqn.
! 209 in ajt_rsp)
!           call rsp_ovlint(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
!                           (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
!                           p_tuples(1)%pdim, w = p_tuples(1)%freq, fock = Fp)

       end if

       if (num_p_tuples <= 2) then

          call rsp_twoint(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
               (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
               p_tuples(1)%pdim, &
               sdf_getdata(D, get_emptypert(), (/1/) ), &
               Fp)

       end if

!        write(*,*) 'Calling xcint (NOTE: XCINT CALL SKIPPED FOR NOW)'

!        call rsp_xcint(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
!                       (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
!                       (/product(p_tuples(1)%pdim)/), 1, &
!                       (/ sdf_getdata(D, get_emptypert(), (/1/)) /), &
!                       Fp)


       ! NOTE: THERE IS NO NEED TO CACHE THE "ALL INNER" CONTRIBUTION
       ! It should be possible to just add it to Fp like already done above
       ! even with the extra complexity from the triangularization 

    end if


    deallocate(ncoutersmall)
    deallocate(ncinnersmall)
    deallocate(pidoutersmall)
    deallocate(o_whichpert)
    deallocate(o_wh_forave)
    deallocate(dens_tuple)
    deallocate(outer_indices)
    deallocate(inner_indices)
    deallocate(tmp)

  end subroutine

end module
