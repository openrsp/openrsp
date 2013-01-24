! Copyright 2012 Magnus Ringholm
!
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

!> @file Contains module rsp_general

!> General response routines. This module organizes, computes and prints
!> response function tensors.
module rsp_general

  use matrix_defop
  use matrix_lowlevel, only: mat_init
  use rsp_contribs
  use rsp_field_tuple
  use rsp_indices_and_addressing
  use rsp_perturbed_matrices
  use rsp_perturbed_sdf
  use rsp_property_caching
  use rsp_sdf_caching
  use interface_xc

  implicit none

public rsp_cfg
  public rsp_prop
  public get_prop
  public rsp_energy
  public get_energy
  public rsp_pulay_kn
  public get_pulay_kn
  public rsp_pulay_lag
  public get_pulay_lag
  public rsp_idem_lag
  public get_idem_lag
  public rsp_scfe_lag
  public get_scfe_lag
  public print_rsp_tensor
  public print_rsp_tensor_stdout
  public print_rsp_tensor_stdout_tr

  type(matrix) :: zeromat

  contains

! Supports supplying F, D, S instances through F_already etc.
! (must then also give zeromat_already) or 
! initializing new ones through F_unpert etc.
! F, D, S instances that are supplied are used blindly, so 
! the user must ensure that the instances contain the correct information w.r.t.
! e.g. the molecular geometry and basis set

  subroutine rsp_prop(pert_unordered, kn, F_unpert, D_unpert, S_unpert, &
                      F_already, D_already, S_already, zeromat_already, file_id)

    implicit none

    character(len=*), optional :: file_id
    type(p_tuple) :: pert, pert_unordered
    type(matrix), optional :: F_unpert, D_unpert, S_unpert, zeromat_already
    type(SDF), optional :: F_already, D_already, S_already
    type(SDF), pointer :: F, D, S
    integer, dimension(2) :: kn
    real :: timing_start, timing_end
    integer :: i, j, num_blks, property_size, nr_ao
    integer, allocatable, dimension(:) :: blk_sizes
    integer, allocatable, dimension(:,:) :: blk_info
    complex(8), allocatable, dimension(:) :: prop

    open(unit=257, file='totterms', status='replace', action='write') 
    write(257,*) 'START'
    close(257)

    open(unit=257, file='cachehit', status='replace', action='write') 
    write(257,*) 'START'
    close(257)

    pert = p_tuple_standardorder(pert_unordered)

! MaR: The get_bestkn function is taken out of use for now until it has been improved
! The choice of k, n will likely be made before this routine is called anyway
! The call (commented out) is kept here for now
!    kn = get_bestkn(pert)
    write(*,*) ' '
    write(*,*) 'General OpenRSP code called'
    write(*,*) ' '
    write(*,*) 'The field labels are', pert%plab
    write(*,*) 'Frequencies are', pert%freq
    write(*,*) 'Dimensions are', pert%pdim
    write(*,*) ' '
    write(*,*) 'Choice of k, n is ', kn(1), ' and ', kn(2)
    write(*,*) ' '

    if (present(S_already) .eqv. .FALSE.) then

! Assumes S_unpert, D_unpert, F_unpert is present

       ! ASSUME CLOSED SHELL
       call mat_init(zeromat, S_unpert%nrow, S_unpert%ncol, &
                     .true., .false., .false., .false., .false.)
       call mat_init_like_and_zero(S_unpert, zeromat)


       call sdf_setup_datatype(S, S_unpert)
       call sdf_setup_datatype(D, D_unpert)
       call sdf_setup_datatype(F, F_unpert)

    else

! Assumes zeromat_already is present

       ! ASSUME CLOSED SHELL
       call mat_init(zeromat, zeromat_already%nrow, zeromat_already%ncol, &
                     .true., .false., .false., .false., .false.)
       call mat_init_like_and_zero(zeromat_already, zeromat)

    end if

    num_blks = get_num_blks(pert)
    allocate(blk_info(num_blks, 3))
    allocate(blk_sizes(num_blks))
    blk_info = get_blk_info(num_blks, pert)
    blk_sizes = get_triangular_sizes(num_blks, blk_info(1:num_blks, 2), &
                                     blk_info(1:num_blks, 3))

    property_size = get_triangulated_size(num_blks, blk_info)

    nr_ao = zeromat%nrow

    allocate(prop(property_size))
    prop = 0.0

    write(*,*) 'Starting clock: About to call get_prop routine'
    write(*,*) ' '
    call cpu_time(timing_start)

    if (present(F_already)) then

       call get_prop(pert, kn, nr_ao, num_blks, blk_sizes, blk_info, &
                     property_size, prop, F_already, D_already, S_already)

    else

       call get_prop(pert, kn, nr_ao, num_blks, blk_sizes, blk_info, &
                     property_size, prop, F, D, S)

   end if


    call cpu_time(timing_end)
    write(*,*) 'Clock stopped: Property was calculated'
    write(*,*) 'Time spent in get_prop:',  timing_end - timing_start, ' seconds'
    write(*,*) ' '

    write(*,*) 'Property was calculated'
    write(*,*) ' '

    if (present(file_id)) then
       open(unit=260, file='rsp_tensor_' // trim(adjustl(file_id)), &
            status='replace', action='write') 
    else
       open(unit=260, file='rsp_tensor', &
            status='replace', action='write') 
    end if

    write(260,*) ' '


    call print_rsp_tensor_tr(1, pert%n_perturbations, pert%pdim, &
    (/ (1, j = 1, (pert%n_perturbations - 1) ) /), num_blks, blk_sizes, &
    blk_info, property_size, prop, 260)

    close(260)

    if (present(file_id)) then
       write(*,*) 'Property was printed to rsp_tensor_' // trim(adjustl(file_id))
    else
       write(*,*) 'Property was printed to rsp_tensor'
    end if

    write(*,*) ' '
    write(*,*) 'End of print'

    open(unit=257, file='totterms', status='old', action='write', position='append') 
    write(257,*) 'END'
    close(257)

    open(unit=257, file='cachehit', status='old', action='write', position='append') 
    write(257,*) 'END'
    close(257)

    deallocate(blk_info)

  end subroutine

  subroutine get_prop(pert, kn, nr_ao, num_blks, blk_sizes, blk_info, &
                      property_size, prop, F, D, S)

    implicit none

    type(SDF) :: F, D, S
    type(p_tuple) :: pert, emptypert
    type(p_tuple), dimension(2) :: emptyp_tuples
    integer :: property_size, nr_ao, num_blks
    integer, dimension(2) :: kn
    integer, dimension(num_blks) :: blk_sizes
    integer, dimension(num_blks,3) :: blk_info
    complex(8), dimension(property_size) :: prop, p_diff
    type(property_cache), pointer :: energy_cache, pulay_kn_cache, &
                                     pulay_lag_cache, idem_cache, scfe_cache

    emptypert%n_perturbations = 0
    allocate(emptypert%pdim(0))    
    allocate(emptypert%plab(0))
    allocate(emptypert%pid(0))
    allocate(emptypert%freq(0))

    emptyp_tuples = (/emptypert, emptypert/)

    ! Get all necessary F, D, S derivatives as dictated by
    ! number of perturbations and kn

    prop = 0.0
!     p_diff = prop

    call rsp_fds(zeromat, pert, kn, F, D, S)
 
    write(*,*) ' '
    write(*,*) 'Finished calculation of perturbed overlap/density/Fock matrices'
    write(*,*) ' '

    call property_cache_allocate(energy_cache)
    call rsp_energy(pert, pert%n_perturbations, kn, 1, (/emptypert/), 0, D, &
                  property_size, energy_cache, prop)

    write(*,*) ' '
    write(*,*) 'Finished calculating HF energy-type contributions'
    write(*,*) ' '

    deallocate(energy_cache)

! write(*,*) 'prop incr'
! write(*,*) real(prop - p_diff)
! 
!     p_diff = prop

    write(*,*) ' '
    write(*,*) 'Calculating exchange/correlation contributions'
    write(*,*) ' '

    call rsp_xcave(nr_ao, pert, kn, num_blks, blk_sizes, blk_info, property_size, prop, D)

    write(*,*) ' '
    write(*,*) 'Finished calculating exchange/correlation contributions'
    write(*,*) ' '

! write(*,*) 'prop incr'
! write(*,*) real(prop - p_diff)
! 
!     p_diff = prop

    call property_cache_allocate(pulay_kn_cache)
    call rsp_pulay_kn(pert, kn, (/emptypert, emptypert/), S, D, F, &
                      property_size, pulay_kn_cache, prop)

    write(*,*) ' '
    write(*,*) 'Finished calculating Pulay k-n type contributions'
    write(*,*) ' '

    deallocate(pulay_kn_cache)

! write(*,*) 'prop incr'
! write(*,*) real(prop - p_diff)
! 
!     p_diff = prop

    call property_cache_allocate(pulay_lag_cache)
    call rsp_pulay_lag(p_tuple_remove_first(pert), kn, &
                       (/p_tuple_getone(pert,1), emptypert/), &
                       S, D, F, property_size, pulay_lag_cache, prop)

    write(*,*) ' '
    write(*,*) 'Finished calculating Pulay lagrangian type contributions' 
    write(*,*) ' '

    deallocate(pulay_lag_cache)

! write(*,*) 'prop incr'
! write(*,*) real(prop - p_diff)
! 
!     p_diff = prop

    call property_cache_allocate(idem_cache)
    call rsp_idem_lag(p_tuple_remove_first(pert), kn, &
                      (/p_tuple_getone(pert,1), emptypert/), &
                      S, D, F, property_size, idem_cache, prop)

    write(*,*) ' '
    write(*,*) 'Finished calculating idempotency lagrangian type contributions'
    write(*,*) ' '

    deallocate(idem_cache)

! write(*,*) 'prop incr'
! write(*,*) real(prop - p_diff)
! 
!     p_diff = prop

    call property_cache_allocate(scfe_cache)
    call rsp_scfe_lag(p_tuple_remove_first(pert), kn, &
                      (/p_tuple_getone(pert,1), emptypert/), &
                      S, D, F, property_size, scfe_cache, prop)

    write(*,*) ' '
    write(*,*) 'Finished calculating SCF lagrangian type contributions'
    write(*,*) ' '

    deallocate(scfe_cache)

! write(*,*) 'prop incr'
! write(*,*) real(prop - p_diff)
! 
!     p_diff = prop

  end subroutine


  ! Calculate and add all the energy contributions

  recursive subroutine rsp_energy(pert, total_num_perturbations, kn, num_p_tuples, &
                                p_tuples, density_order, D, property_size, cache, prop)

    implicit none

    logical :: e_knskip
    type(p_tuple) :: pert
    integer, dimension(2) :: kn
    integer :: num_p_tuples, density_order, i, j, total_num_perturbations, property_size
    type(p_tuple), dimension(num_p_tuples) :: p_tuples, t_new
    type(SDF) :: D
    type(property_cache) :: cache
    complex(8), dimension(property_size) :: prop

    if (pert%n_perturbations >= 1) then

       ! The differentiation can do three things:
       ! 1. Differentiate the energy expression 'directly'

    if (p_tuples(1)%n_perturbations == 0) then

       call rsp_energy(p_tuple_remove_first(pert), total_num_perturbations, &
       kn, num_p_tuples, (/p_tuple_getone(pert,1), p_tuples(2:size(p_tuples))/), &
       density_order, D, property_size, cache, prop)

    else


       call rsp_energy(p_tuple_remove_first(pert), total_num_perturbations,  &
       kn, num_p_tuples, (/p_tuple_extend(p_tuples(1), p_tuple_getone(pert,1)), &
       p_tuples(2:size(p_tuples))/), density_order, D, property_size, cache, prop)

    end if
    
       ! 2. Differentiate all of the contraction densities in turn

       ! Find the number of terms

       do i = 2, num_p_tuples

          t_new = p_tuples

          if (p_tuples(i)%n_perturbations == 0) then

             t_new(i) = p_tuple_getone(pert, 1)

          else

             t_new(i) = p_tuple_extend(t_new(i), p_tuple_getone(pert, 1))

          end if

          call rsp_energy(p_tuple_remove_first(pert), total_num_perturbations, &
          kn, num_p_tuples, t_new, density_order + 1, D, property_size, cache, prop)

       end do

       ! MaR: Since we are only calculating Hartree-Fock type energy terms here,
       ! we don't need to go beyond to perturbed contraction density matrices
       ! (but that is in general needed for XC contributions)
       if (num_p_tuples < 3) then

          ! 3. Chain rule differentiate the energy w.r.t. the density (giving 
          ! a(nother) pert D contraction)

          call rsp_energy(p_tuple_remove_first(pert), total_num_perturbations, &
          kn, num_p_tuples + 1, (/p_tuples(:), p_tuple_getone(pert, 1)/), &
          density_order + 1, D, property_size, cache, prop)

       end if

    ! At the final recursion level: Calculate the contribution (if k,n choice of rule
    ! allows it) or get it from cache if it was already calculated (and if k,n choice 
    ! of rule allows it)

    else

       e_knskip = .FALSE.


!        p_tuples = p_tuples_standardorder(num_p_tuples, p_tuples)

!        write(*,*) 'Getting energy contribution'

       do i = 1, num_p_tuples
 
          if (i > 1) then

!              write(*,*) 'D ', p_tuples(i)%pid

             if(kn_skip(p_tuples(i)%n_perturbations, p_tuples(i)%pid, kn) .EQV. .TRUE.) then

                e_knskip = .TRUE.

             end if
          
          elseif (i == 1) then

!              write(*,*) 'E ', p_tuples(i)%pid

          end if

       end do


       if (e_knskip .EQV. .FALSE.) then

          open(unit=257, file='totterms', status='old', action='write', &
               position='append') 
          write(257,*) 'T'
          close(257)
          
!           write(*,*) 'Evaluating property_cache_already'

          if (property_cache_already(cache, num_p_tuples, &
              p_tuples_standardorder(num_p_tuples, p_tuples)) .EQV. .TRUE.) then

             open(unit=257, file='cachehit', status='old', action='write', &
                  position='append') 
             write(257,*) 'T'
             close(257)

!              write(*,*) 'Getting values from cache'

             ! NOTE (MaR): EVERYTHING MUST BE STANDARD ORDER IN 
             ! THIS CALL (LIKE property_cache_getdata ASSUMES)
             call property_cache_getdata(cache, num_p_tuples, &
                  p_tuples_standardorder(num_p_tuples, p_tuples), property_size, prop)

!              write(*,*) ' '
       
          else


       write(*,*) 'Calculating energy contribution'

       do i = 1, num_p_tuples
 
          if (i > 1) then

             write(*,*) 'D ', p_tuples(i)%pid

!              if(kn_skip(p_tuples(i)%n_perturbations, p_tuples(i)%pid, kn) .EQV. .TRUE.) then
! 
!                 e_knskip = .TRUE.
! 
!              end if
          
          elseif (i == 1) then

             write(*,*) 'E ', p_tuples(i)%pid

          end if

       end do


             call get_energy(num_p_tuples, total_num_perturbations, & 
                  p_tuples_standardorder(num_p_tuples, p_tuples), &
                  density_order, D, property_size, cache, prop)

                  write(*,*) 'Calculated energy contribution'
                  write(*,*) ' '

          end if

       else

!           write(*,*) 'Energy contribution was k-n skipped'
!           write(*,*) ' '

       end if

    end if

  end subroutine

  subroutine get_energy(num_p_tuples, total_num_perturbations, &
                        p_tuples, density_order, D, property_size, cache, prop)

    implicit none

    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    type(p_tuple) :: merged_p_tuple, t_matrix_bra, t_matrix_ket
    type(SDF) :: D
    type(property_cache) :: cache
    type(matrix), allocatable, dimension(:) :: dens_tuple
    type(rsp_field), allocatable, dimension(:) :: nucpot_pert
    integer :: i, j, k, m, n, num_p_tuples, total_num_perturbations, density_order, &
             property_size, offset, dtup_ind, pr_offset, ec_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, merged_nblks
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:) :: ncinnersmall, blk_sizes_merged
    integer, allocatable, dimension(:,:) :: triang_indices_pr, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    integer, dimension(total_num_perturbations) :: ncarray, ncouter, ncinner, pidouter, &
                                              pids_current_contribution, translated_index
    integer, allocatable, dimension(:) :: o_whichpert, o_whichpertbig, o_wh_forave
    integer, allocatable, dimension(:) :: inner_offsets, ncoutersmall, pidoutersmall
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    complex(8), allocatable, dimension(:) :: tmp, contrib, prop_forcache
    complex(8), dimension(property_size) :: prop

!    ncarray = get_ncarray(total_num_perturbations, num_p_tuples, p_tuples)
!    ncouter = nc_only(total_num_perturbations, total_num_perturbations - &
!              p_tuples(1)%n_perturbations, num_p_tuples - 1, &
!              p_tuples(2:num_p_tuples), ncarray)
!    ncinner = nc_only(total_num_perturbations, p_tuples(1)%n_perturbations, 1, &
!                      p_tuples(1), ncarray)

    allocate(dens_tuple(num_p_tuples))
    allocate(nucpot_pert(p_tuples(1)%n_perturbations))
 !   allocate(ncoutersmall(total_num_perturbations - p_tuples(1)%n_perturbations))
 !   allocate(ncinnersmall(p_tuples(1)%n_perturbations))
 !   allocate(pidoutersmall(total_num_perturbations - p_tuples(1)%n_perturbations))

 !   ncoutersmall = nc_onlysmall(total_num_perturbations, total_num_perturbations - &
 !                  p_tuples(1)%n_perturbations, num_p_tuples - 1, &
 !                  p_tuples(2:num_p_tuples), ncarray)
 !   ncinnersmall = nc_onlysmall(total_num_perturbations, p_tuples(1)%n_perturbations, &
 !                  1, p_tuples(1), ncarray)
 !   pidoutersmall = get_pidoutersmall(total_num_perturbations - &
 !                   p_tuples(1)%n_perturbations, num_p_tuples - 1, &
 !                   p_tuples(2:num_p_tuples))

    allocate(o_whichpert(total_num_perturbations))
    allocate(o_wh_forave(total_num_perturbations))
    allocate(nfields(num_p_tuples))
    allocate(nblks_tuple(num_p_tuples))

    do i = 1, num_p_tuples

       nfields(i) = p_tuples(i)%n_perturbations
       nblks_tuple(i) = get_num_blks(p_tuples(i))

    end do

    allocate(blks_tuple_info(num_p_tuples, total_num_perturbations, 3))
    allocate(blks_tuple_triang_size(num_p_tuples))
    allocate(blk_sizes(num_p_tuples, total_num_perturbations))
    allocate(blk_sizes_merged(total_num_perturbations))

    do i = 1, num_p_tuples

       blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p_tuples(i))
       blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                                   blks_tuple_info(i, 1:nblks_tuple(i), :))

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
    allocate(contrib(inner_indices_size))
    allocate(prop_forcache(inner_indices_size*outer_indices_size))

    prop_forcache = 0.0
    contrib = 0.0

!    call sortdimbypid(total_num_perturbations, total_num_perturbations - &
!                      p_tuples(1)%n_perturbations, pidoutersmall, &
!                      ncarray, ncoutersmall, o_whichpert)

    if (total_num_perturbations > p_tuples(1)%n_perturbations) then

       allocate(outer_indices(outer_indices_size,total_num_perturbations - &
                p_tuples(1)%n_perturbations))
       allocate(inner_indices(inner_indices_size,p_tuples(1)%n_perturbations))

       k = 1
    
       do i = 2, num_p_tuples
          do j = 1, p_tuples(i)%n_perturbations
    
             o_wh_forave(p_tuples(i)%pid(j)) = k
             k = k + 1
    
          end do
       end do
    
       k = 1
   
!       do i = 2, num_p_tuples
!          do j = 1, p_tuples(i)%n_perturbations
!   
!             ncoutersmall(k) =  p_tuples(i)%pdim(j)
!             k = k + 1
!   
!          end do
!       end do
   
       do i = 1, num_p_tuples
   
          ! ASSUME CLOSED SHELL
          call mat_init(dens_tuple(i), zeromat%nrow, zeromat%ncol, &
                        .false., .false., .false., .false., .false.)
          call mat_init_like_and_zero(zeromat, dens_tuple(i))
   
       end do
   
       call make_triangulated_tuples_indices(num_p_tuples - 1, total_num_perturbations, & 
            nblks_tuple(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, &
            :, :), blks_tuple_triang_size(2:num_p_tuples), outer_indices)
   
   
       if (p_tuples(1)%n_perturbations > 0) then
   
          call make_triangulated_indices(nblks_tuple(1), blks_tuple_info(1, &
               1:nblks_tuple(1), :), blks_tuple_triang_size(1), inner_indices)
   
       end if
   
       do i = 1, size(outer_indices, 1)
   
          dtup_ind = 0
   
          do j = 2, num_p_tuples
   
             call sdf_getdata_s(D, p_tuples(j), outer_indices(i, &
                  dtup_ind+1:dtup_ind + p_tuples(j)%n_perturbations), dens_tuple(j))
   
             dtup_ind = dtup_ind + p_tuples(j)%n_perturbations
   
          end do
   
          tmp = 0.0
          contrib = 0.0
   
          if (num_p_tuples == 1) then
   
             call rsp_oneave_tr(p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                            (/ (1, j = 1, p_tuples(1)%n_perturbations) /), & 
                            p_tuples(1)%pdim, sdf_getdata(D, get_emptypert(), (/1/)), &
                            nblks_tuple(1),  blks_tuple_info(1, 1:nblks_tuple(1), :), &
                            blk_sizes(1, 1:nblks_tuple(1)), inner_indices_size, contrib)
   
          elseif (num_p_tuples == 2) then
   
             call rsp_oneave_tr(p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                            (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                            p_tuples(1)%pdim, dens_tuple(2), &
                            nblks_tuple(1),  blks_tuple_info(1, 1:nblks_tuple(1), :), &
                            blk_sizes(1, 1:nblks_tuple(1)), inner_indices_size, contrib)
   
          end if
   
          tmp = tmp + contrib
          contrib = 0.0
   
          if (num_p_tuples == 1) then
   
             t_matrix_bra = get_emptypert()
             t_matrix_ket = get_emptypert()

             call rsp_ovlave_t_matrix(p_tuples(1)%n_perturbations, p_tuples(1), &
                                      t_matrix_bra, t_matrix_ket, &
                                      sdf_getdata(D, get_emptypert(), (/1/)), &
                                      inner_indices_size, contrib)

   
          elseif (num_p_tuples == 2) then

             t_matrix_bra = get_emptypert()
             t_matrix_ket = get_emptypert()

             call rsp_ovlave_t_matrix(p_tuples(1)%n_perturbations, p_tuples(1), &
                                      t_matrix_bra, t_matrix_ket, &
                                      dens_tuple(2), inner_indices_size, contrib)
   
          end if
   
          tmp = tmp + contrib
          contrib = 0.0
   
          if (num_p_tuples == 1) then
   
             call rsp_twoave_tr(p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                             (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                             p_tuples(1)%pdim, sdf_getdata(D, get_emptypert(), &
                             (/1/)), sdf_getdata(D, get_emptypert(), (/1/)) , &
                             inner_indices_size, contrib)
   
          elseif (num_p_tuples == 2) then
   
             call rsp_twoave_tr(p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                             (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                             p_tuples(1)%pdim, dens_tuple(2), &
                             sdf_getdata(D, get_emptypert(), (/1/)) , &
                             inner_indices_size, contrib)
    
          elseif (num_p_tuples == 3) then
    
             call rsp_twoave_tr(p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                             (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                             p_tuples(1)%pdim, dens_tuple(2), dens_tuple(3), &
                             inner_indices_size, contrib)
    
          end if
    
          tmp = tmp + contrib
    
          if (p_tuples(1)%n_perturbations > 0) then
    
             do j = 1, size(inner_indices, 1)
    
                offset = get_triang_blks_tuple_offset(num_p_tuples, total_num_perturbations, &
                         nblks_tuple, (/ (p_tuples(k)%n_perturbations, k = 1, num_p_tuples) /), &
                         blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                         (/inner_indices(j, :), outer_indices(i, :) /)) 
    
                prop_forcache(offset) = prop_forcache(offset) + tmp(j)
    
             end do
    
          else
    
             offset = get_triang_blks_tuple_offset(num_p_tuples - 1, total_num_perturbations,  &
                      nblks_tuple(2:num_p_tuples), &
                      (/ (p_tuples(k)%n_perturbations, k = 2, num_p_tuples) /), &
                      blks_tuple_info(2:num_p_tuples, :, :), blk_sizes(2:num_p_tuples,:), & 
                      blks_tuple_triang_size(2:num_p_tuples), (/outer_indices(i, :) /)) 
    
             prop_forcache(offset) = prop_forcache(offset) + tmp(1)
    
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
   
             ! This can be problematic - consider rewriting merge_p_tuple as subroutine
             merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))
   
          end do
   
       end if
   
       merged_p_tuple = p_tuple_standardorder(merged_p_tuple)
   
       ! MR: NOT DOING THE FOLLOWING FOR THE MERGED PERT ASSUMES THAT 
       ! PIDS ARE IN STANDARD ORDER? FIND OUT
   
       k = 1
       do i = 1, num_p_tuples
          do j = 1, p_tuples(i)%n_perturbations
             pids_current_contribution(k) = p_tuples(i)%pid(j)
             k = k + 1
          end do
       end do
    
       merged_nblks = get_num_blks(merged_p_tuple)
       
       allocate(merged_blk_info(1, merged_nblks, 3))
       
       merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)
       blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
       merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
       merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)
    
       allocate(triang_indices_pr(merged_triang_size, sum(merged_blk_info(1, :,2))))
    
       call make_triangulated_indices(merged_nblks, merged_blk_info, & 
                                      merged_triang_size, triang_indices_pr)
    
       do i = 1, size(triang_indices_pr, 1)
    
          pr_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
                      (/sum(nfields)/), &
                      (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
                      (/triang_indices_pr(i, :) /))
    
          do j = 1, total_num_perturbations
    
             translated_index(j) = triang_indices_pr(i,pids_current_contribution(j))
    
          end do
    
          if (p_tuples(1)%n_perturbations > 0) then
    
             ec_offset = get_triang_blks_tuple_offset(num_p_tuples, &
                         total_num_perturbations, nblks_tuple, &
                         nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                         (/ translated_index(:) /))
    
          else
    
             ec_offset = get_triang_blks_tuple_offset(num_p_tuples - 1, &
                         total_num_perturbations, nblks_tuple(2:num_p_tuples), &
                         nfields(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, :, :), &
                         blk_sizes(2:num_p_tuples,:), blks_tuple_triang_size(2:num_p_tuples), &
                         (/ translated_index(:) /))
    
          end if
   
          prop(pr_offset) = prop(pr_offset) + prop_forcache(ec_offset)
   
       end do
   
       call property_cache_add_element(cache, num_p_tuples, p_tuples, &
                                       inner_indices_size * outer_indices_size, prop_forcache)   
   
       deallocate(triang_indices_pr)
       deallocate(outer_indices)
       deallocate(inner_indices)

    else

       ! MR: THIS IS THE CASE 'ALL INDICES ARE INNER INDICES'
       ! THIS TERM OCCURS ONLY ONCE AND DOES NOT NEED TO BE CACHED

       do i = 1, p_tuples(1)%n_perturbations

          nucpot_pert(i) = rsp_field(p_tuples(1)%plab(i), p_tuples(1)%freq(i), 1, &
                                     p_tuples(1)%pdim(i))

       end do

       tmp = 0.0
       contrib = 0.0

       call rsp_nucpot_tr(nucpot_pert, property_size, contrib) 

       tmp = tmp + contrib
       contrib = 0.0

       call rsp_oneave_tr(p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                       (/ (1, j = 1, p_tuples(1)%n_perturbations) /), p_tuples(1)%pdim, &
                       sdf_getdata(D, get_emptypert(), (/1/)) , &
                       nblks_tuple(1),  blks_tuple_info(1, 1:nblks_tuple(1), :), &
                       blk_sizes(1, 1:nblks_tuple(1)), property_size, contrib)

       tmp = tmp + contrib
       contrib = 0.0

       t_matrix_bra = get_emptypert()
       t_matrix_ket = get_emptypert()

       call rsp_ovlave_t_matrix(p_tuples(1)%n_perturbations, p_tuples(1), &
                                t_matrix_bra, t_matrix_ket, &
                                sdf_getdata(D, get_emptypert(), (/1/)), &
                                inner_indices_size, contrib)

       tmp = tmp + contrib
       contrib = 0.0

       call rsp_twoave_tr(p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                       (/ (1, j = 1, p_tuples(1)%n_perturbations) /), p_tuples(1)%pdim, &
                       sdf_getdata(D, get_emptypert(), (/1/)) , &
                       sdf_getdata(D, get_emptypert(), (/1/)) , property_size, contrib)

       tmp = tmp + 0.5*(contrib)

       prop =  prop + tmp

       call p_tuple_p1_cloneto_p2(p_tuples(1), merged_p_tuple)

       do i = 2, num_p_tuples

          ! This can be problematic - consider rewriting merge_p_tuple as subroutine
          merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))

       end do

       merged_p_tuple = p_tuple_standardorder(merged_p_tuple)
       merged_nblks = get_num_blks(merged_p_tuple)
       allocate(merged_blk_info(1, merged_nblks, 3))
       merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)
       blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
       merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))

    end if

!  write(*,*) 'energy contribution'
!  call print_rsp_tensor_stdout_tr(1, total_num_perturbations, merged_p_tuple%pdim, &
!  (/ (1, j = 1, (merged_p_tuple%n_perturbations - 1) ) /), merged_nblks, blk_sizes_merged, &
!  merged_blk_info, property_size, prop_forcache)

    deallocate(merged_blk_info)
    deallocate(nfields)
    deallocate(nblks_tuple)
    deallocate(blks_tuple_info)
    deallocate(blks_tuple_triang_size)
    deallocate(blk_sizes)
    deallocate(blk_sizes_merged)
    deallocate(nucpot_pert)
    deallocate(dens_tuple)
!    deallocate(ncoutersmall)
!    deallocate(ncinnersmall)
!    deallocate(pidoutersmall)
    deallocate(o_whichpert)
    deallocate(o_wh_forave)
    deallocate(tmp)
    deallocate(contrib)
    deallocate(prop_forcache)

  end subroutine


  recursive subroutine rsp_pulay_kn(pert, kn, p12, S, D, F, property_size, cache, prop)

    implicit none

    type(p_tuple) :: pert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    integer :: property_size, i
    integer, dimension(2) :: kn
    complex(8), dimension(property_size) :: prop
    
    if (pert%n_perturbations > 0) then

       call rsp_pulay_kn(p_tuple_remove_first(pert), kn, &
       (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), S, D, F, property_size, &
       cache, prop)

       call rsp_pulay_kn(p_tuple_remove_first(pert), kn, &
       (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), S, D, F, property_size, &
       cache, prop)

    else

       if (kn_skip(p12(2)%n_perturbations, p12(2)%pid, kn) .EQV. .FALSE.) then



          open(unit=257, file='totterms', status='old', action='write', &
               position='append')
          write(257,*) 'T'
          close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

!              write(*,*) 'Getting values from cache'
!              write(*,*) ' '

             open(unit=257, file='cachehit', status='old', action='write', &
                  position='append') 
             write(257,*) 'T'
             close(257)

             call property_cache_getdata(cache, 2, p12, property_size, prop)
       
          else

          write(*,*) 'Calculating Pulay k-n contribution:'
          write(*,*) 'S', p12(1)%pid
          write(*,*) 'W', p12(2)%pid


             call get_pulay_kn((/ (p_tuple_standardorder(p12(i)) , i = 1, 2)  /), & 
                               kn, F, D, S, property_size, cache, prop)

             write(*,*) 'Calculated Pulay k-n contribution'
             write(*,*) ' '

          end if

       else

!           write(*,*) 'Pulay k-n contribution was k-n skipped:'
!           write(*,*) 'S ', p12(1)%pid 
!           write(*,*) 'W ', p12(2)%pid 
!           write(*,*) ' '

       end if 

    end if

  end subroutine


  subroutine get_pulay_kn(p12, kn, F, D, S, property_size, cache, prop)

    implicit none

    type(p_tuple) :: pert, emptypert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structb
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    type(matrix) :: W
    integer :: i, j, k, sstr_incr, offset, total_num_perturbations, &
               property_size, dtup_ind, pr_offset, ca_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, merged_nblks
    integer, dimension(p12(1)%n_perturbations + p12(2)%n_perturbations) :: & 
    pids_current_contribution, translated_index
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:) :: ncinnersmall, blk_sizes_merged
    integer, allocatable, dimension(:,:) :: triang_indices_pr, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    integer :: d_supsize
    integer, dimension(2) :: kn
    integer, allocatable, dimension(:) :: ncarray, ncinner, inner_offsets, &
                                          which_index_is_pid
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    complex(8), allocatable, dimension(:) :: tmp, prop_forcache
    complex(8), dimension(property_size) :: prop

    d_supsize = derivative_superstructure_getsize(p12(2), kn, .FALSE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/))

    allocate(deriv_structb(d_supsize, 3))

    sstr_incr = 0

    call derivative_superstructure(p12(2), kn, .FALSE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize, sstr_incr, deriv_structb)

    allocate(nfields(2))
    allocate(nblks_tuple(2))
    
    
    
    do i = 1, 2
    
       nfields(i) = p12(i)%n_perturbations
       nblks_tuple(i) = get_num_blks(p12(i))
    
    end do
    
    total_num_perturbations = sum(nfields)
    
    allocate(blks_tuple_info(2, total_num_perturbations, 3))
    allocate(blks_tuple_triang_size(2))
    allocate(blk_sizes(2, total_num_perturbations))
    allocate(blk_sizes_merged(total_num_perturbations))
    
    do i = 1, 2
    
       blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p12(i))
    
       blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                                   blks_tuple_info(i, 1:nblks_tuple(i), :))
    
       blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
       blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))

    end do
    
    if (p12(2)%n_perturbations == 0) then
    
       outer_indices_size = 1
    
    else
    
       outer_indices_size = blks_tuple_triang_size(2)
    
    end if
    
    inner_indices_size = blks_tuple_triang_size(1)
    
    allocate(prop_forcache(inner_indices_size * outer_indices_size))
    allocate(ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(ncinner(p12(1)%n_perturbations))
    allocate(tmp(inner_indices_size))
    allocate(inner_offsets(inner_indices_size))
    allocate(outer_indices(outer_indices_size, p12(2)%n_perturbations))
    allocate(inner_indices(inner_indices_size, p12(1)%n_perturbations))
    allocate(which_index_is_pid(p12(1)%n_perturbations + p12(2)%n_perturbations))

    prop_forcache = 0.0
    
    ncarray = get_ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations, 2, p12)
    ncinner = nc_onlysmall(p12(1)%n_perturbations + p12(2)%n_perturbations, &
                           p12(1)%n_perturbations, 1, p12(1), ncarray)
    
    which_index_is_pid = 0
    
    do i = 1, p12(2)%n_perturbations
    
       which_index_is_pid(p12(2)%pid(i)) = i
    
    end do
    
    
    if (p12(2)%n_perturbations > 0) then
    
       call make_triangulated_indices(nblks_tuple(2), blks_tuple_info(2, &
            1:nblks_tuple(2), :), blks_tuple_triang_size(2), outer_indices)
    
    end if
    
    call make_triangulated_indices(nblks_tuple(1), blks_tuple_info(1, &
         1:nblks_tuple(1), :), blks_tuple_triang_size(1), inner_indices)
    
    ! ASSUME CLOSED SHELL
    call mat_init(W, zeromat%nrow, zeromat%ncol, &
                  .false., .false., .false., .false., .false.)
    call mat_init_like_and_zero(zeromat, W)
    
    do i = 1, size(outer_indices, 1)
    
       tmp = 0.0
       W = zeromat
       call mat_ensure_alloc(W)
 
       call rsp_get_matrix_w(zeromat, d_supsize, deriv_structb, p12(1)%n_perturbations + &
                             p12(2)%n_perturbations, which_index_is_pid, &
                             p12(2)%n_perturbations, outer_indices(i,:), F, D, S, W)
 
       call rsp_ovlave_tr(p12(1)%n_perturbations, p12(1)%plab, &
                          (/ (j/j, j = 1, p12(1)%n_perturbations) /), &
                          p12(1)%pdim, nblks_tuple(1), blks_tuple_info(1, &
                          1:nblks_tuple(1), :), blk_sizes(1, 1:nblks_tuple(1)), &
                          size(tmp), W, tmp)
 
       do j = 1, size(inner_indices, 1)
    
          if (p12(2)%n_perturbations > 0) then
    
             offset = get_triang_blks_tuple_offset(2, total_num_perturbations, nblks_tuple, &
                      (/ p12(1)%n_perturbations, p12(2)%n_perturbations /), &
                      blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                      (/inner_indices(j, :), outer_indices(i, :) /)) 
          else
    
             offset = get_triang_blks_tuple_offset(1, total_num_perturbations, &
                      nblks_tuple(1), (/ p12(1)%n_perturbations /), &
                      blks_tuple_info(1,:,:), blk_sizes(1,:), blks_tuple_triang_size(1), &
                      (/inner_indices(j, :) /)) 

          end if
   
          prop_forcache(offset) = prop_forcache(offset) + tmp(j)

       end do
    
    end do
    
    call p_tuple_p1_cloneto_p2(p12(1), merged_p_tuple)
    
    if (p12(2)%n_perturbations > 0) then
    
       merged_p_tuple = merge_p_tuple(merged_p_tuple, p12(2))
    
    end if
    
    merged_p_tuple = p_tuple_standardorder(merged_p_tuple)
    
    ! MaR: NOT DOING THE FOLLOWING FOR THE MERGED PERT ASSUMES THAT 
    ! PIDS ARE IN STANDARD ORDER? FIND OUT

    k = 1
    do i = 1, 2
       do j = 1, p12(i)%n_perturbations
          pids_current_contribution(k) = p12(i)%pid(j)
          k = k + 1
       end do
    end do
    
    merged_nblks = get_num_blks(merged_p_tuple)

    allocate(merged_blk_info(1, merged_nblks, 3))

    merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)
    blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
    merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
    merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)
    
    allocate(triang_indices_pr(merged_triang_size, sum(merged_blk_info(1, :,2))))
    
    call make_triangulated_indices(merged_nblks, merged_blk_info, & 
         merged_triang_size, triang_indices_pr)
    
    do i = 1, size(triang_indices_pr, 1)
    
       pr_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
                   (/sum(nfields)/), &
                   (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
                   (/triang_indices_pr(i, :) /))
   
       do j = 1, total_num_perturbations
    
          translated_index(j) = triang_indices_pr(i,pids_current_contribution(j))
    
       end do
    
       if (p12(2)%n_perturbations > 0) then
    
          ca_offset = get_triang_blks_tuple_offset(2, &
                      total_num_perturbations, nblks_tuple, &
                      nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                      (/ translated_index(:) /))
    
       else
    
          ca_offset = get_triang_blks_tuple_offset(1, &
                      total_num_perturbations, nblks_tuple(1), &
                      nfields(1), blks_tuple_info(1, :, :), &
                      blk_sizes(1,:), blks_tuple_triang_size(1), & 
                      (/ translated_index(:) /))
    
       end if
    
       prop(pr_offset) = prop(pr_offset) + prop_forcache(ca_offset)
    
    end do
    

!     write(*,*) 'pulay kn contribution'
!     call print_rsp_tensor_stdout_tr(1, total_num_perturbations, merged_p_tuple%pdim, &
!     (/ (1, j = 1, (merged_p_tuple%n_perturbations - 1) ) /), merged_nblks, blk_sizes_merged, &
!     merged_blk_info, property_size, prop_forcache)

    call property_cache_add_element(cache, 2, p12,  &
         inner_indices_size * outer_indices_size, prop_forcache)    
   
    W = 0 
    deallocate(deriv_structb)
    deallocate(ncarray)
    deallocate(ncinner)
    deallocate(tmp)
    deallocate(inner_offsets)
    deallocate(outer_indices)
    deallocate(inner_indices)
    deallocate(which_index_is_pid)
    deallocate(nfields)
    deallocate(nblks_tuple)
    deallocate(blks_tuple_info)
    deallocate(blks_tuple_triang_size)
    deallocate(blk_sizes)
    deallocate(blk_sizes_merged)
    deallocate(prop_forcache)

  end subroutine

  recursive subroutine rsp_pulay_lag(pert, kn, p12, S, D, F, property_size, cache, prop)

    implicit none

    type(p_tuple) :: pert
    type(p_tuple), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    integer :: property_size, i
    integer, dimension(2) :: kn
    complex(8), dimension(property_size) :: prop
    
    if (pert%n_perturbations > 0) then

       call rsp_pulay_lag(p_tuple_remove_first(pert), kn, &
       (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), &
       S, D, F, property_size, cache, prop)
       call rsp_pulay_lag(p_tuple_remove_first(pert), kn, &
       (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), &
       S, D, F, property_size, cache, prop)

    else

       ! At lowest level:
       if (kn_skip(p12(1)%n_perturbations, p12(1)%pid, kn) .EQV. .FALSE.) then



       open(unit=257, file='totterms', status='old', action='write', position='append') 
       write(257,*) 'T'
       close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

!              write(*,*) 'Getting values from cache'
!              write(*,*) ' '
       
             open(unit=257, file='cachehit', status='old', action='write', &
             position='append') 
             write(257,*) 'T'
             close(257)

             call property_cache_getdata(cache, 2, p12, property_size, prop)

          else

       write(*,*) 'Calculating Pulay lagrange contribution:'
       write(*,*) 'S', p12(1)%pid
       write(*,*) 'W', p12(2)%pid, 'primed', kn(2)

             call get_pulay_lag((/ (p_tuple_standardorder(p12(i)) , i = 1, 2) /), & 
                               kn, F, D, S, property_size, cache, prop)

             write(*,*) 'Calculated Pulay lagrange contribution'
             write(*,*) ' '

          end if

       else

!           write(*,*) 'Pulay lagrange contribution was k-n skipped:'
!           write(*,*) 'S', p12(1)%pid 
!           write(*,*) 'W', p12(2)%pid, 'primed', kn(2)
!           write(*,*) ' '

       end if

    end if

  end subroutine


  subroutine get_pulay_lag(p12, kn, F, D, S, property_size, cache, prop)

    implicit none

    type(p_tuple) :: pert, emptypert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structb
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    type(matrix) :: W
    integer :: i, j, k ,m, incr
    integer :: d_supsize, total_num_perturbations, &
             property_size, offset, dtup_ind, pr_offset, ca_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, merged_nblks
    integer, dimension(p12(1)%n_perturbations + p12(2)%n_perturbations) :: & 
    pids_current_contribution, translated_index
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:) :: ncinnersmall, blk_sizes_merged
    integer, allocatable, dimension(:,:) :: triang_indices_pr, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    integer, dimension(2) :: kn
    integer, allocatable, dimension(:) :: ncarray, ncinner, inner_offsets, &
                                          which_index_is_pid
    integer, allocatable, dimension(:) :: outer_ind_b_large
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    complex(8), allocatable, dimension(:) :: tmp, prop_forcache
    complex(8), dimension(property_size) :: prop

    d_supsize = derivative_superstructure_getsize(p12(2), kn, .TRUE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/))
   
    allocate(deriv_structb(d_supsize, 3))

    incr = 0

    call derivative_superstructure(p12(2), kn, .TRUE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize, incr, deriv_structb)


    allocate(nfields(2))
    allocate(nblks_tuple(2))

    do i = 1, 2

       nfields(i) = p12(i)%n_perturbations
       nblks_tuple(i) = get_num_blks(p12(i))

    end do

    total_num_perturbations = sum(nfields)

    allocate(blks_tuple_info(2, total_num_perturbations, 3))
    allocate(blks_tuple_triang_size(2))
    allocate(blk_sizes(2, total_num_perturbations))
    allocate(blk_sizes_merged(total_num_perturbations))

    do i = 1, 2

       blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p12(i))
       blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                                   blks_tuple_info(i, 1:nblks_tuple(i), :))
       blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
       blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))

    end do

    if (p12(2)%n_perturbations == 0) then

       outer_indices_size = 1

    else

       outer_indices_size = blks_tuple_triang_size(2)

    end if

    inner_indices_size = blks_tuple_triang_size(1)

    allocate(prop_forcache(inner_indices_size * outer_indices_size))
    allocate(ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(ncinner(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_ind_b_large(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(tmp(inner_indices_size))
    allocate(inner_offsets(inner_indices_size))
    allocate(outer_indices(outer_indices_size, p12(2)%n_perturbations))
    allocate(inner_indices(inner_indices_size, p12(1)%n_perturbations))
    allocate(which_index_is_pid(p12(1)%n_perturbations + p12(2)%n_perturbations))

    prop_forcache = 0.0

    ncarray = get_ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations, 2, p12)
    ncinner = nc_only(p12(1)%n_perturbations + p12(2)%n_perturbations, &
              p12(1)%n_perturbations, 1, p12(1), ncarray)

    which_index_is_pid = 0

    do i = 1, p12(2)%n_perturbations

       which_index_is_pid(p12(2)%pid(i)) = i

    end do

    call make_triangulated_indices(nblks_tuple(2), blks_tuple_info(2, &
         1:nblks_tuple(2), :), blks_tuple_triang_size(2), outer_indices)


    call make_triangulated_indices(nblks_tuple(1), blks_tuple_info(1, &
         1:nblks_tuple(1), :), blks_tuple_triang_size(1), inner_indices)

    ! ASSUME CLOSED SHELL
    call mat_init(W, zeromat%nrow, zeromat%ncol, &
                  .false., .false., .false., .false., .false.)
    call mat_init_like_and_zero(zeromat, W)

    do i = 1, size(outer_indices, 1)

       tmp = 0.0
       W = zeromat
       call mat_ensure_alloc(W)

       call rsp_get_matrix_w(zeromat, d_supsize, deriv_structb, p12(1)%n_perturbations + &
                            p12(2)%n_perturbations, which_index_is_pid, &
                            p12(2)%n_perturbations, outer_indices(i,:), F, D, S, W)

       call rsp_ovlave_tr(p12(1)%n_perturbations, p12(1)%plab, &
                       (/ (j/j, j = 1, p12(1)%n_perturbations) /), &
                       p12(1)%pdim, nblks_tuple(1), blks_tuple_info(1, &
                       1:nblks_tuple(1), :), blk_sizes(1, 1:nblks_tuple(1)), &
                       size(tmp), W, tmp)

       do j = 1, size(inner_indices, 1)

          offset = get_triang_blks_tuple_offset(2, total_num_perturbations, nblks_tuple, &
                   (/ p12(1)%n_perturbations, p12(2)%n_perturbations /), &
                   blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                   (/inner_indices(j, :), outer_indices(i, :) /)) 

          prop_forcache(offset) = prop_forcache(offset) + tmp(j)

       end do

    end do

    call p_tuple_p1_cloneto_p2(p12(1), merged_p_tuple)

    if (p12(2)%n_perturbations > 0) then

       merged_p_tuple = merge_p_tuple(merged_p_tuple, p12(2))

    end if

    merged_p_tuple = p_tuple_standardorder(merged_p_tuple)

    ! MR: NOT DOING THE FOLLOWING FOR THE MERGED PERT ASSUMES THAT 
    ! PIDS ARE IN STANDARD ORDER? FIND OUT

    k = 1
    do i = 1, 2
       do j = 1, p12(i)%n_perturbations
          pids_current_contribution(k) = p12(i)%pid(j)
       k = k + 1
       end do
    end do

    merged_nblks = get_num_blks(merged_p_tuple)

    allocate(merged_blk_info(1, merged_nblks, 3))

    merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)
    blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
    merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
    merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)

    allocate(triang_indices_pr(merged_triang_size, sum(merged_blk_info(1, :,2))))

    call make_triangulated_indices(merged_nblks, merged_blk_info, & 
         merged_triang_size, triang_indices_pr)

    do i = 1, size(triang_indices_pr, 1)

       pr_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
                   (/sum(nfields)/), &
                   (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
                   (/triang_indices_pr(i, :) /))

       do j = 1, total_num_perturbations

          translated_index(j) = triang_indices_pr(i,pids_current_contribution(j))

       end do

       if (p12(2)%n_perturbations > 0) then

          ca_offset = get_triang_blks_tuple_offset(2, &
                      total_num_perturbations, nblks_tuple, &
                      nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                      (/ translated_index(:) /))

       else

          ca_offset = get_triang_blks_tuple_offset(1, &
                      total_num_perturbations, nblks_tuple(1), &
                      nfields(1), blks_tuple_info(1, :, :), &
                      blk_sizes(1,:), blks_tuple_triang_size(1), & 
                      (/ translated_index(:) /))

       end if

       prop(pr_offset) = prop(pr_offset) + prop_forcache(ca_offset)

    end do


!     write(*,*) 'pulay lag contribution'
!     call print_rsp_tensor_stdout_tr(1, total_num_perturbations, merged_p_tuple%pdim, &
!     (/ (1, j = 1, (merged_p_tuple%n_perturbations - 1) ) /), merged_nblks, blk_sizes_merged, &
!     merged_blk_info, property_size, prop_forcache)

    call property_cache_add_element(cache, 2, p12,  &
         inner_indices_size * outer_indices_size, prop_forcache)    

    deallocate(deriv_structb)
    deallocate(ncarray)
    deallocate(ncinner)
    deallocate(outer_ind_b_large)
    deallocate(tmp)
    deallocate(inner_offsets)
    deallocate(outer_indices)
    deallocate(inner_indices)
    deallocate(which_index_is_pid)
    deallocate(nfields)
    deallocate(nblks_tuple)
    deallocate(blks_tuple_info)
    deallocate(blks_tuple_triang_size)
    deallocate(blk_sizes)
    deallocate(blk_sizes_merged)
    deallocate(triang_indices_pr)
    deallocate(merged_blk_info)
    deallocate(prop_forcache)
    W = 0

  end subroutine

  recursive subroutine rsp_idem_lag(pert, kn, p12, S, D, F, &
                                    property_size, cache, prop)

    implicit none

    type(p_tuple) :: pert
    type(p_tuple), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    integer :: property_size, i
    integer, dimension(2) :: kn
    complex(8), dimension(property_size) :: prop
    
    if (pert%n_perturbations > 0) then

       call rsp_idem_lag(p_tuple_remove_first(pert), kn, &
       (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), S, D, F, property_size, &
       cache, prop)
       call rsp_idem_lag(p_tuple_remove_first(pert), kn, &
       (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), S, D, F, property_size, &
       cache, prop)

    else

       if (kn_skip(p12(1)%n_perturbations, p12(1)%pid, kn) .EQV. .FALSE.) then



          open(unit=257, file='totterms', status='old', action='write', &
               position='append') 
          write(257,*) 'T'
          close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

!              write(*,*) 'Getting values from cache'
!              write(*,*) ' '

             open(unit=257, file='cachehit', status='old', action='write', &
                  position='append')
             write(257,*) 'T'
             close(257)

             call property_cache_getdata(cache, 2, p12, property_size, prop)
      
          else

          write(*,*) 'Calculating idempotency lagrange contribution'
          write(*,*) 'Zeta', p12(1)%pid
          write(*,*) 'Z', p12(2)%pid, 'primed', kn(2)

             ! At lowest level:
             call get_idem_lag((/ (p_tuple_standardorder(p12(i)) , i = 1, 2) /), & 
                               kn, F, D, S, property_size, cache, prop)

             write(*,*) 'Calculated idempotency lagrange contribution'
             write(*,*) ' '

          end if

       else

!           write(*,*) 'Idempotency lagrange contribution was k-n skipped:'
!           write(*,*) 'Zeta', p12(1)%pid 
!           write(*,*) 'Z', p12(2)%pid, 'primed', kn(2)
!           write(*,*) ' '

       end if

    end if

  end subroutine


  subroutine get_idem_lag(p12, kn, F, D, S, property_size, cache, prop)

    implicit none

    type(p_tuple) :: pert, emptypert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structa, deriv_structb
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    type(matrix) :: Zeta, Z
    integer :: i, j, k, m, n, p, incr1, incr2, total_num_perturbations, &
             property_size, offset, dtup_ind, pr_offset, ca_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, merged_nblks
    integer, dimension(p12(1)%n_perturbations + p12(2)%n_perturbations) :: & 
    pids_current_contribution, translated_index
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:) :: ncinnersmall, blk_sizes_merged
    integer, allocatable, dimension(:,:) :: triang_indices_pr, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    integer, dimension(2) :: kn, d_supsize
    integer, allocatable, dimension(:) :: ncarray, ncinner, which_index_is_pid1, &
                                          which_index_is_pid2
    integer, allocatable, dimension(:) :: outer_ind_a_large, outer_ind_b_large
    integer, allocatable, dimension(:,:) :: outer_indices_a, outer_indices_b
    complex(8), dimension(property_size) :: prop
    complex(8), allocatable, dimension(:) :: prop_forcache

    d_supsize = 0
    d_supsize(1) = derivative_superstructure_getsize(p_tuple_remove_first(p12(1)), &
                   kn, .FALSE., (/get_emptypert(), get_emptypert(), get_emptypert()/))
    d_supsize(2) = derivative_superstructure_getsize(p12(2), &
                   kn, .TRUE., (/get_emptypert(), get_emptypert(), get_emptypert()/))

    allocate(deriv_structa(d_supsize(1), 3))
    allocate(deriv_structb(d_supsize(2), 3))

    incr1 = 0
    incr2 = 0

    call derivative_superstructure(p_tuple_remove_first(p12(1)), kn, .FALSE., & 
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize(1), incr1, deriv_structa)
    call derivative_superstructure(p12(2), kn, .TRUE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize(2), incr2, deriv_structb)

    allocate(nfields(2))
    allocate(nblks_tuple(2))

    do i = 1, 2

       nfields(i) = p12(i)%n_perturbations
       nblks_tuple(i) = get_num_blks(p12(i))

    end do

    total_num_perturbations = sum(nfields)

    allocate(blks_tuple_info(2, total_num_perturbations, 3))
    allocate(blks_tuple_triang_size(2))
    allocate(blk_sizes(2, total_num_perturbations))
    allocate(blk_sizes_merged(total_num_perturbations))

    do i = 1, 2

       blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p12(i))
       blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                                   blks_tuple_info(i, 1:nblks_tuple(i), :))
       blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
       blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))

    end do

    if (p12(2)%n_perturbations == 0) then

       outer_indices_size = 1

    else

       outer_indices_size = blks_tuple_triang_size(2)

    end if

    inner_indices_size = blks_tuple_triang_size(1)

    allocate(prop_forcache(inner_indices_size * outer_indices_size))
    allocate(ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(ncinner(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_ind_a_large(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_ind_b_large(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_indices_a(inner_indices_size, p12(1)%n_perturbations))
    allocate(outer_indices_b(outer_indices_size, p12(2)%n_perturbations))
    allocate(which_index_is_pid1(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(which_index_is_pid2(p12(1)%n_perturbations + p12(2)%n_perturbations))

    prop_forcache = 0.0

    ncarray = get_ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations, 2, p12)
    ncinner = nc_only(p12(1)%n_perturbations + p12(2)%n_perturbations, &
                      p12(1)%n_perturbations, 1, p12(1), ncarray)

    which_index_is_pid1 = 0

    do i = 1, p12(1)%n_perturbations

       which_index_is_pid1(p12(1)%pid(i)) = i

    end do

    which_index_is_pid2 = 0

    do i = 1, p12(2)%n_perturbations

       which_index_is_pid2(p12(2)%pid(i)) = i

    end do

    call make_triangulated_indices(nblks_tuple(1), blks_tuple_info(1, &
         1:nblks_tuple(1), :), blks_tuple_triang_size(1), outer_indices_a)


    call make_triangulated_indices(nblks_tuple(2), blks_tuple_info(2, &
         1:nblks_tuple(2), :), blks_tuple_triang_size(2), outer_indices_b)

    offset = 0.0

    ! ASSUME CLOSED SHELL
    call mat_init(Z, zeromat%nrow, zeromat%ncol, &
                  .false., .false., .false., .false., .false.)
    call mat_init_like_and_zero(zeromat, Z)

    ! ASSUME CLOSED SHELL
    call mat_init(Zeta, zeromat%nrow, zeromat%ncol, &
                  .false., .false., .false., .false., .false.)
    call mat_init_like_and_zero(zeromat, Zeta)

    do i = 1, size(outer_indices_a, 1)

       Zeta = zeromat
       call mat_ensure_alloc(Zeta)

       call rsp_get_matrix_zeta(zeromat, p_tuple_getone(p12(1), 1), kn, d_supsize(1), &
            deriv_structa, p12(1)%n_perturbations + p12(2)%n_perturbations, &
            which_index_is_pid1, p12(1)%n_perturbations, outer_indices_a(i,:), &
            F, D, S, Zeta)

       do j = 1, size(outer_indices_b, 1)

          Z = zeromat
          call mat_ensure_alloc(Z)

          call rsp_get_matrix_z(zeromat, d_supsize(2), deriv_structb, kn, &
               p12(1)%n_perturbations + p12(2)%n_perturbations, which_index_is_pid2, &
               p12(2)%n_perturbations, outer_indices_b(j,:), F, D, S, Z)

          offset = get_triang_blks_tuple_offset(2, total_num_perturbations, nblks_tuple, &
                   (/ p12(1)%n_perturbations, p12(2)%n_perturbations /), &
                   blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                   (/outer_indices_a(i, :), outer_indices_b(j, :) /)) 

          prop_forcache(offset) = prop_forcache(offset) -trace(Zeta, Z)

       end do

    end do

    call p_tuple_p1_cloneto_p2(p12(1), merged_p_tuple)

    if (p12(2)%n_perturbations > 0) then

       merged_p_tuple = merge_p_tuple(merged_p_tuple, p12(2))

    end if

    merged_p_tuple = p_tuple_standardorder(merged_p_tuple)

! MR: NOT DOING THE FOLLOWING FOR THE MERGED PERT ASSUMES THAT 
! PIDS ARE IN STANDARD ORDER? FIND OUT

    k = 1
    do i = 1, 2
       do j = 1, p12(i)%n_perturbations
          pids_current_contribution(k) = p12(i)%pid(j)
          k = k + 1
       end do
    end do

    merged_nblks = get_num_blks(merged_p_tuple)

    allocate(merged_blk_info(1, merged_nblks, 3))

    merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)
    blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
    merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
    merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)

    allocate(triang_indices_pr(merged_triang_size, sum(merged_blk_info(1, :,2))))

    call make_triangulated_indices(merged_nblks, merged_blk_info, & 
         merged_triang_size, triang_indices_pr)

    do i = 1, size(triang_indices_pr, 1)

       pr_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
                   (/sum(nfields)/), &
                   (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
                   (/triang_indices_pr(i, :) /))

       do j = 1, total_num_perturbations

          translated_index(j) = triang_indices_pr(i,pids_current_contribution(j))

       end do

       if (p12(2)%n_perturbations > 0) then

          ca_offset = get_triang_blks_tuple_offset(2, &
                      total_num_perturbations, nblks_tuple, &
                      nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                      (/ translated_index(:) /))

       else

          ca_offset = get_triang_blks_tuple_offset(1, &
                      total_num_perturbations, nblks_tuple(1), &
                      nfields(1), blks_tuple_info(1, :, :), &
                      blk_sizes(1,:), blks_tuple_triang_size(1), & 
                      (/ translated_index(:) /))

       end if

       prop(pr_offset) = prop(pr_offset) + prop_forcache(ca_offset)

    end do

!     write(*,*) 'idempotency contribution'
!  call print_rsp_tensor_stdout_tr(1, total_num_perturbations, merged_p_tuple%pdim, &
!  (/ (1, j = 1, (merged_p_tuple%n_perturbations - 1) ) /), merged_nblks, blk_sizes_merged, &
!  merged_blk_info, property_size, prop_forcache)

    call property_cache_add_element(cache, 2, p12,  &
         inner_indices_size * outer_indices_size, prop_forcache)    

    deallocate(deriv_structa)
    deallocate(deriv_structb)
    deallocate(ncarray)
    deallocate(ncinner)
    deallocate(outer_ind_a_large)
    deallocate(outer_ind_b_large)
    deallocate(outer_indices_a)
    deallocate(outer_indices_b)
    deallocate(which_index_is_pid1)
    deallocate(which_index_is_pid2)
    deallocate(nfields)
    deallocate(nblks_tuple)
    deallocate(blks_tuple_info)
    deallocate(blks_tuple_triang_size)
    deallocate(blk_sizes)
    deallocate(blk_sizes_merged)
    deallocate(prop_forcache)
    Zeta = 0
    Z = 0

  end subroutine



  recursive subroutine rsp_scfe_lag(pert, kn, p12, S, D, F, &
                                    property_size, cache, prop)

    implicit none

    type(p_tuple) :: pert
    type(p_tuple), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    integer :: property_size, i
    integer, dimension(2) :: kn
    complex(8), dimension(property_size) :: prop
    
    if (pert%n_perturbations > 0) then

       call rsp_scfe_lag(p_tuple_remove_first(pert), kn, &
            (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), &
            S, D, F, property_size, cache, prop)
       call rsp_scfe_lag(p_tuple_remove_first(pert), kn, &
            (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), &
            S, D, F, property_size, cache, prop)

    else

       if (kn_skip(p12(1)%n_perturbations, p12(1)%pid, kn) .EQV. .FALSE.) then



          open(unit=257, file='totterms', status='old', action='write', &
               position='append') 
          write(257,*) 'T'
          close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

             open(unit=257, file='cachehit', status='old', action='write', &
                  position='append') 
             write(257,*) 'T'
             close(257)

!              write(*,*) 'Getting values from cache'
!              write(*,*) ' '

             call property_cache_getdata(cache, 2, p12, property_size, prop)
       
          else

          write(*,*) 'Calculating scfe lagrange contribution'
          write(*,*) 'Lambda', p12(1)%pid
          write(*,*) 'Y', p12(2)%pid, 'primed', kn(2)

             ! At lowest level:
             call get_scfe_lag((/ (p_tuple_standardorder(p12(i)) , i = 1, 2) /), &
             kn, F, D, S, property_size, cache, prop)

             write(*,*) 'Calculated scfe lagrange contribution'
             write(*,*) ' '

          end if

       else

!           write(*,*) 'scfe lagrange contribution was k-n skipped:'
!           write(*,*) 'Lambda', p12(1)%pid 
!           write(*,*) 'Y', p12(2)%pid, 'primed', kn(2)
!           write(*,*) ' '

       end if

    end if

  end subroutine


  subroutine get_scfe_lag(p12, kn, F, D, S, property_size, cache, prop)

    implicit none

    type(p_tuple) :: pert, emptypert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structa, deriv_structb
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    type(matrix) :: L, Y
    integer :: i, j, k, m, n, p, incr1, incr2, total_num_perturbations, &
             property_size, offset, dtup_ind, pr_offset, ca_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, merged_nblks
    integer, dimension(p12(1)%n_perturbations + p12(2)%n_perturbations) :: & 
    pids_current_contribution, translated_index
    integer, allocatable, dimension(:) :: ncinnersmall, blk_sizes_merged, nfields, nblks_tuple
    integer, allocatable, dimension(:) :: blks_tuple_triang_size
    integer, allocatable, dimension(:,:) :: triang_indices_pr, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    integer, dimension(2) :: kn, d_supsize
    integer, allocatable, dimension(:) :: ncarray, ncinner, which_index_is_pid1, which_index_is_pid2
    integer, allocatable, dimension(:) :: outer_ind_a_large, outer_ind_b_large
    integer, allocatable, dimension(:,:) :: outer_indices_a, outer_indices_b
    complex(8), dimension(property_size) :: prop
    complex(8), allocatable, dimension(:) :: prop_forcache

    d_supsize = 0

    d_supsize(1) = derivative_superstructure_getsize(p_tuple_remove_first(p12(1)), &
                   kn, .FALSE., (/get_emptypert(), get_emptypert(), get_emptypert()/))
    d_supsize(2) = derivative_superstructure_getsize(p12(2), &
                   kn, .TRUE., (/get_emptypert(), get_emptypert(), get_emptypert()/))

    allocate(deriv_structa(d_supsize(1), 3))
    allocate(deriv_structb(d_supsize(2), 3))

    incr1 = 0
    incr2 = 0

    call derivative_superstructure(p_tuple_remove_first(p12(1)), kn, .FALSE., &
                    (/get_emptypert(), get_emptypert(), get_emptypert()/), & 
                    d_supsize(1), incr1, deriv_structa)
    call derivative_superstructure(p12(2), kn, .TRUE., &
                    (/get_emptypert(), get_emptypert(), get_emptypert()/), &
                    d_supsize(2), incr2, deriv_structb)

    allocate(nfields(2))
    allocate(nblks_tuple(2))

    do i = 1, 2

       nfields(i) = p12(i)%n_perturbations
       nblks_tuple(i) = get_num_blks(p12(i))

    end do

    total_num_perturbations = sum(nfields)

    allocate(blks_tuple_info(2, total_num_perturbations, 3))
    allocate(blks_tuple_triang_size(2))
    allocate(blk_sizes(2, total_num_perturbations))
    allocate(blk_sizes_merged(total_num_perturbations))

    do i = 1, 2

       blks_tuple_info(i, :, :) = get_blk_info(nblks_tuple(i), p12(i))
       blks_tuple_triang_size(i) = get_triangulated_size(nblks_tuple(i), &
                                   blks_tuple_info(i, 1:nblks_tuple(i), :))
       blk_sizes(i, 1:nblks_tuple(i)) = get_triangular_sizes(nblks_tuple(i), &
       blks_tuple_info(i,1:nblks_tuple(i),2), blks_tuple_info(i,1:nblks_tuple(i),3))

    end do

    if (p12(2)%n_perturbations == 0) then

       outer_indices_size = 1

    else

       outer_indices_size = blks_tuple_triang_size(2)

    end if

    inner_indices_size = blks_tuple_triang_size(1)

    allocate(prop_forcache(inner_indices_size * outer_indices_size))
    allocate(ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(ncinner(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_indices_a(inner_indices_size, p12(1)%n_perturbations))
    allocate(outer_indices_b(outer_indices_size, p12(2)%n_perturbations))
    allocate(outer_ind_a_large(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_ind_b_large(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(which_index_is_pid1(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(which_index_is_pid2(p12(1)%n_perturbations + p12(2)%n_perturbations))

    prop_forcache = 0.0

    ncarray = get_ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations, 2, p12)
    ncinner = nc_only(p12(1)%n_perturbations + p12(2)%n_perturbations, &
              p12(1)%n_perturbations, 1, p12(1), ncarray)

    which_index_is_pid1 = 0

    do i = 1, p12(1)%n_perturbations

       which_index_is_pid1(p12(1)%pid(i)) = i

    end do

    which_index_is_pid2 = 0

    do i = 1, p12(2)%n_perturbations

       which_index_is_pid2(p12(2)%pid(i)) = i

    end do

    call make_triangulated_indices(nblks_tuple(1), blks_tuple_info(1, &
         1:nblks_tuple(1), :), blks_tuple_triang_size(1), outer_indices_a)

    call make_triangulated_indices(nblks_tuple(2), blks_tuple_info(2, &
         1:nblks_tuple(2), :), blks_tuple_triang_size(2), outer_indices_b)

    offset = 0

    ! ASSUME CLOSED SHELL
    call mat_init(Y, zeromat%nrow, zeromat%ncol, &
                  .false., .false., .false., .false., .false.)
    call mat_init_like_and_zero(zeromat, Y)

    ! ASSUME CLOSED SHELL
    call mat_init(L, zeromat%nrow, zeromat%ncol, &
                  .false., .false., .false., .false., .false.)
    call mat_init_like_and_zero(zeromat, L)

    do i = 1, size(outer_indices_a, 1)

       L = zeromat
       call mat_ensure_alloc(L)

       call rsp_get_matrix_lambda(zeromat, p_tuple_getone(p12(1), 1), d_supsize(1), &
            deriv_structa, p12(1)%n_perturbations + p12(2)%n_perturbations, &
            which_index_is_pid1, p12(1)%n_perturbations, outer_indices_a(i,:), D, S, L)

       do j = 1, size(outer_indices_b, 1)

          Y = zeromat
          call mat_ensure_alloc(Y)

          call rsp_get_matrix_y(zeromat, d_supsize(2), deriv_structb, &
               p12(1)%n_perturbations + p12(2)%n_perturbations, which_index_is_pid2, &
               p12(2)%n_perturbations, outer_indices_b(j,:), F, D, S, Y)


          offset = get_triang_blks_tuple_offset(2, total_num_perturbations, nblks_tuple, &
                   (/ p12(1)%n_perturbations, p12(2)%n_perturbations /), &
                   blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                   (/outer_indices_a(i, :), outer_indices_b(j, :) /)) 

          prop_forcache(offset) = prop_forcache(offset) - trace(L, Y)

       end do

    end do

    call p_tuple_p1_cloneto_p2(p12(1), merged_p_tuple)

    if (p12(2)%n_perturbations > 0) then

       merged_p_tuple = merge_p_tuple(merged_p_tuple, p12(2))

    end if

    merged_p_tuple = p_tuple_standardorder(merged_p_tuple)

    ! MR: NOT DOING THE FOLLOWING FOR THE MERGED PERT ASSUMES THAT 
    ! PIDS ARE IN STANDARD ORDER? FIND OUT

    k = 1
    do i = 1, 2
       do j = 1, p12(i)%n_perturbations
          pids_current_contribution(k) = p12(i)%pid(j)
          k = k + 1
       end do
    end do

    merged_nblks = get_num_blks(merged_p_tuple)

    allocate(merged_blk_info(1, merged_nblks, 3))

    merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)
    blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
    merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
    merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)

    allocate(triang_indices_pr(merged_triang_size, sum(merged_blk_info(1, :,2))))

    call make_triangulated_indices(merged_nblks, merged_blk_info, & 
         merged_triang_size, triang_indices_pr)

    do i = 1, size(triang_indices_pr, 1)

       pr_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
                   (/sum(nfields)/), &
                   (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
                   (/triang_indices_pr(i, :) /))

       do j = 1, total_num_perturbations

          translated_index(j) = triang_indices_pr(i,pids_current_contribution(j))

       end do

       if (p12(2)%n_perturbations > 0) then

          ca_offset = get_triang_blks_tuple_offset(2, &
                      total_num_perturbations, nblks_tuple, &
                      nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
                      (/ translated_index(:) /))

       else

          ca_offset = get_triang_blks_tuple_offset(1, &
                      total_num_perturbations, nblks_tuple(1), &
                      nfields(1), blks_tuple_info(1, :, :), &
                      blk_sizes(1,:), blks_tuple_triang_size(1), & 
                      (/ translated_index(:) /))

       end if

       prop(pr_offset) = prop(pr_offset) + prop_forcache(ca_offset)

    end do

!     write(*,*) 'scfe contribution'
!     call print_rsp_tensor_stdout_tr(1, total_num_perturbations, merged_p_tuple%pdim, &
!     (/ (1, j = 1, (merged_p_tuple%n_perturbations - 1) ) /), merged_nblks, blk_sizes_merged, &
!     merged_blk_info, property_size, prop_forcache)

    call property_cache_add_element(cache, 2, p12,  &
         inner_indices_size * outer_indices_size, prop_forcache)    

    deallocate(deriv_structa)
    deallocate(deriv_structb)
    deallocate(ncarray)
    deallocate(ncinner)
    deallocate(outer_indices_a)
    deallocate(outer_indices_b)
    deallocate(outer_ind_a_large)
    deallocate(outer_ind_b_large)
    deallocate(which_index_is_pid1)
    deallocate(which_index_is_pid2)
    deallocate(nfields)
    deallocate(nblks_tuple)
    deallocate(blks_tuple_info)
    deallocate(blks_tuple_triang_size)
    deallocate(blk_sizes)
    deallocate(blk_sizes_merged)
    deallocate(prop_forcache)
    L = 0
    Y = 0

  end subroutine

  recursive subroutine print_rsp_tensor(npert, lvl, pdim, prop, offset)

    implicit none

    integer :: npert, i, j, offset, lvl, new_offset
    integer, dimension(npert) :: pdim
    complex(8), dimension(product(pdim)) :: prop

    if (lvl > 1) then

    do i = 1, pdim(npert - lvl + 1)

       new_offset = offset + (i - 1)*product(pdim(npert - lvl + 1:npert))/ &
                                             pdim(npert - lvl + 1)

       call print_rsp_tensor(npert, lvl - 1, pdim, prop, new_offset)

    end do

    open(unit=260, file='rsp_tensor', status='old', action='write', &
         position='append') 
    write(260,*) ' '
    close(260)

    else

    open(unit=260, file='rsp_tensor', status='old', action='write', &
         position='append') 
    write(260,*) real(prop(offset:offset+pdim(npert) - 1))
    close(260)

    end if

  end subroutine


  recursive subroutine print_rsp_tensor_tr(lvl, npert, pdim, ind, &
                       nblks, blk_sizes, blk_info, propsize, prop, print_id)

    implicit none

    integer :: lvl, npert, nblks, i, propsize, print_id
    integer, dimension(npert) :: pdim
    integer, dimension(npert - 1) :: ind, new_ind
    integer, dimension(nblks) :: blk_sizes
    integer, dimension(nblks, 3) :: blk_info
    complex(8), dimension(pdim(npert)) :: line_for_print
    complex(8), dimension(propsize) :: prop

    if (lvl == npert) then

       do i = 1, pdim(npert)

          line_for_print(i) = prop(get_triang_blks_offset(nblks, npert, &
                              blk_info, blk_sizes, (/ ind(:), i  /)))

       end do

!        open(unit=260, file='rsp_tensor', status='old', action='write', &
!             position='append') 
       write(print_id,*) real(line_for_print)
!        close(260)

    else

       new_ind = ind

       do i = 1, pdim(lvl)

          new_ind = ind
          new_ind(lvl) = i

          call print_rsp_tensor_tr(lvl + 1, npert, pdim, new_ind, &
                                   nblks, blk_sizes, blk_info, propsize, prop, print_id)

       end do

!        open(unit=260, file='rsp_tensor', status='old', action='write', &
!             position='append') 
       write(print_id,*) ' '
!        close(260)

    end if

  end subroutine


  recursive subroutine print_rsp_tensor_stdout_tr(lvl, npert, pdim, ind, &
                       nblks, blk_sizes, blk_info, propsize, prop)

    implicit none

    integer :: lvl, npert, nblks, i, propsize
    integer, dimension(npert) :: pdim
    integer, dimension(npert - 1) :: ind, new_ind
    integer, dimension(nblks) :: blk_sizes
    integer, dimension(nblks, 3) :: blk_info
    complex(8), dimension(pdim(npert)) :: line_for_print
    complex(8), dimension(propsize) :: prop

    if (lvl == npert) then

       do i = 1, pdim(npert)

          line_for_print(i) = prop(get_triang_blks_offset(nblks, npert, &
                                   blk_info, blk_sizes, (/ ind(:), i  /)))

       end do

       write(*,*) real(line_for_print)

    else

       new_ind = ind

       do i = 1, pdim(lvl)

          new_ind = ind
          new_ind(lvl) = i

          call print_rsp_tensor_stdout_tr(lvl + 1, npert, pdim, new_ind, &
                                          nblks, blk_sizes, blk_info, propsize, prop)

       end do

       write(*,*) ' '

    end if

  end subroutine


  recursive subroutine print_rsp_tensor_stdout(npert, lvl, pdim, prop, offset)

    implicit none

    integer :: npert, i, j, offset, lvl, new_offset
    integer, dimension(npert) :: pdim
    complex(8), dimension(product(pdim)) :: prop

    if (lvl > 1) then

       do i = 1, pdim(npert - lvl + 1)

          new_offset = offset + (i - 1)*product(pdim(npert - lvl + 1:npert))/ &
                                                pdim(npert - lvl + 1)

          call print_rsp_tensor_stdout(npert, lvl - 1, pdim, prop, new_offset)

       end do

       write(*,*) ' '

    else

       write(*,*) real(prop(offset:offset+pdim(npert) - 1))

    end if

  end subroutine


end module
