! Copyright 2012 Magnus Ringholm
!
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

!> @file Contains module rsp_general

!> General response routines. This module organizes, computes and prints
!> response function tensors.
module rsp_general

  use matrix_defop
  use rsp_contribs
!   use rsp_equations

  use rsp_field_tuple
  use rsp_indices_and_addressing
  use rsp_perturbed_matrices
  use rsp_perturbed_sdf
  use rsp_property_caching
  use rsp_sdf_caching

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

!   type rsp_cfg
! 
!      type(matrix) :: zeromat
! 
!   end type



  contains



  subroutine rsp_prop(pert_unordered, kn, F_unperturbed, D_unperturbed, S_unperturbed)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert, pert_unordered
    type(matrix) :: F_unperturbed, D_unperturbed, S_unperturbed
    type(SDF), pointer :: F, D, S
    integer, dimension(2) :: kn
    integer :: i, j, num_blks, property_size
    integer, allocatable, dimension(:,:) :: blk_info
    complex(8), allocatable, dimension(:) :: prop

    open(unit=257, file='totterms', status='replace', action='write') 
    write(257,*) 'START'
    close(257)

    open(unit=257, file='cachehit', status='replace', action='write') 
    write(257,*) 'START'
    close(257)

    pert = p_tuple_standardorder(pert_unordered)

!     call test_making_triangulated_indices(pert)

! The get_bestkn function is taken out of use for now until it has been improved
! The choice of k, n will likely be made before this routine is called anyway
! The call (commented out) is kept here for now
!    kn = get_bestkn(pert)
   write(*,*) ' '
   write(*,*) 'Choice of k, n is ', kn(1), ' and ', kn(2)
   write(*,*) ' '

    allocate(S)
    S%next => S
    S%last = .TRUE.
    S%perturb%n_perturbations = 0
    allocate(S%perturb%pdim(0))
    allocate(S%perturb%plab(0))
    allocate(S%perturb%pid(0))
    allocate(S%perturb%freq(0))
    allocate(S%data(1))
    S%data = S_unperturbed

    allocate(D)
    D%next => D
    D%last = .TRUE.
    D%perturb%n_perturbations = 0
    allocate(D%perturb%pdim(0))
    allocate(D%perturb%plab(0))
    allocate(D%perturb%pid(0))
    allocate(D%perturb%freq(0))
    allocate(D%data(1))
    D%data = D_unperturbed

    allocate(F)
    F%next => F
    F%last = .TRUE.
    F%perturb%n_perturbations = 0
    allocate(F%perturb%pdim(0))
    allocate(F%perturb%plab(0))
    allocate(F%perturb%pid(0))
    allocate(F%perturb%freq(0))
    allocate(F%data(1))
    F%data = F_unperturbed


    num_blks = get_num_blks(pert)


    allocate(blk_info(num_blks, 3))

    blk_info = get_blk_info(num_blks, pert)

!     write(*,*) 'blk_info is', blk_info

    property_size = get_triangulated_size(num_blks, blk_info)


    allocate(prop(product(pert%pdim)))
    prop = 0.0

    call get_prop(mol, pert, kn, property_size, prop, F, D, S)

    write(*,*) 'Property was calculated and printed to rsp_tensor'
    write(*,*) ' '
    open(unit=260, file='rsp_tensor', status='replace', action='write') 
    write(260,*) ' '
    close(260)

    call print_rsp_tensor(size(pert%pdim),size(pert%pdim),pert%pdim, prop, 1)

    write(*,*) 'End of print'

    open(unit=257, file='totterms', status='old', action='write', position='append') 
    write(257,*) 'END'
    close(257)

    open(unit=257, file='cachehit', status='old', action='write', position='append') 
    write(257,*) 'END'
    close(257)

    deallocate(blk_info)

  end subroutine



  subroutine get_prop(mol, pert, kn, property_size, prop, F, D, S)

    implicit none

    type(SDF) :: F, D, S
    type(rsp_cfg) :: mol
    type(p_tuple) :: pert, emptypert
    type(p_tuple), dimension(2) :: emptyp_tuples
    integer :: property_size
    integer, dimension(2) :: kn
    complex(8), dimension(property_size) :: prop
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

    call rsp_fds(mol, pert, kn, F, D, S)
 
    write(*,*) ' '
    write(*,*) 'Finished calculation of perturbed overlap/density/Fock matrices'
    write(*,*) ' '

    call property_cache_allocate(energy_cache)
    call rsp_energy(mol, pert, pert%n_perturbations, kn, 1, (/emptypert/), 0, D, &
                  property_size, energy_cache, prop)

    write(*,*) ' '
    write(*,*) 'Finished calculating energy-type contributions'
    write(*,*) ' '

    deallocate(energy_cache)

!     call property_cache_allocate(pulay_kn_cache)
!     call rsp_pulay_kn(mol, pert, kn, (/emptypert, emptypert/), S, D, F, &
!                       property_size, pulay_kn_cache, prop)
! 
!     write(*,*) ' '
!     write(*,*) 'Finished calculating Pulay k-n type contributions'
!     write(*,*) ' '
! 
!     deallocate(pulay_kn_cache)
! 
!     call property_cache_allocate(pulay_lag_cache)
!     call rsp_pulay_lag(mol, p_tuple_remove_first(pert), kn, &
!                        (/p_tuple_getone(pert,1), emptypert/), &
!                        S, D, F, property_size, pulay_lag_cache, prop)
! 
!     write(*,*) ' '
!     write(*,*) 'Finished calculating Pulay lagrangian type contributions' 
!     write(*,*) ' '
! 
!     deallocate(pulay_lag_cache)
! 
!     call property_cache_allocate(idem_cache)
!     call rsp_idem_lag(mol, p_tuple_remove_first(pert), kn, &
!                       (/p_tuple_getone(pert,1), emptypert/), &
!                       S, D, F, property_size, idem_cache, prop)
! 
!     write(*,*) ' '
!     write(*,*) 'Finished calculating idempotency lagrangian type contributions'
!     write(*,*) ' '
! 
!     deallocate(idem_cache)
! 
!     call property_cache_allocate(scfe_cache)
!     call rsp_scfe_lag(mol, p_tuple_remove_first(pert), kn, &
!                       (/p_tuple_getone(pert,1), emptypert/), &
!                       S, D, F, property_size, scfe_cache, prop)
! 
!     write(*,*) ' '
!     write(*,*) 'Finished calculating SCF lagrangian type contributions'
!     write(*,*) ' '
! 
!     deallocate(scfe_cache)

  end subroutine


  ! Calculate and add all the energy contributions

  recursive subroutine rsp_energy(mol, pert, total_num_perturbations, kn, num_p_tuples, &
                                p_tuples, density_order, D, property_size, cache, prop)

    implicit none

    logical :: e_knskip
    type(rsp_cfg) :: mol
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

       call rsp_energy(mol, p_tuple_remove_first(pert), total_num_perturbations, &
       kn, num_p_tuples, (/p_tuple_getone(pert,1), p_tuples(2:size(p_tuples))/), &
       density_order, D, property_size, cache, prop)

    else


       call rsp_energy(mol, p_tuple_remove_first(pert), total_num_perturbations,  &
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

          call rsp_energy(mol, p_tuple_remove_first(pert), total_num_perturbations, &
          kn, num_p_tuples, t_new, density_order + 1, D, property_size, cache, prop)

       end do


       ! 3. Chain rule differentiate the energy w.r.t. the density (giving 
       ! a(nother) pert D contraction)

       call rsp_energy(mol, p_tuple_remove_first(pert), total_num_perturbations, &
       kn, num_p_tuples + 1, (/p_tuples(:), p_tuple_getone(pert, 1)/), &
       density_order + 1, D, property_size, cache, prop)


    ! At the final recursion level: Calculate the contribution (if k,n choice of rule
    ! allows it) or get it from cache if it was already calculated (and if k,n choice 
    ! of rule allow it)

    else


    e_knskip = .FALSE.

       write(*,*) 'Getting energy contribution'

       do i = 1, num_p_tuples
 
          if (i > 1) then

             write(*,*) 'D ', p_tuples(i)%pid

             if(kn_skip(p_tuples(i)%n_perturbations, p_tuples(i)%pid, kn) .EQV. .TRUE.) then

                e_knskip = .TRUE.

             end if
          
          elseif (i == 1) then

             write(*,*) 'E ', p_tuples(i)%pid

          end if


       end do

       if (e_knskip .EQV. .FALSE.) then

          open(unit=257, file='totterms', status='old', action='write', &
               position='append') 
          write(257,*) 'T'
          close(257)
          
          write(*,*) 'Evaluating property_cache_already'

          if (property_cache_already(cache, num_p_tuples, p_tuples) .EQV. .TRUE.) then

             open(unit=257, file='cachehit', status='old', action='write', &
                  position='append') 
             write(257,*) 'T'
             close(257)

             write(*,*) 'Getting values from cache'


             ! NOTE (MaR): EVERYTHING IS IN STANDARD ORDER IN 
             ! THIS CALL (LIKE property_cache_getdata ASSUMES)
             call property_cache_getdata(cache, num_p_tuples, &
                  p_tuples_standardorder(num_p_tuples, p_tuples), property_size, prop)

             write(*,*) ' '
       
          else

             call get_energy(mol, num_p_tuples, total_num_perturbations, & 
                  (/ (p_tuple_standardorder(p_tuples(i)) , i = 1, num_p_tuples ) /), &
                  density_order, D, property_size, cache, prop)

                  write(*,*) 'Calculated energy contribution'
                  write(*,*) ' '

          end if

       else

          write(*,*) 'Energy contribution was k-n skipped'
          write(*,*) ' '

       end if

    end if

  end subroutine


! MR: WORKING ON IMPLEMENTING TRIANGULATED INDICES AND STORAGE HERE

  subroutine get_energy(mol, num_p_tuples, total_num_perturbations, &
                        p_tuples, density_order, D, property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple), dimension(num_p_tuples) :: p_tuples
    type(p_tuple) :: merged_p_tuple
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
    integer, dimension(total_num_perturbations) :: ncarray, ncouter, ncinner, pidouter
    integer, allocatable, dimension(:) :: o_whichpert, o_whichpertbig, o_wh_forave
    integer, allocatable, dimension(:) :: inner_offsets, ncoutersmall, pidoutersmall
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    complex(8), allocatable, dimension(:) :: tmp, contrib, prop_forcache
    complex(8), dimension(property_size) :: prop
!     complex(8), dimension(property_size) :: prop_forcache



    ncarray = get_ncarray(total_num_perturbations, num_p_tuples, p_tuples)
    ncouter = nc_only(total_num_perturbations, total_num_perturbations - &
              p_tuples(1)%n_perturbations, num_p_tuples - 1, &
              p_tuples(2:num_p_tuples), ncarray)
    ncinner = nc_only(total_num_perturbations, p_tuples(1)%n_perturbations, 1, &
                      p_tuples(1), ncarray)

    allocate(dens_tuple(num_p_tuples))
    allocate(nucpot_pert(p_tuples(1)%n_perturbations))
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



allocate(nfields(num_p_tuples))
allocate(nblks_tuple(num_p_tuples))




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
    



!     if (p_tuples(1)%n_perturbations > 0) then

       allocate(tmp(inner_indices_size))
       allocate(contrib(inner_indices_size))
       allocate(prop_forcache(inner_indices_size*outer_indices_size))




!     else

!        allocate(tmp(1))
!        allocate(contrib(1))
! 
!     end if
    prop_forcache = 0.0
    contrib = 0.0

    call sortdimbypid(total_num_perturbations, total_num_perturbations - &
                      p_tuples(1)%n_perturbations, pidoutersmall, &
                      ncarray, ncoutersmall, o_whichpert)


    if (total_num_perturbations > p_tuples(1)%n_perturbations) then


    allocate(outer_indices(outer_indices_size,size(ncoutersmall)))
    allocate(inner_indices(inner_indices_size,size(ncinnersmall)))

    do i = 1, size(o_whichpert)

       if (.NOT.(o_whichpert(i) == 0)) then

       o_wh_forave(o_whichpert(i)) = i

       end if
  
    end do

    k = 1

    do i = 2, num_p_tuples

       do j = 1, p_tuples(i)%n_perturbations

          ncoutersmall(k) =  p_tuples(i)%pdim(j)
          k = k + 1

       end do

    end do

    do i = 1, num_p_tuples

       dens_tuple(i) = mat_alloc_like(mol%zeromat)
       dens_tuple(i) = mat_zero_like(mol%zeromat)
       call mat_ensure_alloc(dens_tuple(i))

    end do

       call make_triangulated_tuples_indices(num_p_tuples - 1, total_num_perturbations, & 
            nblks_tuple(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, &
            :, :), blks_tuple_triang_size(2:num_p_tuples), outer_indices)


! write(*,*) 'oind:'
! 
! do i = 1, outer_indices_size
! 
! write(*,*) outer_indices(i,:)
! 
! end do


       if (p_tuples(1)%n_perturbations > 0) then

          call make_triangulated_indices(nblks_tuple(1), blks_tuple_info(1, &
               1:nblks_tuple(1), :), blks_tuple_triang_size(1), inner_indices)

       end if

    do i = 1, size(outer_indices, 1)

       dtup_ind = 0

       do j = 2, num_p_tuples

          call sdf_getdata_s(D, p_tuples(j), outer_indices(i, &
               dtup_ind+1:dtup_ind + p_tuples(j)%n_perturbations), dens_tuple(j))

! if (p_tuples(j)%n_perturbations == 3) then
! 
! write(*,*) 'looking at found density matrix at indices', outer_indices(i, &
!                dtup_ind+1:dtup_ind + p_tuples(j)%n_perturbations )
! 
! write(*,*) dens_tuple(j)%elms
! 
! end if

          dtup_ind = dtup_ind + p_tuples(j)%n_perturbations




       end do

       tmp = 0.0
       contrib = 0.0

       if (num_p_tuples == 1) then

          call rsp_oneave_tr(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                         (/ (1, j = 1, p_tuples(1)%n_perturbations) /), & 
                         p_tuples(1)%pdim, sdf_getdata(D, get_emptypert(), (/1/)), &
                         inner_indices_size, contrib)

       elseif (num_p_tuples == 2) then

          call rsp_oneave_tr(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                         (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                         p_tuples(1)%pdim, dens_tuple(2), inner_indices_size, contrib)

       end if

       tmp = tmp + contrib

! if (i == 1) then
! write(*,*) 'AFTER ONEAVE'
! write(*,*) real(tmp)
! write(*,*) ' '
! end if

       contrib = 0.0

       if (num_p_tuples == 1) then

          call rsp_twoave_tr(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                          (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                          p_tuples(1)%pdim, sdf_getdata(D, get_emptypert(), &
                          (/1/)), sdf_getdata(D, get_emptypert(), (/1/)) , &
                          inner_indices_size, contrib)

       elseif (num_p_tuples == 2) then

          call rsp_twoave_tr(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                          (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                          p_tuples(1)%pdim, dens_tuple(2), &
                          sdf_getdata(D, get_emptypert(), (/1/)) , &
                          inner_indices_size, contrib)

       elseif (num_p_tuples == 3) then

          call rsp_twoave_tr(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                          (/ (1, j = 1, p_tuples(1)%n_perturbations) /), &
                          p_tuples(1)%pdim, dens_tuple(2), dens_tuple(3), &
                          inner_indices_size, contrib)

       end if

       tmp = tmp + contrib

! if (i == 1) then
! write(*,*) 'AFTER TWOAVE'
! write(*,*) real(tmp)
! write(*,*) ' '
! end if

! NOTE (MaR): XCAVE CALL REMOVED FOR NOW
! 
!        contrib = 0.0
! 
!        call rsp_xcave(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
!                      (/ (1, j = 1, p_tuples(1)%n_perturbations) /), p_tuples(1)%pdim, &
!                      num_p_tuples, (/ sdf_getdata(D, get_emptypert(), (/1/)), &
!                      (dens_tuple(k), k = 2, num_p_tuples) /), contrib)
! 
! 
!        tmp = tmp + contrib

       if (p_tuples(1)%n_perturbations > 0) then

          do j = 1, size(inner_indices, 1)


! write(*,*) 'indices:', (/inner_indices(j, :), outer_indices(i, :) /)

offset = get_triang_blks_tuple_offset(num_p_tuples, total_num_perturbations, nblks_tuple, &
         (/ (p_tuples(k)%n_perturbations, k = 1, num_p_tuples) /), &
         blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
         (/inner_indices(j, :), outer_indices(i, :) /)) 

!              prop(offset) = prop(offset) + tmp(j)
             prop_forcache(offset) = prop_forcache(offset) + tmp(j)

          end do

       else

offset = get_triang_blks_tuple_offset(num_p_tuples - 1, total_num_perturbations,  &
         nblks_tuple(2:num_p_tuples), &
         (/ (p_tuples(k)%n_perturbations, k = 2, num_p_tuples) /), &
         blks_tuple_info(2:num_p_tuples, :, :), blk_sizes(2:num_p_tuples,:), & 
         blks_tuple_triang_size(2:num_p_tuples), &
         (/outer_indices(i, :) /)) 

!           prop(offset) = prop(offset) + tmp(1)
          prop_forcache(offset) = prop_forcache(offset) + tmp(1)

       end if

       end do





if (p_tuples(1)%n_perturbations > 0) then

call p_tuple_p1_cloneto_p2(p_tuples(1), merged_p_tuple)

do i = 2, num_p_tuples

! This can be problematic - consider rewriting merge_p_tuple as subroutine

merged_p_tuple = merge_p_tuple(merged_p_tuple, p_tuples(i))

end do
! write(*,*) '7'
else

call p_tuple_p1_cloneto_p2(p_tuples(2), merged_p_tuple)

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

allocate(triang_indices_pr(merged_triang_size, sum(merged_blk_info(1, :,2))))

! write(*,*) 'made allocation'

call make_triangulated_indices(merged_nblks, merged_blk_info, & 
     merged_triang_size, triang_indices_pr)

! write(*,*) '8b'
do i = 1, size(triang_indices_pr, 1)

pr_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
            (/sum(nfields)/), &
            (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
            (/triang_indices_pr(i, :) /))
! write(*,*) 'indices are', triang_indices_pr(i, :)

! write(*,*) 'pr offset is', pr_offset

if (p_tuples(1)%n_perturbations > 0) then

ec_offset = get_triang_blks_tuple_offset(num_p_tuples, &
            total_num_perturbations, nblks_tuple, &
            nfields, blks_tuple_info, blk_sizes, blks_tuple_triang_size, &
            (/triang_indices_pr(i, :) /))


else

ec_offset = get_triang_blks_tuple_offset(num_p_tuples - 1, &
            total_num_perturbations, nblks_tuple(2:num_p_tuples), &
            nfields(2:num_p_tuples), blks_tuple_info(2:num_p_tuples, :, :), &
            blk_sizes(2:num_p_tuples,:), blks_tuple_triang_size(2:num_p_tuples), & 
            (/triang_indices_pr(i, :) /))



end if


! write(*,*) 'lo offset is', lo_offset 
! write(*,*) '8d'
prop(pr_offset) = prop(pr_offset) + prop_forcache(ec_offset)
! write(*,*) '8e'
end do






    call property_cache_add_element(cache, num_p_tuples, p_tuples, &
                                    inner_indices_size * outer_indices_size, prop_forcache)   


deallocate(merged_blk_info)
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

!        write(*,*) 'all indices inner'

       tmp = 0.0
       contrib = 0.0

       call rsp_nucpot_tr(nucpot_pert, property_size, contrib) 
       tmp = tmp + contrib

! write(*,*) 'AFTER NUCPOT'
! write(*,*) real(contrib(1))
! write(*,*) ' '

       contrib = 0.0

       call rsp_oneave_tr(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                       (/ (1, j = 1, p_tuples(1)%n_perturbations) /), p_tuples(1)%pdim, &
                       sdf_getdata(D, get_emptypert(), (/1/)) , property_size, contrib)

       tmp = tmp + contrib

!  write(*,*) 'AFTER ONEAVE'

       contrib = 0.0

       call rsp_twoave_tr(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
                       (/ (1, j = 1, p_tuples(1)%n_perturbations) /), p_tuples(1)%pdim, &
                       sdf_getdata(D, get_emptypert(), (/1/)) , &
                       sdf_getdata(D, get_emptypert(), (/1/)) , property_size, contrib)

       tmp = tmp + 0.5*(contrib)

! write(*,*) 'AFTER TWOAVE'

! NOTE (MaR): XCAVE CALL REMOVED FOR NOW
 
!        contrib = 0.0
!
!        call rsp_xcave(mol, p_tuples(1)%n_perturbations, p_tuples(1)%plab, &
!                      (/ (1, j = 1, p_tuples(1)%n_perturbations) /), p_tuples(1)%pdim, &
!                      1, (/ sdf_getdata(D, get_emptypert(), (/1/)) /), contrib)
! 
!        tmp = tmp + contrib

       prop =  prop + tmp
!        prop_forcache = prop_forcache + tmp

    end if

 

!  write(*,*) 'energy contribution'
!  call print_rsp_tensor_stdout(total_num_perturbations,total_num_perturbations, &
!                               ncarray, prop_forcache, 1)


deallocate(nfields)
deallocate(nblks_tuple)
deallocate(blks_tuple_info)
deallocate(blks_tuple_triang_size)


deallocate(blk_sizes)
deallocate(blk_sizes_merged)

    deallocate(nucpot_pert)
    deallocate(dens_tuple)
    deallocate(ncoutersmall)
    deallocate(ncinnersmall)
    deallocate(pidoutersmall)
    deallocate(o_whichpert)
    deallocate(o_wh_forave)
!     deallocate(outer_indices)
!     deallocate(inner_indices)
    deallocate(tmp)
    deallocate(contrib)

  end subroutine




  recursive subroutine rsp_pulay_kn(mol, pert, kn, p12, S, D, F, property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    integer :: property_size, i
    integer, dimension(2) :: kn
    complex(8), dimension(property_size) :: prop
    
    if (pert%n_perturbations > 0) then

       call rsp_pulay_kn(mol, p_tuple_remove_first(pert), kn, &
       (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), S, D, F, property_size, &
       cache, prop)

       call rsp_pulay_kn(mol, p_tuple_remove_first(pert), kn, &
       (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), S, D, F, property_size, &
       cache, prop)

    else

       if (kn_skip(p12(2)%n_perturbations, p12(2)%pid, kn) .EQV. .FALSE.) then


          write(*,*) 'Getting Pulay k-n contribution:'
          write(*,*) 'S', p12(1)%pid
          write(*,*) 'W', p12(2)%pid

          open(unit=257, file='totterms', status='old', action='write', &
               position='append')
          write(257,*) 'T'
          close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

             write(*,*) 'Getting values from cache'
             write(*,*) ' '

             open(unit=257, file='cachehit', status='old', action='write', &
                  position='append') 
             write(257,*) 'T'
             close(257)

             call property_cache_getdata(cache, 2, p12, property_size, prop)
       
          else

             call get_pulay_kn(mol, (/ (p_tuple_standardorder(p12(i)) , i = 1, 2)  /), & 
                               kn, F, D, S, property_size, cache, prop)

             write(*,*) 'Calculated Pulay k-n contribution'
             write(*,*) ' '

          end if

       else

          write(*,*) 'Pulay k-n contribution was k-n skipped:'
          write(*,*) 'S ', p12(1)%pid 
          write(*,*) 'W ', p12(2)%pid 
          write(*,*) ' '

       end if 

    end if

  end subroutine


  subroutine get_pulay_kn(mol, p12, kn, F, D, S, property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert, emptypert
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structb
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    type(matrix) :: W
    integer :: i, j, sstr_incr, offset, &
             property_size, dtup_ind, pr_offset, ec_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, merged_nblks
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:) :: ncinnersmall, blk_sizes_merged
    integer, allocatable, dimension(:,:) :: triang_indices_pr, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    integer :: d_supsize
    integer, dimension(2) :: kn
    integer, allocatable, dimension(:) :: ncarray, ncinner, inner_offsets, &
                                          which_index_is_pid
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    complex(8), allocatable, dimension(:) :: tmp
    complex(8), dimension(property_size) :: prop
    complex(8), dimension(property_size) :: prop_forcache

    prop_forcache = 0.0

    d_supsize = derivative_superstructure_getsize(mol, p12(2), kn, .FALSE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/))

    allocate(deriv_structb(d_supsize, 3))

    sstr_incr = 0

    call derivative_superstructure(mol, p12(2), kn, .FALSE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize, sstr_incr, deriv_structb)

    allocate(ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(ncinner(p12(1)%n_perturbations))
    allocate(tmp(product(p12(1)%pdim)))
    allocate(inner_offsets(product(p12(1)%pdim)))
    allocate(outer_indices(product(p12(2)%pdim), p12(2)%n_perturbations))
    allocate(inner_indices(product(p12(1)%pdim), p12(1)%n_perturbations))
    allocate(which_index_is_pid(p12(1)%n_perturbations + p12(2)%n_perturbations))

    ncarray = get_ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations, 2, p12)
    ncinner = nc_onlysmall(p12(1)%n_perturbations + p12(2)%n_perturbations, &
                      p12(1)%n_perturbations, 1, p12(1), ncarray)

    which_index_is_pid = 0

    do i = 1, p12(2)%n_perturbations

       which_index_is_pid(p12(2)%pid(i)) = i

    end do

    call make_indices(p12(1)%n_perturbations, 1, p12(1)%pdim, 0, inner_indices)
    call make_indices(p12(2)%n_perturbations, 1, p12(2)%pdim, 0, outer_indices)

    W = mat_alloc_like(mol%zeromat)
    W = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(W)

    do i = 1, size(outer_indices, 1)

       tmp = 0.0
       W = mat_zero_like(mol%zeromat)

! call sdf_getdata_s(D, deriv_structb(1,1), get_fds_data_index(deriv_structb(1,1), &
!      p12(1)%n_perturbations + p12(2)%n_perturbations, which_index_is_pid, &
!      p12(2)%n_perturbations, outer_indices(i,:)), W)
! 
! write(*,*) 'got D test at i = ', i
! write(*,*) W%elms
! 
! call sdf_getdata_s(F, deriv_structb(1,1), get_fds_data_index(deriv_structb(1,1), &
!      p12(1)%n_perturbations + p12(2)%n_perturbations, which_index_is_pid, &
!      p12(2)%n_perturbations, outer_indices(i,:)), W)
! 
! write(*,*) 'got F test at i = ', i
! write(*,*) W%elms



       call rsp_get_matrix_w(mol, d_supsize, deriv_structb, p12(1)%n_perturbations + &
                            p12(2)%n_perturbations, which_index_is_pid, &
                            p12(2)%n_perturbations, outer_indices(i,:), F, D, S, W)

! write(*,*) 'got W at i = ', i
! write(*,*) W%elms

!MR: TEMPORARILY DISABLED
!        call rsp_ovlave(mol, p12(1)%n_perturbations, p12(1)%plab, &
!                       (/ (j/j, j = 1, p12(1)%n_perturbations) /), p12(1)%pdim, W, tmp)

       do j = 1, size(inner_indices, 1)

          offset = get_one_tensor_offset(p12(1)%n_perturbations + &
                   p12(2)%n_perturbations, (/inner_indices(j,:), &
                   outer_indices(i,:) /), (/ p12(1)%pid(:), p12(2)%pid(:) /), ncarray)

          prop(offset) = prop(offset) + tmp(j)
          prop_forcache(offset) = prop_forcache(offset) + tmp(j)

       end do

    end do

!     write(*,*) 'pulay kn contribution'
! 
!  call print_rsp_tensor_stdout(p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               ncarray, prop_forcache, 1)


    call property_cache_add_element(cache, 2, p12, property_size, prop_forcache)    

    deallocate(deriv_structb)
    deallocate(ncarray)
    deallocate(ncinner)
    deallocate(tmp)
    deallocate(inner_offsets)
    deallocate(outer_indices)
    deallocate(inner_indices)
    deallocate(which_index_is_pid)
    W = 0

  end subroutine




  recursive subroutine rsp_pulay_lag(mol, pert, kn, p12, S, D, F, property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    type(p_tuple), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    integer :: property_size, i
    integer, dimension(2) :: kn
    complex(8), dimension(property_size) :: prop
    
    if (pert%n_perturbations > 0) then

       call rsp_pulay_lag(mol, p_tuple_remove_first(pert), kn, &
       (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), &
       S, D, F, property_size, cache, prop)
       call rsp_pulay_lag(mol, p_tuple_remove_first(pert), kn, &
       (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), &
       S, D, F, property_size, cache, prop)

    else

       ! At lowest level:
       if (kn_skip(p12(1)%n_perturbations, p12(1)%pid, kn) .EQV. .FALSE.) then

       write(*,*) 'Getting Pulay lagrange contribution:'
       write(*,*) 'S', p12(1)%pid
       write(*,*) 'W', p12(2)%pid, 'primed', kn(2)

       open(unit=257, file='totterms', status='old', action='write', position='append') 
       write(257,*) 'T'
       close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

             write(*,*) 'Getting values from cache'
             write(*,*) ' '
       
             open(unit=257, file='cachehit', status='old', action='write', &
             position='append') 
             write(257,*) 'T'
             close(257)

             call property_cache_getdata(cache, 2, p12, property_size, prop)

          else

             call get_pulay_lag(mol, (/ (p_tuple_standardorder(p12(i)) , i = 1, 2) /), & 
                               kn, F, D, S, property_size, cache, prop)

             write(*,*) 'Calculated Pulay lagrange contribution'
             write(*,*) ' '

          end if

       else

          write(*,*) 'Pulay lagrange contribution was k-n skipped:'
          write(*,*) 'S', p12(1)%pid 
          write(*,*) 'W', p12(2)%pid, 'primed', kn(2)
          write(*,*) ' '

       end if

    end if

  end subroutine


  subroutine get_pulay_lag(mol, p12, kn, F, D, S, property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert, emptypert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structb
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    type(matrix) :: W
    integer :: i, j, k ,m, incr
    integer :: d_supsize, &
             property_size, offset, dtup_ind, pr_offset, ec_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, merged_nblks
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size
    integer, allocatable, dimension(:) :: ncinnersmall, blk_sizes_merged
    integer, allocatable, dimension(:,:) :: triang_indices_pr, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    integer, dimension(2) :: kn
    integer, allocatable, dimension(:) :: ncarray, ncinner, inner_offsets, &
                                          which_index_is_pid
    integer, allocatable, dimension(:) :: outer_ind_b_large
    integer, allocatable, dimension(:,:) :: outer_indices, inner_indices
    complex(8), allocatable, dimension(:) :: tmp
    complex(8), dimension(property_size) :: prop
    complex(8), dimension(property_size) :: prop_forcache

    prop_forcache = 0.0

    d_supsize = derivative_superstructure_getsize(mol, p12(2), kn, .TRUE., &
                (/get_emptypert(), get_emptypert(), get_emptypert()/))
   
    allocate(deriv_structb(d_supsize, 3))

    incr = 0

    call derivative_superstructure(mol, p12(2), kn, .TRUE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize, incr, deriv_structb)

    allocate(ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(ncinner(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_ind_b_large(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(tmp(product(p12(1)%pdim)))
    allocate(inner_offsets(product(p12(1)%pdim)))
    allocate(outer_indices(product(p12(2)%pdim), p12(2)%n_perturbations))
    allocate(inner_indices(product(p12(1)%pdim), p12(1)%n_perturbations))
    allocate(which_index_is_pid(p12(1)%n_perturbations + p12(2)%n_perturbations))

    ncarray = get_ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations, 2, p12)
    ncinner = nc_only(p12(1)%n_perturbations + p12(2)%n_perturbations, &
              p12(1)%n_perturbations, 1, p12(1), ncarray)

    which_index_is_pid = 0

    do i = 1, p12(2)%n_perturbations

       which_index_is_pid(p12(2)%pid(i)) = i

    end do

    call make_indices(p12(2)%n_perturbations, 1, p12(2)%pdim, 0, outer_indices)
    call make_indices(p12(1)%n_perturbations, 1, p12(1)%pdim, 0, inner_indices)

    W = mat_alloc_like(mol%zeromat)
    W = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(W)

    do i = 1, size(outer_indices, 1)

       tmp = 0.0
       W = mat_zero_like(mol%zeromat)

       call rsp_get_matrix_w(mol, d_supsize, deriv_structb, p12(1)%n_perturbations + &
                            p12(2)%n_perturbations, which_index_is_pid, &
                            p12(2)%n_perturbations, outer_indices(i,:), F, D, S, W)

!MR: TEMPORARILY DISABLED
!        call rsp_ovlave(mol, p12(1)%n_perturbations, p12(1)%plab, &
!                        (/ (j/j, j = 1, p12(1)%n_perturbations) /), &
!                        p12(1)%pdim, W, tmp)

       do j = 1, size(inner_indices, 1)

          offset = get_one_tensor_offset(p12(1)%n_perturbations + &
                   p12(2)%n_perturbations, (/inner_indices(j,:), &
                   outer_indices(i,:) /), (/ p12(1)%pid(:), p12(2)%pid(:) /), ncarray)

          prop(offset) = prop(offset) + tmp(j)
          prop_forcache(offset) = prop_forcache(offset) + tmp(j)

       end do

    end do
!     write(*,*) 'pulay lag contribution'
! 
!  call print_rsp_tensor_stdout(p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               ncarray, prop_forcache, 1)

    call property_cache_add_element(cache, 2, p12, property_size, prop_forcache)

    deallocate(deriv_structb)
    deallocate(ncarray)
    deallocate(ncinner)
    deallocate(outer_ind_b_large)
    deallocate(tmp)
    deallocate(inner_offsets)
    deallocate(outer_indices)
    deallocate(inner_indices)
    deallocate(which_index_is_pid)
    W = 0

  end subroutine




  recursive subroutine rsp_idem_lag(mol, pert, kn, p12, S, D, F, &
                                    property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    type(p_tuple), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    integer :: property_size, i
    integer, dimension(2) :: kn
    complex(8), dimension(property_size) :: prop
    
    if (pert%n_perturbations > 0) then

       call rsp_idem_lag(mol, p_tuple_remove_first(pert), kn, &
       (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), S, D, F, property_size, &
       cache, prop)
       call rsp_idem_lag(mol, p_tuple_remove_first(pert), kn, &
       (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), S, D, F, property_size, &
       cache, prop)

    else

       if (kn_skip(p12(1)%n_perturbations, p12(1)%pid, kn) .EQV. .FALSE.) then

          write(*,*) 'Getting idempotency lagrange contribution'
          write(*,*) 'Zeta', p12(1)%pid
          write(*,*) 'Z', p12(2)%pid, 'primed', kn(2)

          open(unit=257, file='totterms', status='old', action='write', &
               position='append') 
          write(257,*) 'T'
          close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

             write(*,*) 'Getting values from cache'
             write(*,*) ' '

             open(unit=257, file='cachehit', status='old', action='write', &
                  position='append')
             write(257,*) 'T'
             close(257)

             call property_cache_getdata(cache, 2, p12, property_size, prop)
      
          else

             ! At lowest level:
             call get_idem_lag(mol, (/ (p_tuple_standardorder(p12(i)) , i = 1, 2) /), & 
                               kn, F, D, S, property_size, cache, prop)

             write(*,*) 'Calculated idempotency lagrange contribution'
             write(*,*) ' '

          end if

       else

          write(*,*) 'Idempotency lagrange contribution was k-n skipped:'
          write(*,*) 'Zeta', p12(1)%pid 
          write(*,*) 'Z', p12(2)%pid, 'primed', kn(2)
          write(*,*) ' '

       end if

    end if

  end subroutine


  subroutine get_idem_lag(mol, p12, kn, F, D, S, property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert, emptypert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structa, deriv_structb
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    type(matrix) :: Zeta, Z
    integer :: i, j, k, m, n, p, incr1, incr2, &
             property_size, offset, dtup_ind, pr_offset, ec_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, merged_nblks
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
    complex(8), dimension(property_size) :: prop_forcache

    prop_forcache = 0.0

    d_supsize = 0

    d_supsize(1) = derivative_superstructure_getsize(mol, p_tuple_remove_first(p12(1)), &
                   kn, .FALSE., (/get_emptypert(), get_emptypert(), get_emptypert()/))
    d_supsize(2) = derivative_superstructure_getsize(mol, p12(2), &
                   kn, .TRUE., (/get_emptypert(), get_emptypert(), get_emptypert()/))

    allocate(deriv_structa(d_supsize(1), 3))
    allocate(deriv_structb(d_supsize(2), 3))

    incr1 = 0
    incr2 = 0

    call derivative_superstructure(mol, p_tuple_remove_first(p12(1)), kn, .FALSE., & 
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize(1), incr1, deriv_structa)
    call derivative_superstructure(mol, p12(2), kn, .TRUE., &
         (/get_emptypert(), get_emptypert(), get_emptypert()/), &
         d_supsize(2), incr2, deriv_structb)


    allocate(ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(ncinner(p12(1)%n_perturbations + p12(2)%n_perturbations))
!     allocate(ncprod(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_ind_a_large(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_ind_b_large(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_indices_a(product(p12(1)%pdim), p12(1)%n_perturbations))
    allocate(outer_indices_b(product(p12(2)%pdim), p12(2)%n_perturbations))
    allocate(which_index_is_pid1(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(which_index_is_pid2(p12(1)%n_perturbations + p12(2)%n_perturbations))

    ncarray = get_ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations, 2, p12)
    ncinner = nc_only(p12(1)%n_perturbations + p12(2)%n_perturbations, &
                      p12(1)%n_perturbations, 1, p12(1), ncarray)

!     do i = 1, size(ncarray)
! 
!        ncprod(i) = product(ncarray(i:size(ncarray)))/ncarray(i)
! 
!     end do

    which_index_is_pid1 = 0

    do i = 1, p12(1)%n_perturbations

       which_index_is_pid1(p12(1)%pid(i)) = i

    end do

    which_index_is_pid2 = 0

    do i = 1, p12(2)%n_perturbations

       which_index_is_pid2(p12(2)%pid(i)) = i

    end do

    call make_indices(p12(1)%n_perturbations, 1, p12(1)%pdim, 0, outer_indices_a)
    call make_indices(p12(2)%n_perturbations, 1, p12(2)%pdim, 0, outer_indices_b)

    offset = 0.0

    Z = mat_alloc_like(mol%zeromat)
    Z = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(Z)

    Zeta = mat_alloc_like(mol%zeromat)
    Zeta = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(Zeta)

    do i = 1, size(outer_indices_a, 1)

       Zeta = mat_zero_like(mol%zeromat)

       call rsp_get_matrix_zeta(mol, p_tuple_getone(p12(1), 1), kn, d_supsize(1), &
            deriv_structa, p12(1)%n_perturbations + p12(2)%n_perturbations, &
            which_index_is_pid1, p12(1)%n_perturbations, outer_indices_a(i,:), &
            F, D, S, Zeta)

       do j = 1, size(outer_indices_b, 1)

          Z = mat_zero_like(mol%zeromat)

          call rsp_get_matrix_z(mol, d_supsize(2), deriv_structb, kn, &
               p12(1)%n_perturbations + p12(2)%n_perturbations, which_index_is_pid2, &
               p12(2)%n_perturbations, outer_indices_b(j,:), F, D, S, Z)

          offset = get_one_tensor_offset(p12(1)%n_perturbations + &
                   p12(2)%n_perturbations, (/outer_indices_a(i,:), &
                   outer_indices_b(j,:) /), (/ p12(1)%pid(:), p12(2)%pid(:) /), ncarray)

          prop(offset) = prop(offset) -tr(Zeta, Z)
          prop_forcache(offset) = prop_forcache(offset) -tr(Zeta, Z)

       end do

    end do

!     write(*,*) 'idempotency contribution'
! 
!  call print_rsp_tensor_stdout(p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               ncarray, prop_forcache, 1)

    call property_cache_add_element(cache, 2, p12, property_size, prop_forcache) 

    deallocate(deriv_structa)
    deallocate(deriv_structb)
    deallocate(ncarray)
    deallocate(ncinner)
!     deallocate(ncprod)
    deallocate(outer_ind_a_large)
    deallocate(outer_ind_b_large)
    deallocate(outer_indices_a)
    deallocate(outer_indices_b)
    deallocate(which_index_is_pid1)
    deallocate(which_index_is_pid2)
    Zeta = 0
    Z = 0

  end subroutine



  recursive subroutine rsp_scfe_lag(mol, pert, kn, p12, S, D, F, &
                                    property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert
    type(p_tuple), dimension(2) :: p12
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    integer :: property_size, i
    integer, dimension(2) :: kn
    complex(8), dimension(property_size) :: prop
    
    if (pert%n_perturbations > 0) then

       call rsp_scfe_lag(mol, p_tuple_remove_first(pert), kn, &
            (/p_tuple_extend(p12(1), p_tuple_getone(pert, 1)), p12(2)/), &
            S, D, F, property_size, cache, prop)
       call rsp_scfe_lag(mol, p_tuple_remove_first(pert), kn, &
            (/p12(1), p_tuple_extend(p12(2), p_tuple_getone(pert, 1))/), &
            S, D, F, property_size, cache, prop)

    else

       if (kn_skip(p12(1)%n_perturbations, p12(1)%pid, kn) .EQV. .FALSE.) then

          write(*,*) 'Getting scfe lagrange contribution'
          write(*,*) 'Lambda', p12(1)%pid
          write(*,*) 'Y', p12(2)%pid, 'primed', kn(2)

          open(unit=257, file='totterms', status='old', action='write', &
               position='append') 
          write(257,*) 'T'
          close(257)

          if (property_cache_already(cache, 2, p12) .EQV. .TRUE.) then

             open(unit=257, file='cachehit', status='old', action='write', &
                  position='append') 
             write(257,*) 'T'
             close(257)

             write(*,*) 'Getting values from cache'
             write(*,*) ' '

             call property_cache_getdata(cache, 2, p12, property_size, prop)
       
          else

             ! At lowest level:
             call get_scfe_lag(mol, (/ (p_tuple_standardorder(p12(i)) , i = 1, 2) /), &
             kn, F, D, S, property_size, cache, prop)

             write(*,*) 'Calculated scfe lagrange contribution'
             write(*,*) ' '

          end if

       else

          write(*,*) 'scfe lagrange contribution was k-n skipped:'
          write(*,*) 'Lambda', p12(1)%pid 
          write(*,*) 'Y', p12(2)%pid, 'primed', kn(2)
          write(*,*) ' '

       end if

    end if

  end subroutine




  subroutine get_scfe_lag(mol, p12, kn, F, D, S, property_size, cache, prop)

    implicit none

    type(rsp_cfg) :: mol
    type(p_tuple) :: pert, emptypert, merged_p_tuple
    type(p_tuple), dimension(2) :: p12
    type(p_tuple), dimension(:,:), allocatable :: deriv_structa, deriv_structb
    type(SDF) :: S, D, F
    type(property_cache) :: cache
    type(matrix) :: L, Y
    integer :: i, j, k, m, n, p, incr1, incr2, &
             property_size, offset, dtup_ind, pr_offset, ec_offset, inner_indices_size, &
               outer_indices_size, merged_triang_size, merged_nblks
    integer, allocatable, dimension(:) :: ncinnersmall, blk_sizes_merged
    integer, allocatable, dimension(:,:) :: triang_indices_pr, blk_sizes
    integer, allocatable, dimension(:,:,:) :: merged_blk_info, blks_tuple_info
    integer, dimension(2) :: kn, d_supsize
    integer, allocatable, dimension(:) :: ncarray, ncinner, which_index_is_pid1, which_index_is_pid2
    integer, allocatable, dimension(:) :: outer_ind_a_large, outer_ind_b_large
    integer, allocatable, dimension(:,:) :: outer_indices_a, outer_indices_b
    complex(8), dimension(property_size) :: prop
    complex(8), dimension(property_size) :: prop_forcache

    prop_forcache = 0.0
    d_supsize = 0

    d_supsize(1) = derivative_superstructure_getsize(mol, p_tuple_remove_first(p12(1)), &
                   kn, .FALSE., (/get_emptypert(), get_emptypert(), get_emptypert()/))
    d_supsize(2) = derivative_superstructure_getsize(mol, p12(2), &
                   kn, .TRUE., (/get_emptypert(), get_emptypert(), get_emptypert()/))

    allocate(deriv_structa(d_supsize(1), 3))
    allocate(deriv_structb(d_supsize(2), 3))

    incr1 = 0
    incr2 = 0

    call derivative_superstructure(mol, p_tuple_remove_first(p12(1)), kn, .FALSE., &
                    (/get_emptypert(), get_emptypert(), get_emptypert()/), & 
                    d_supsize(1), incr1, deriv_structa)
    call derivative_superstructure(mol, p12(2), kn, .TRUE., &
                    (/get_emptypert(), get_emptypert(), get_emptypert()/), &
                    d_supsize(2), incr2, deriv_structb)

    allocate(ncarray(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(ncinner(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_indices_a(product(p12(1)%pdim), p12(1)%n_perturbations))
    allocate(outer_indices_b(product(p12(2)%pdim), p12(2)%n_perturbations))
    allocate(outer_ind_a_large(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(outer_ind_b_large(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(which_index_is_pid1(p12(1)%n_perturbations + p12(2)%n_perturbations))
    allocate(which_index_is_pid2(p12(1)%n_perturbations + p12(2)%n_perturbations))

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

    call make_indices(p12(1)%n_perturbations, 1, p12(1)%pdim, 0, outer_indices_a)
    call make_indices(p12(2)%n_perturbations, 1, p12(2)%pdim, 0, outer_indices_b)

    offset = 0

    Y = mat_alloc_like(mol%zeromat)
    Y = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(Y)

    L = mat_alloc_like(mol%zeromat)
    L = mat_zero_like(mol%zeromat)
    call mat_ensure_alloc(L)

    do i = 1, size(outer_indices_a, 1)

       L = mat_zero_like(mol%zeromat)

       call rsp_get_matrix_lambda(mol, p_tuple_getone(p12(1), 1), d_supsize(1), &
            deriv_structa, p12(1)%n_perturbations + p12(2)%n_perturbations, &
            which_index_is_pid1, p12(1)%n_perturbations, outer_indices_a(i,:), D, S, L)

       do j = 1, size(outer_indices_b, 1)

          Y = mat_zero_like(mol%zeromat)

          call rsp_get_matrix_y(mol, d_supsize(2), deriv_structb, &
               p12(1)%n_perturbations + p12(2)%n_perturbations, which_index_is_pid2, &
               p12(2)%n_perturbations, outer_indices_b(j,:), F, D, S, Y)

          offset = get_one_tensor_offset(p12(1)%n_perturbations + &
                   p12(2)%n_perturbations, (/outer_indices_a(i,:), &
          outer_indices_b(j,:) /), (/ p12(1)%pid(:), p12(2)%pid(:) /), ncarray)

          prop(offset) = prop(offset) - tr(L, Y)
          prop_forcache(offset) = prop_forcache(offset) - tr(L, Y)

       end do

    end do

!     write(*,*) 'scfe contribution'
!  
! call print_rsp_tensor_stdout(p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               p12(1)%n_perturbations + p12(2)%n_perturbations, &
!                               ncarray, prop_forcache, 1)

    call property_cache_add_element(cache, 2, p12, property_size, prop_forcache)

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
