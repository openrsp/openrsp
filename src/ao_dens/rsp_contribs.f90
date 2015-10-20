! Copyright 2012      Gao Bin
!           2012      Radovan Bast
!           2009-2012 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

!> @file Contains module rsp_contribs

!> This module contains routines for calculating contributions
!> to molecular properties (1st order, linear response, etc.),
!> and perturbed Fock matrices.
module rsp_contribs

!  use matrix_defop, matrix => openrsp_matrix
!  use matrix_lowlevel, only: mat_init, mat_zero_like
!  use interface_molecule, only: get_num_atoms
!  use interface_xc, only: rsp_xcint_interface
!  use interface_1el, only: interface_1el_ovlave_tr,        &
!                           interface_1el_ovlint_tr,        &
!                           interface_1el_oneave_tr,        &
!                           interface_1el_oneint_tr,        &
!                           interface_1el_ovlave_half_diff, &
!                           interface_1el_ovlint_half_diff
!  use interface_nuclear, only: nuclear_potential
!  use basis_set,  only: cgto
  use rsp_field_tuple, only: p_tuple, p_tuple_remove_first, p_tuple_getone, &
                             p_tuple_standardorder, merge_p_tuple, &
                             p_tuple_p1_cloneto_p2, p_tuple_to_external_tuple
  use rsp_indices_and_addressing
  use qcmatrix_f
                                        

  implicit none

!  public rsp_field
!  public rsp_field_bas
!  public rsp_field_dim
!  public rsp_field_anti
!  public rsp_field_ordering
!
!  ! MaR: ROUTINES RETURNING NON-REDUNDANT VALUES W.R.T TENSOR SYMMETRY
!  ! MaR: THESE COULD REPLACE THE CORRESPONDING ABOVE ROUTINES WHEN APPROPRIATE
!  public rsp_nucpot
!  public rsp_ovlave
  public rsp_ovlave_t_matrix_2014
!  public rsp_oneave
!  public rsp_ovlint
  public rsp_ovlint_t_matrix_2014
!  public rsp_oneint
!  public rsp_xcint_adapt
!  public rsp_pe



  private

contains


!  ! MaR: Seems to work properly, but memory usage is not tensor symmetry nonredundant
!  !> Contribution from nuclear repulsion and nuclei--field interaction
!  !> to response functions. Fields (type rsp_field) are here in arbitrary order.
!  !> (in normal mode) Fields are sorted, component ranges extended to what
!  !> rsp_backend expects (currently only full ranges). Then the call is then relayed
!  !> to rsp_backend's nuclear_potential, which computes and returns the requested real
!  !> tensor in standard order. The requested ranges of this tensor is then reordered
!  !> and added to rspfunc
!  subroutine rsp_nucpot(fields, propsize, rspfunc_output)
!    !> field descriptors (label freq comp ncomp)
!    type(rsp_field), intent(in)    :: fields(:)
!    !> output tensor, to which nuclear contribution is *ADDED*
!    integer                        :: propsize
!    complex(8),      intent(inout) :: rspfunc_output(propsize)
!    !> tmp tensor, to which nuclear contribution is *ADDED*
!    complex(8) :: rspfunc(product(fields%ncomp))
!    !---------------------------------------------------------------
!    integer      nf, ncor, ngeo, ext_ncomp, i
!    integer      order(size(fields)), tcomp(size(fields))
!    character(4) ext_label(2)
!    logical      nonz
!    rspfunc = 0.0
!    ! prepare and determine ordering. This also validates comp/ncomp
!    call count_and_prepare(fields, .true., .false., order, tcomp, nucpot = nonz)
!    ! early return if zero
!    if (.not.nonz) return
!    ! find number of GEO, and two, one or none external fields
!    nf = size(fields)
!    ngeo = count(fields%label == 'GEO ')
!    ext_label(:) = (/'NONE','NONE'/)
!    ext_ncomp = 1
!    if (ngeo == nf-1) then
!       ext_label(1) = fields(order(nf))%label
!       ext_ncomp    = tcomp(nf)
!    else if (ngeo == nf-2) then
!       ext_label(1) = fields(order(nf-1))%label
!       ext_label(2) = fields(order(nf))%label
!       ext_ncomp    = tcomp(nf-1) * tcomp(nf)
!    end if
!    ! use inner to avoid allocate'ing 'nucpot' below
!    ncor = 1
!    if (ngeo > 0) ncor = tcomp(1)
!    call inner
!  contains
!    subroutine inner
!      real(8) nucpot(ncor**ngeo * ext_ncomp)
!      integer h, i, j, k, m, n, p
!      call nuclear_potential(ngeo, ncor, ext_label, ext_ncomp, nucpot)
!      ! add requested component ranges to rspfunc
!      call permute_selcomp_add((1d0,0d0), nf, order, fields(:)%comp, &
!                               tcomp, fields(:)%ncomp, nucpot, rspfunc)
!
!! MaR: SIMPLE LOOPS TO ASSIGN VALUES WITH ONLY GEOMETRIC PERTURBATIONS
!
!      if (ngeo == nf) then
!
!         if (ngeo == 1) then
!
!            rspfunc_output = rspfunc
!
!         else if (ngeo == 2) then
!
!            h = 0
!            do i = 1, ncor
!               do j = i, ncor
!                  h = h + 1
!                  rspfunc_output(h) = rspfunc((i - 1)*ncor + j)
!               end do
!            end do
!
!         else if (ngeo == 3) then
!
!            h = 0
!            do i = 1, ncor
!               do j = i, ncor
!                  do k = j, ncor
!                     h = h + 1
!                     rspfunc_output(h) = rspfunc((i - 1)*ncor**2 + (j - 1)*ncor + k)
!                  end do
!               end do
!            end do
!
!         else if (ngeo == 4) then
!
!            h = 0
!            do i = 1, ncor
!               do j = i, ncor
!                  do k = j, ncor
!                     do m = k, ncor
!                        h = h + 1
!                        rspfunc_output(h) = rspfunc((i - 1)*ncor**3 + (j - 1)*ncor**2 + &
!                                            (k - 1)*ncor + m)
!                     end do
!                  end do
!               end do
!            end do
!
!         else if (ngeo == 5) then
!
!            h = 0
!            do i = 1, ncor
!               do j = i, ncor
!                  do k = j, ncor
!                     do m = k, ncor
!                        do n = m, ncor
!                           h = h + 1
!                           rspfunc_output(h) = rspfunc((i - 1)*ncor**4 + (j - 1)*ncor**3 + &
!                                                       (k - 1)*ncor**2 + (m - 1)*ncor + n)
!                        end do
!                     end do
!                  end do
!               end do
!            end do
!
!         else if (ngeo == 6) then
!
!            h = 0
!            do i = 1, ncor
!               do j = i, ncor
!                  do k = j, ncor
!                     do m = k, ncor
!                        do n = m, ncor
!                           do p = n, ncor
!                              h = h + 1
!                              rspfunc_output(h) = rspfunc((i - 1)*ncor**5 + &
!                              (j - 1)*ncor**4 + (k - 1)*ncor**3 + &
!                              (m - 1)*ncor**2 + (n - 1)*ncor + p)
!                           end do
!                        end do
!                     end do
!                  end do
!               end do
!            end do
!
!
!         else if (ngeo > 6) then
!
!            write(*,*) 'rsp_nucpot error: No support for ngeo > 6 yet'
!            call quit('rsp_nucpot error: No support for ngeo > 6 yet')
!
!         else
!
!            write(*,*) 'rsp_nucpot error: Unknown field setup'
!            call quit('rsp_nucpot error: Unknown field setup')
!
!         end if
!
!         else if (ngeo == (nf - 1)) then
!
!            if (nf == 2) then
!
!               write(*,*) 'rsp_nucpot warning: Generally untested support for one non-geo. field with nf = 2'
!
!               h = 0
!               do i = 1, ncor
!                  do j = 1, 3
!                     h = h + 1
!                     rspfunc_output(h) = rspfunc(ncor * (j - 1) + i)
!                  end do
!               end do
!
!            else if (nf == 1) then
!
!               write(*,*) 'rsp_nucpot warning: Generally untested support for one non-geo. field with nf = 1'
!
!               rspfunc_output = rspfunc
!                              
!            else
!
!               rspfunc_output = 0.0
!
!            end if
!
!         else
!
!            rspfunc_output = 0.0
!
!         end if
!
!       end subroutine
!  end subroutine
!
!
!  !> average f-perturbed overlap integrals with perturbed density D
!  !> and energy-weighted density DFD
!  subroutine rsp_ovlave(nf, f, c, nc, nblks, blk_info, & 
!                                      blk_sizes, propsize, DFD, ave)
!    !> number of fields
!    integer,       intent(in)  :: nf, propsize
!    integer :: nblks
!    integer, dimension(nblks) :: blk_sizes
!    integer, dimension(nblks, 3) :: blk_info
!    !> field labels in std order
!    character(4),  intent(in)  :: f(nf)
!    !> first and number of- components in each field
!    integer,       intent(in)  :: c(nf), nc(nf)
!    !> energy-weighted density matrix
!    type(matrix),  intent(in) :: DFD
!    !> output average
!    complex(8),    intent(inout) :: ave(propsize)
!
!    !----------------------------------------------
!
!    call interface_1el_ovlave_tr(nf, f, c, nc, nblks, blk_info, & 
!                                 blk_sizes, propsize, ave = ave, DFD = DFD)
!
!  end subroutine



  recursive subroutine rsp_ovlave_t_matrix_2014(num_fields, fields, bra, &
                                           ket, D, get_ovl_exp, propsize, ave)

    implicit none

    integer :: num_fields, propsize, i, j, k, under_construction, merged_nblks
    integer :: ave_offset, tmp_result_offset, merged_triang_size, tmp_ave_size
    integer :: total_num_perturbations, np_bra, np_ket
    integer, dimension(0) :: noc
    type(QcMat) :: D
    type(p_tuple) :: fields, bra, ket, merged_p_tuple, bra_static, ket_static
    external :: get_ovl_exp
    integer, dimension(bra%n_perturbations + ket%n_perturbations) :: pids_current_contribution
    complex(8), dimension(propsize) :: ave
    complex(8), allocatable, dimension(:) :: tmp_ave
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size, &
                                          blk_sizes_merged, translated_index, pert_ext_bra, &
                                          pert_ext_ket
    integer, allocatable, dimension(:,:) :: blk_sizes, merged_indices
    integer, allocatable, dimension(:,:,:) :: blks_tuple_info, merged_blk_info


under_construction = 0

if (under_construction == 1) then

! write(*,*) 'Called for T matrix contribution - currently unfinished so nothing was added'

else
    

    if (num_fields > 0) then

       call rsp_ovlave_t_matrix_2014(num_fields - 1, p_tuple_remove_first(fields), &
            merge_p_tuple(bra, p_tuple_getone(fields, 1)), ket, D, get_ovl_exp, propsize, ave)

       call rsp_ovlave_t_matrix_2014(num_fields - 1, p_tuple_remove_first(fields), &
            bra, merge_p_tuple(ket, p_tuple_getone(fields, 1)), D, get_ovl_exp, propsize, ave)

    else
   

      call p_tuple_to_external_tuple(bra, np_bra, pert_ext_bra)
      call p_tuple_to_external_tuple(ket, np_ket, pert_ext_ket)


       merged_p_tuple = p_tuple_standardorder(merge_p_tuple(bra, ket))

       ! Make frequency independent blocks for bra/ket
       call p_tuple_p1_cloneto_p2(bra, bra_static)
       call p_tuple_p1_cloneto_p2(ket, ket_static)

       bra_static%freq = 0.0
       ket_static%freq = 0.0

       allocate(nfields(2))
       allocate(nblks_tuple(2))
    
       nfields(1) = bra_static%n_perturbations
       nblks_tuple(1) = get_num_blks(bra_static)

       nfields(2) = ket_static%n_perturbations
       nblks_tuple(2) = get_num_blks(ket_static)
               
       total_num_perturbations = sum(nfields)

       allocate(blks_tuple_info(2, total_num_perturbations, 3))
       allocate(blks_tuple_triang_size(2))
       allocate(blk_sizes(2, total_num_perturbations))
       allocate(blk_sizes_merged(total_num_perturbations))
    
       blks_tuple_info(1, :, :) = get_blk_info(nblks_tuple(1), bra_static)
       blks_tuple_triang_size(1) = get_triangulated_size(nblks_tuple(1), &
                                   blks_tuple_info(1, 1:nblks_tuple(1), :))
       blk_sizes(1, 1:nblks_tuple(1)) = get_triangular_sizes(nblks_tuple(1), &
       blks_tuple_info(1,1:nblks_tuple(1),2), blks_tuple_info(1,1:nblks_tuple(1),3))

       blks_tuple_info(2, :, :) = get_blk_info(nblks_tuple(2), ket_static)
       blks_tuple_triang_size(2) = get_triangulated_size(nblks_tuple(2), &
                                   blks_tuple_info(2, 1:nblks_tuple(2), :))
       blk_sizes(2, 1:nblks_tuple(2)) = get_triangular_sizes(nblks_tuple(2), &
       blks_tuple_info(2,1:nblks_tuple(2),2), blks_tuple_info(2,1:nblks_tuple(2),3))


       ! Also make frequency dependent blocks for bra/ket
       ! The latter blocks will be larger or the same size as the former
       ! Loop over indices of the latter block, apply frequency factors and access
       ! elements of the former block (related to the interface/integral routines)
      
       ! MaR: Not completely sure if this is the correct size
       tmp_ave_size = product(blks_tuple_triang_size)

       allocate(tmp_ave(tmp_ave_size))
       tmp_ave = 0.0

       call get_ovl_exp(np_bra, pert_ext_bra, np_ket, pert_ext_ket, &
                        0, noc, 1, (/D/), tmp_ave_size, tmp_ave)
         
       
       k = 1
       do j = 1, bra_static%n_perturbations
          pids_current_contribution(k) = bra_static%pid(j)
          k = k + 1
       end do

       do j = 1, ket_static%n_perturbations
          pids_current_contribution(k) = ket_static%pid(j)
          k = k + 1
       end do

       merged_nblks = get_num_blks(merged_p_tuple)

       allocate(merged_blk_info(1, merged_nblks, 3))


       merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)
       blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
       merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
       merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)

       allocate(merged_indices(merged_triang_size, total_num_perturbations))
       allocate(translated_index(total_num_perturbations))

       call make_triangulated_indices(merged_nblks, merged_blk_info, & 
            merged_triang_size, merged_indices)

       do i = 1, size(merged_indices, 1)

          ave_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
                   (/sum(nfields)/), &
                   (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
                   (/merged_indices(i, :) /))
   
          do j = 1, total_num_perturbations
    
             translated_index(j) = merged_indices(i,pids_current_contribution(j))
    
          end do

          ! MaR: Will there ever be a completely unperturbed case? If so, it should be zero anyway
          ! because of no frequencies.
    
          if (bra_static%n_perturbations == 0) then
    
             tmp_result_offset = get_triang_blks_tuple_offset(1, &
                                 total_num_perturbations, nblks_tuple(2), &
                                 nfields(2), blks_tuple_info(2, :, :), &
                                 blk_sizes(2,:), blks_tuple_triang_size(2), & 
                                 (/ translated_index(:) /))
    
          elseif (ket_static%n_perturbations == 0) then
    
             tmp_result_offset = get_triang_blks_tuple_offset(1, &
                                 total_num_perturbations, nblks_tuple(1), &
                                 nfields(1), blks_tuple_info(1, :, :), &
                                 blk_sizes(1,:), blks_tuple_triang_size(1), & 
                                 (/ translated_index(:) /))

          else

             tmp_result_offset = get_triang_blks_tuple_offset(2, &
                                 total_num_perturbations, nblks_tuple, &
                                 nfields, blks_tuple_info, blk_sizes, &
                                 blks_tuple_triang_size, &
                                 (/ translated_index(:) /))
    
          end if

          ave(ave_offset) = ave(ave_offset) + &
          0.5 * (sum(bra%freq) - sum(ket%freq)) * tmp_ave(tmp_result_offset)

       end do

       deallocate(pert_ext_bra)
       deallocate(pert_ext_ket)
       
       deallocate(tmp_ave)
       deallocate(merged_indices)
       deallocate(translated_index)
       deallocate(nfields)
       deallocate(nblks_tuple)
       deallocate(blks_tuple_info)
       deallocate(blks_tuple_triang_size)
       deallocate(blk_sizes)
       deallocate(blk_sizes_merged)
       deallocate(merged_blk_info)

    end if

  end if

  end subroutine


!   ! MaR: This routine not tested - awaiting development in integral code
!   !> Compute half-differentiated overlap contribution to Fock matrices
  recursive subroutine rsp_ovlint_t_matrix_2014(num_fields, fields, bra, &
                                           ket, get_ovl_mat, propsize, fock)

    implicit none

    integer :: nr_ao, num_fields, propsize, i, j, k, under_construction, merged_nblks
    integer :: fock_offset, int_result_offset, merged_triang_size, tmp_fock_size
    integer :: total_num_perturbations, np_bra, np_ket
    integer, dimension(0) :: noc
    type(p_tuple) :: fields, bra, ket, merged_p_tuple, tester, bra_static, ket_static
    integer, dimension(bra%n_perturbations + ket%n_perturbations) :: pids_current_contribution
    type(QcMat), dimension(propsize) :: fock
    type(QcMat), allocatable, dimension(:) :: tmp_fock
    external :: get_ovl_mat
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size, &
                                          blk_sizes_merged, translated_index, pert_ext_bra, &
                                          pert_ext_ket
                                          
    integer, allocatable, dimension(:,:) :: blk_sizes, merged_indices
    integer, allocatable, dimension(:,:,:) :: blks_tuple_info, merged_blk_info


under_construction = 0

if (under_construction == 1) then

write(*,*) 'Called for T matrix contribution - currently unfinished so nothing was added'

else
    

    if (num_fields > 0) then

!        tester = p_tuple_getone(fields, 1)

       call rsp_ovlint_t_matrix_2014(num_fields - 1, p_tuple_remove_first(fields), &
            merge_p_tuple(bra, p_tuple_getone(fields, 1)), ket, get_ovl_mat, propsize, fock)

!        tester = p_tuple_remove_first(fields)
!        tester = p_tuple_getone(fields, 1)
!        tester = merge_p_tuple(ket, p_tuple_getone(fields, 1))

       call rsp_ovlint_t_matrix_2014(num_fields - 1, p_tuple_remove_first(fields), &
            bra, merge_p_tuple(ket, p_tuple_getone(fields, 1)), get_ovl_mat, propsize, fock)

! tester = p_tuple_getone(fields, 1)

    else

      call p_tuple_to_external_tuple(bra, np_bra, pert_ext_bra)
      call p_tuple_to_external_tuple(ket, np_ket, pert_ext_ket)


       merged_p_tuple = p_tuple_standardorder(merge_p_tuple(bra, ket))

       ! Make frequency independent blocks for bra/ket
       call p_tuple_p1_cloneto_p2(bra, bra_static)
       call p_tuple_p1_cloneto_p2(ket, ket_static)

       bra_static%freq = 0.0
       ket_static%freq = 0.0

       allocate(nfields(2))
       allocate(nblks_tuple(2))
    
       nfields(1) = bra_static%n_perturbations
       nblks_tuple(1) = get_num_blks(bra_static)

       nfields(2) = ket_static%n_perturbations
       nblks_tuple(2) = get_num_blks(ket_static)
               
       total_num_perturbations = sum(nfields)

       allocate(blks_tuple_info(2, total_num_perturbations, 3))
       allocate(blks_tuple_triang_size(2))
       allocate(blk_sizes(2, total_num_perturbations))
       allocate(blk_sizes_merged(total_num_perturbations))
    
       blks_tuple_info(1, :, :) = get_blk_info(nblks_tuple(1), bra_static)
       blks_tuple_triang_size(1) = get_triangulated_size(nblks_tuple(1), &
                                   blks_tuple_info(1, 1:nblks_tuple(1), :))
       blk_sizes(1, 1:nblks_tuple(1)) = get_triangular_sizes(nblks_tuple(1), &
       blks_tuple_info(1,1:nblks_tuple(1),2), blks_tuple_info(1,1:nblks_tuple(1),3))

       blks_tuple_info(2, :, :) = get_blk_info(nblks_tuple(2), ket_static)
       blks_tuple_triang_size(2) = get_triangulated_size(nblks_tuple(2), &
                                   blks_tuple_info(2, 1:nblks_tuple(2), :))
       blk_sizes(2, 1:nblks_tuple(2)) = get_triangular_sizes(nblks_tuple(2), &
       blks_tuple_info(2,1:nblks_tuple(2),2), blks_tuple_info(2,1:nblks_tuple(2),3))


       ! Also make frequency dependent blocks for bra/ket
       ! The latter blocks will be larger or the same size as the former
       ! Loop over indices of the latter block, apply frequency factors and access
       ! elements of the former block (related to the interface/integral routines)
      
       ! MaR: Not completely sure if this is the correct size
       tmp_fock_size = product(blks_tuple_triang_size)

       allocate(tmp_fock(tmp_fock_size))

       do i = 1, tmp_fock_size
          
          call QcMatInit(tmp_fock(i))

       end do

       call get_ovl_mat(np_bra, pert_ext_bra, np_ket, pert_ext_ket, &
                        0, noc, tmp_fock_size, tmp_fock)
       
       k = 1
       do j = 1, bra_static%n_perturbations
          pids_current_contribution(k) = bra_static%pid(j)
          k = k + 1
       end do

       do j = 1, ket_static%n_perturbations
          pids_current_contribution(k) = ket_static%pid(j)
          k = k + 1
       end do

       merged_nblks = get_num_blks(merged_p_tuple)

       allocate(merged_blk_info(1, merged_nblks, 3))


       merged_blk_info(1, :, :) = get_blk_info(merged_nblks, merged_p_tuple)
       blk_sizes_merged(1:merged_nblks) = get_triangular_sizes(merged_nblks, &
       merged_blk_info(1,1:merged_nblks,2), merged_blk_info(1,1:merged_nblks,3))
       merged_triang_size = get_triangulated_size(merged_nblks, merged_blk_info)

       allocate(merged_indices(merged_triang_size, total_num_perturbations))
       allocate(translated_index(total_num_perturbations))

       call make_triangulated_indices(merged_nblks, merged_blk_info, & 
            merged_triang_size, merged_indices)

       do i = 1, size(merged_indices, 1)

          fock_offset = get_triang_blks_tuple_offset(1, merged_nblks, (/merged_nblks/), &
                   (/sum(nfields)/), &
                   (/merged_blk_info/), blk_sizes_merged, (/merged_triang_size/), &
                   (/merged_indices(i, :) /))
   
          do j = 1, total_num_perturbations
    
             translated_index(j) = merged_indices(i,pids_current_contribution(j))
    
          end do

          ! MaR: Will there ever be a completely unperturbed case? If so, it should be zero anyway
          ! because of no frequencies.
    
          if (bra_static%n_perturbations == 0) then
    
             int_result_offset = get_triang_blks_tuple_offset(1, &
                                 total_num_perturbations, nblks_tuple(2), &
                                 nfields(2), blks_tuple_info(2, :, :), &
                                 blk_sizes(2,:), blks_tuple_triang_size(2), & 
                                 (/ translated_index(:) /))
    
          elseif (ket_static%n_perturbations == 0) then
    
             int_result_offset = get_triang_blks_tuple_offset(1, &
                                 total_num_perturbations, nblks_tuple(1), &
                                 nfields(1), blks_tuple_info(1, :, :), &
                                 blk_sizes(1,:), blks_tuple_triang_size(1), & 
                                 (/ translated_index(:) /))

          else

             int_result_offset = get_triang_blks_tuple_offset(2, &
                                 total_num_perturbations, nblks_tuple, &
                                 nfields, blks_tuple_info, blk_sizes, &
                                 blks_tuple_triang_size, &
                                 (/ translated_index(:) /))
    
          end if


          ! MaR: Frequency factor should be complex in general
          call QcMatrAXPY(dreal(0.5 * (sum(bra%freq) - sum(ket%freq))), tmp_fock(int_result_offset), fock(fock_offset))

!           fock(fock_offset) = fock(fock_offset) + &
!           0.5 * (sum(bra%freq) - sum(ket%freq)) * tmp_fock(int_result_offset)

       end do

       deallocate(pert_ext_bra)
       deallocate(pert_ext_ket)
       
       deallocate(merged_indices)
       deallocate(translated_index)
       deallocate(nfields)
       deallocate(nblks_tuple)
       deallocate(blks_tuple_info)
       deallocate(blks_tuple_triang_size)
       deallocate(blk_sizes)
       deallocate(blk_sizes_merged)
       deallocate(merged_blk_info)


       do i = 1, tmp_fock_size
          
          call QcMatDst(tmp_fock(i))

       end do

       deallocate(tmp_fock)

    end if

end if

  end subroutine


!
!
!  subroutine rsp_xcint_adapt(nr_ao, nf, f, c, nc, D, propsize, xcint)
!
!    integer :: i, j, k, nr_ao, propsize
!    !> number of fields
!    integer,       intent(in)    :: nf
!    !> field labels in std order
!    character(4),  intent(in)    :: f(nf)
!    !> first and number of- components in each field
!    integer,       intent(in)    :: c(nf), nc(nf)
!    !> output perturbed integrals
!    type(matrix),  intent(inout) :: xcint(propsize)
!type(matrix) :: Db, Fb
!    type(matrix), allocatable, dimension(:,:) :: tmp_xcint
!    type(matrix), dimension(:) :: D
!    !--------------------------------------------------
!
!    if (nf == 0) then
!      call rsp_xcint_interface(f, D, F=xcint(1))
!    else if (nf == 1) then
!       if (all(f==(/'GEO '/))) then
!          call rsp_xcint_interface(f, D, Fg=xcint)
!       else
!       end if
!    else if (nf == 2) then
!       if (all(f==(/'GEO ','GEO '/))) then
!          allocate(tmp_xcint(nc(1), nc(1)))
!          do i = 1, nc(1)
!             do j = 1, nc(1)
!!                 call mat_init(tmp_xcint(i,j), nr_ao, nr_ao)
!!        call mat_init_like_and_zero(xcint(1), tmp_xcint(i,j))
!
!             end do
!          end do
!          call rsp_xcint_interface(f, D, Fgg=tmp_xcint)
!          do i = 1, nc(1)
!             do j = 1, i
!                k = get_triang_blks_offset(1, 2, (/1, 2, nc(1)/), &
!                                           (/propsize/), (/i, j/))
!
!                xcint(k) = xcint(k) + tmp_xcint(i,j)
!                tmp_xcint(i, j) = 0
!             end do
!          end do
!          deallocate(tmp_xcint)
!       else
!       end if
!    else
!    end if
!
!  end subroutine
!
!
!
!  function index_of_field(f) result(i)
!    character(4) :: f
!    integer      :: i
!    do i = 1, size(all_known_fields)
!        if (all_known_fields(i)%label == f) return
!    end do
!    call quit('Field not found: ' // f)
!  end function
!
!
!  function rsp_field_anti(f)
!    character(4), intent(in) :: f(:)
!    logical :: rsp_field_anti(size(f))
!    integer :: i
!    rsp_field_anti = (/(all_known_fields(index_of_field(f(i)))%anti, i=1,size(f))/)
!  end function
!
!
!!  !> shape (dimensions) of property for fields f(:)
!!  function rsp_field_dim(f)
!!    character(4),  intent(in) :: f(:)
!!    integer :: rsp_field_dim(size(f)), i
!!    rsp_field_dim = (/(all_known_fields(index_of_field(f(i)))%ncomp, i=1,size(f))/)
!!    ! loop through mol-dependent
!!    do i = 1, size(f)
!!       if (rsp_field_dim(i) /= -1) then
!!          ! cycle
!!       else if (f(i) == 'GEO ') then
!!          rsp_field_dim(i) = 3 * get_num_atoms()
!!       else
!!          call quit('rsp_field_dim error: Number of comp. unknown for ' // f(i))
!!       end if
!!    end do
!!  end function
!
!
!  !> which fields are basis-perturbing
!  function rsp_field_bas(f)
!    character(4), intent(in) :: f(:)
!    logical :: rsp_field_bas(size(f))
!    integer :: i
!    rsp_field_bas = (/(all_known_fields(index_of_field(f(i)))%bas, i=1,size(f))/)
!  end function
!
!
!  function prefix_zeros(n, l)
!    integer, intent(in) :: n, l !number, length
!    character(l)        :: prefix_zeros !resulting n in ascii
!    character(1), parameter :: char0to9(0:9) &
!          = (/'0','1','2','3','4','5','6','7','8','9'/)
!    integer :: i, k
!    k = n
!    do i = l, 1, -1
!       prefix_zeros(i:i) = char0to9(mod(k,10))
!       k = k / 10
!    end do
!    if (k /= 0) call quit('prefix_zeros error: Argument integer does not fit ' &
!                       // 'in the specified number of ASCII caracters')
!  end function
!
!
!
!  ! determine 'canonical' ordering 'o' of fields 'f'. Also, validate fields%comp/ncomp
!  subroutine count_and_prepare(f, geo_first, bas_first, o, tot_ncomp, &
!                               nucpot, overlap, hamilt)
!    type(rsp_field),   intent(in)  :: f(:)
!    logical,           intent(in)  :: geo_first
!    logical,           intent(in)  :: bas_first
!    integer,           intent(out) :: o(size(f))
!    integer,           intent(out) :: tot_ncomp(size(f))
!    logical, optional, intent(out) :: nucpot, overlap, hamilt
!    integer i, j, k, jj, kk, nel, nmag, nbas, ngeo, fld_idx(size(f))
!    ! look up field indices
!    fld_idx(:) = (/(index_of_field(f(i)%label), i=1,size(f))/)
!    ! initialize order to unity
!    o(:) = (/(i,i=1,size(f))/)
!    ! selection sort, select field to be in i'th place from i:size(f)
!    do i = 1, size(f)
!       jj = i
!       j  = o(i)
!       do kk = i+1, size(f)
!          k = o(kk)
!          ! firstly, if geo_first, select any 'GEO ' ahead of others
!          if (geo_first .and. f(k)%label /= 'GEO ' &
!                        .and. f(j)%label == 'GEO ') cycle
!          if (.not.geo_first .or. (f(k)%label == 'GEO ' &
!                             .eqv. f(j)%label == 'GEO ')) then
!             ! secondly, if bas_first, those with field_stats%bas
!             if (bas_first .and. .not.all_known_fields(fld_idx(k))%bas &
!                           .and.      all_known_fields(fld_idx(j))%bas) cycle
!             if (.not.bas_first .or. (all_known_fields(fld_idx(k))%bas &
!                                .eqv. all_known_fields(fld_idx(j))%bas)) then
!                ! thirdly, pick highest fld_idx
!                if (fld_idx(k) <  fld_idx(j)) cycle
!                if (fld_idx(k) == fld_idx(j)) then
!                   ! fourthly, choose field with highest ncomp
!                   if (f(k)%ncomp <  f(j)%ncomp) cycle
!                   if (f(k)%ncomp == f(j)%ncomp) then
!                      ! fifthly, lowest comp
!                      if (f(k)%comp >  f(j)%comp) cycle
!                      if (f(k)%comp == f(j)%comp) then
!                         ! sixthly, highest abs(re(freq))
!                         if (abs(real(f(k)%freq)) <  abs(real(f(j)%freq))) cycle
!                         if (abs(real(f(k)%freq)) == abs(real(f(j)%freq))) then
!                            ! seventhly, lowest sign(re(freq))
!                            if (real(f(k)%freq) >  real(f(j)%freq)) cycle
!                            if (real(f(k)%freq) == real(f(j)%freq)) then
!                               ! eigthly, highest abs(im(freq))
!                               if (abs(aimag(f(k)%freq)) <  abs(aimag(f(j)%freq))) cycle
!                               if (abs(aimag(f(k)%freq)) == abs(aimag(f(j)%freq))) then
!                                  ! ninethly, lowest sign(im(freq))
!                                  if (aimag(f(k)%freq) >  aimag(f(j)%freq)) cycle
!                                  if (aimag(f(k)%freq) == aimag(f(j)%freq)) then
!                                     ! tenthly, if fields j and k are *identical*, go for the
!                                     ! one with lowest input position
!                                     if (k > j) cycle
!                                  end if
!                               end if
!                            end if
!                         end if
!                      end if
!                   end if
!                end if
!             end if
!          end if
!          jj = kk
!          j  = k  !new minimum
!       end do
!       !swap entries i and j
!       k     = o(i)
!       o(i)  = j
!       o(jj) = k
!    end do
!    ! place the full number of components in tot_ncomp
!    do i = 1, size(f)
!       tot_ncomp(i) = all_known_fields(fld_idx(o(i)))%ncomp
!       ! ajt FIXME -1 means ncomp is uninitialized (molecule-dependent)
!       !           This init should be done in a separate routine, together
!       !           with all the other fields' ncomp
!       if (tot_ncomp(i) == -1 .and. f(o(i))%label == 'GEO ') then
!          tot_ncomp(i) = 3*get_num_atoms()
!          all_known_fields(fld_idx(o(i)))%ncomp = tot_ncomp(i)
!       end if
!    end do
!
!    ! radovan: initialize
!    nmag = 0
!    nel  = 0
!    nbas = 0
!    ngeo = 0
!
!    ! if requested, determine whether the nuclear potential, the overlap
!    ! integrals, or the Hamiltonian, are preturbed by this tuple of fields
!    if (present(nucpot)) &
!       ngeo = count(f(:)%label == 'GEO ')
!    if (present(overlap) .or. present(hamilt)) &
!       nbas = count((/(all_known_fields(fld_idx(i))%bas, i=1,size(f))/))
!    if (present(nucpot) .or. present(hamilt)) then
!       nel  = count((/(all_known_fields(fld_idx(i))%lin,  i=1,size(f))/))
!       nmag = count((/(all_known_fields(fld_idx(i))%quad, i=1,size(f))/))
!    end if
!    if (present(nucpot)) &
!       nucpot = (((nel <= 1 .and. nmag == 0) &
!             .or. (nel == 0 .and. nmag <= 2)) &
!            .and. ngeo + nel + nmag == size(f))
!    if (present(overlap)) &
!       overlap = (nbas == size(f))
!    if (present(hamilt)) &
!       nucpot = (((nel <= 1 .and. nmag == 0) &
!             .or. (nel == 0 .and. nmag <= 2)) &
!            .and. nbas + nel + nmag == size(f))
!  end subroutine
!
!
!
!  function rsp_field_ordering(f) result(o)
!    type(rsp_field), intent(in) :: f(:)
!    integer                     :: o(size(f))
!    integer tc(size(f))
!    call count_and_prepare(f, .false., .false., o, tc)
!  end function
!
!
!  !> Add the average contribution 'ave', with dimensions 'dima' to response tensor
!  !> 'rsp', with dimensions 'dimr', but where dimensions are permuted by 'perm'.
!  !> 'idxr' selects (if >0) which starting components in 'ave' for each dim in 'rsp',
!  !> or (if <0) index of a density matrix dimension.
!  !> @prefac  sign or prefactor: rsp = rsp + prefac * ave
!  !> @ndim    number of dimensions
!  !> @perm    ordering of dimensions. For each dim in ave, the corresponding dim in rsp
!  !> @idxr    starting component in ave for each dim in rsp, or if negative, the index
!  !>          of a density dimension in rsp
!  !> @dima    dimensions of ave. Density dimensions must have dim 1
!  !> @dimr    dimensions of rsp
!  !> @ave     real array of integral averages, as from integral program
!  !> @rsp     complex respons function tensor
!  subroutine permute_selcomp_add(prefac, ndim, perm, idxr, dima, dimr, ave, rsp)
!    complex(8), intent(in)    :: prefac
!    integer,    intent(in)    :: ndim, perm(ndim), idxr(ndim)
!    integer,    intent(in)    :: dima(ndim), dimr(ndim)
!    real(8),    intent(in)    :: ave(product(dima))
!    complex(8), intent(inout) :: rsp(product(dimr))
!    integer i, ia, ir, stpr(ndim), stpa(ndim), ii(ndim), dd(ndim)
!    ! calculate dimension steps in ave (cumulative products of
!    ! dimensions), as well as offset due to starting indices idxr (ia)
!    ia = 0
!    do i = 1, ndim
!       if (i==1) stpa(i) = 1
!       if (i/=1) stpa(i) = stpa(i-1) * dima(i-1)
!       if (idxr(perm(i)) > 0) & !positive means starting comp
!          ia = ia + stpa(i) * (idxr(perm(i)) - 1)
!    end do
!    ! calculate (permuted) dimension steps in rsp, and offset due to
!    ! density indices (ir), and permuted dimensions (dd)
!    ir = 0
!    do i = 1, ndim
!       if (i==1) stpr(perm(i)) = 1
!       if (i/=1) stpr(perm(i)) = stpr(perm(i-1)) * dimr(i-1)
!       if (idxr(i) <= 0) then !negative means density index
!          ir = ir + stpr(perm(i)) * (-idxr(i) - 1)
!          dd(perm(i)) = 1
!       else
!          dd(perm(i)) = dimr(i)
!       end if
!    end do
!    ! loop over indices in ave and rsp
!    ii = 0 !indices from zero to dd-1
!    do
!       rsp(ir+1) = rsp(ir+1) + prefac * ave(ia+1)
!       ! increment indices
!       do i = 1, ndim
!          ii(i) = ii(i) + 1
!          ia = ia + stpa(i)
!          ir = ir + stpr(i)
!          if (ii(i) /= dd(i)) exit
!          ii(i) = 0
!          ia = ia - stpa(i) * dd(i)
!          ir = ir - stpr(i) * dd(i)
!       end do
!       if (i == ndim+1) exit
!    end do
!  end subroutine
!
!  subroutine rsp_pe(nr_ao, nf, f, c, nc, dens, propsize, fock)
!
!    use pe_variables, only: peqm, pe_gspol
!    use interface_pelib
!
!    integer                     :: i, j, k, nr_ao, propsize
!    !> number of fields
!    integer, intent(in)         :: nf
!    !> field labels in std order
!    character(4), intent(in)    :: f(nf)
!    !> first and number of- components in each field
!    integer, intent(in)         :: c(nf), nc(nf)
!    !> output perturbed integrals
!    type(matrix), intent(inout) :: fock(propsize)
!    !> perturbed density matrix
!    type(matrix), intent(in)    :: dens
!
!    real(8), allocatable        :: pe_dmat(:,:)
!    real(8), allocatable        :: pe_fmat(:,:)
!
!    if (.not. peqm .or. pe_gspol) return
!
!    if (any(f /= 'EL  ')) then
!        stop 'ERROR: PE-OpenRSP not implemented for other than EL.'
!    end if
!
!    if (any(f == 'EL  ')) then
!       return
!    end if
!
!    if (nf == 0) then
!        if (propsize /= 1) stop 'ERROR: propsize /= 1'
!        allocate(pe_dmat(nr_ao,nr_ao))
!        allocate(pe_fmat(nr_ao,nr_ao))
!        pe_fmat = 0.0d0
!        pe_dmat = 0.0d0
!        call daxpy(nr_ao*nr_ao, 1.0d0, dens%elms, 1, pe_dmat, 1)
!        call pe_response_operator(pe_dmat, pe_fmat, nr_ao, propsize)
!        deallocate(pe_dmat)
!        do i = 1, propsize
!            call daxpy(nr_ao*nr_ao, 2.0d0, pe_fmat, 1, fock(i)%elms, 1)
!        end do
!        deallocate(pe_fmat)
!    end if
!
!  end subroutine rsp_pe

end module
