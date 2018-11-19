! Copyright 2012      Gao Bin
!           2012      Radovan Bast
!           2009-2012 Andreas J. Thorvaldsen
! This source code form is subject to the terms of the
! GNU Lesser General Public License, version 2.1.
! If a copy of the GNU LGPL v2.1 was not distributed with this
! code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

!> @file Contains module rsp_contribs

!> This module contains routines for calculating contributions
!> to molecular properties (1st order, linear response, etc.),
!> and perturbed Fock matrices.

! APR 2016: Not yet adapted for new host framework

module rsp_contribs

  use rsp_field_tuple, only: p_tuple, p_tuple_remove_first, p_tuple_getone, &
                             p_tuple_standardorder, merge_p_tuple, &
                             p1_cloneto_p2, p_tuple_to_external_tuple
  use rsp_indices_and_addressing
  use qcmatrix_f
                                        
  implicit none

  public rsp_ovlave_t_matrix_2014
  public rsp_ovlint_t_matrix_2014

  private

  contains

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
    integer, dimension(bra%npert + ket%npert) :: pids_current_contribution
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
       call p1_cloneto_p2(bra, bra_static)
       call p1_cloneto_p2(ket, ket_static)

       bra_static%freq = 0.0
       ket_static%freq = 0.0

       allocate(nfields(2))
       allocate(nblks_tuple(2))
    
       nfields(1) = bra_static%npert
       nblks_tuple(1) = get_num_blks(bra_static)

       nfields(2) = ket_static%npert
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
       do j = 1, bra_static%npert
          pids_current_contribution(k) = bra_static%pid(j)
          k = k + 1
       end do

       do j = 1, ket_static%npert
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
    
          if (bra_static%npert == 0) then
    
             tmp_result_offset = get_triang_blks_tuple_offset(1, &
                                 total_num_perturbations, nblks_tuple(2), &
                                 nfields(2), blks_tuple_info(2, :, :), &
                                 blk_sizes(2,:), blks_tuple_triang_size(2), & 
                                 (/ translated_index(:) /))
    
          elseif (ket_static%npert == 0) then
    
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
    integer, dimension(bra%npert + ket%npert) :: pids_current_contribution
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
       call p1_cloneto_p2(bra, bra_static)
       call p1_cloneto_p2(ket, ket_static)

       bra_static%freq = 0.0
       ket_static%freq = 0.0

       allocate(nfields(2))
       allocate(nblks_tuple(2))
    
       nfields(1) = bra_static%npert
       nblks_tuple(1) = get_num_blks(bra_static)

       nfields(2) = ket_static%npert
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
          
          call QcMatInit(tmp_fock(i), fock(1))

       end do
       call get_ovl_mat(np_bra, pert_ext_bra, np_ket, pert_ext_ket, &
                        0, noc, tmp_fock_size, tmp_fock)
       k = 1
       do j = 1, bra_static%npert
          pids_current_contribution(k) = bra_static%pid(j)
          k = k + 1
       end do

       do j = 1, ket_static%npert
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
    
          if (bra_static%npert == 0) then
    
             int_result_offset = get_triang_blks_tuple_offset(1, &
                                 total_num_perturbations, nblks_tuple(2), &
                                 nfields(2), blks_tuple_info(2, :, :), &
                                 blk_sizes(2,:), blks_tuple_triang_size(2), & 
                                 (/ translated_index(:) /))
    
          elseif (ket_static%npert == 0) then
    
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

end module
