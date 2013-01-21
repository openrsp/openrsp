! commented out by Bin Gao, 2012-11-05
!#ifdef PRG_DIRAC
#define GRCONT_NOT_AVAILABLE
!#endif

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

  use matrix_defop
  use matrix_lowlevel, only: mat_init
  use interface_molecule
  use interface_io
  use interface_xc
  use interface_f77_memory
  use interface_1el
  use interface_scf
  use interface_basis
  use dalton_ifc
  use nuc_contributions
  use basis_set,  only: cgto
  use rsp_field_tuple
  use rsp_sdf_caching

  ! MaR: QUICK-FIX USE STATEMENT TO GET SUPPORT FOR DUMMY rsp_cfg TYPE
  use rsp_perturbed_matrices



  implicit none

  public rsp_nucpot
  public rsp_ovlave
  public rsp_oneave
  public rsp_twoave
  public rsp_ovlint
  public rsp_oneint
  public rsp_twoint
  public rsp_field
  public rsp_field_bas
  public rsp_field_dim
  public rsp_field_anti
  public rsp_field_ordering

  ! MaR: ROUTINES RETURNING NON-REDUNDANT VALUES W.R.T TENSOR SYMMETRY
  ! MaR: THESE COULD REPLACE THE CORRESPONDING ABOVE ROUTINES WHEN APPROPRIATE
  public rsp_nucpot_tr
  public rsp_ovlave_tr
  public rsp_ovlave_t_matrix
  public rsp_oneave_tr
  public rsp_twoave_tr
  public rsp_ovlint_tr
  public rsp_ovlint_t_matrix
  public rsp_oneint_tr
  public rsp_twoint_tr
  public rsp_xcint_tr_adapt

  !> Type describing a single field in a response function
  !> or response equation. A response equation (or density)
  !> corresponds to an array of prop_field. Similarly
  !> a response function corresponds to an array of prop_field
  !> whose freqs sum to zero.
  type rsp_field
     sequence
     !> 4-char pert label
     character(4) :: label
     !> frequency
     complex(8)   :: freq
     !> first component
     integer      :: comp
     !> number of components
     integer      :: ncomp
  end type

  !> private struct to collect properties of perturbing "fields"
  type field_stats
     !> four-letter abbreviation
     character(4)  :: label
     !> long name
     character(64) :: name
     !> number of components (when known, -1 otherwise)
     integer       :: ncomp
     !> anti-symmetric (1,3,5th ord.) perturbed integrals
     logical       :: anti
     !> basis dependent (sa. GEO and MAG)
     logical       :: bas
     !> one-electron operator linear in field strength (EL)
     logical       :: lin
     !> one-electron operator quadratic in field strength (MAGO)
     logical       :: quad
  end type


  ! to compactify the table below
  logical, parameter :: T = .true.
  logical, parameter :: F = .false.


  !> ajt nov09: AUX0..AUX9 are 10 configurable basis-independent 1-electron
  !>            perturbations, configured by setting the corresponding
  !>            HERMIT integral label in prop_auxlab(0:9).
  !> ajt jan10: EXCI is a ZERO (no) perturbation, and is introduced to
  !>            allow the same code to contract response functions and
  !>            "generalized transition moments".
  !> ajt may10: FREQ is also a ZERO (no) perturbation, and is introduced to
  !>            allow the same code to contract response functions and
  !>            frequency-differentiated response functions.
  type(field_stats) :: all_known_fields(12) = &                  !nc an ba ln qu
     (/field_stats('EXCI', 'Generalized "excitation" field'      , 1, F, F, T, T), &
       field_stats('FREQ', 'Generalized "freqency" field'        , 1, F, F, T, T), &
       field_stats('EL  ', 'Electric field'                      , 3, F, F, T, F), &
       field_stats('VEL ', 'Velocity'                            , 3, T, F, T, F), &
       field_stats('MAGO', 'Magnetic field w/o. London orbitals' , 3, T, F, F, T), &
       field_stats('MAG ', 'Magnetic field with London orbitals' , 3, T, T, F, F), &
       field_stats('ELGR', 'Electric field gradient'             , 6, F, F, T, F), &
       field_stats('VIBM', 'Displacement along vibrational modes',-1, F, T, F, F), &
       field_stats('GEO ', 'Nuclear coordinates'                 ,-1, F, T, F, F), & !-1=mol-dep
       field_stats('NUCM', 'Nuclear magnetic moment'             ,-1, F, T, F, T), & !-1=mol-dep
       field_stats('AOCC', 'AO contraction coefficients'         ,-1, F, T, F, F), & !-1=mol-dep
       field_stats('AOEX', 'AO exponents'                        ,-1, F, T, F, F)/)  !-1=mol-dep

  character(1), parameter :: xyz(3) = (/'X','Y','Z'/)

  private

contains


  !> Contribution from nuclear repulsion and nuclei--field interaction
  !> to response functions. Fields (type rsp_field) are here in arbitrary order.
  !> (in normal mode) Fields are sorted, component ranges extended to what
  !> rsp_backend expects (currently only full ranges). Then the call is then relayed
  !> to rsp_backend's nuclear_potential, which computes and returns the requested real
  !> tensor in standard order. The requested ranges of this tensor is then reordered
  !> and added to rspfunc
  subroutine rsp_nucpot(fields, rspfunc)
    !> field descriptors (label freq comp ncomp)
    type(rsp_field), intent(in)    :: fields(:)
    !> output tensor, to which nuclear contribution is *ADDED*
    complex(8),      intent(inout) :: rspfunc(product(fields%ncomp))
    !---------------------------------------------------------------
    integer      nf, ncor, ngeo, ext_ncomp, i
    integer      order(size(fields)), tcomp(size(fields))
    character(4) ext_label(2)
    logical      nonz
    ! prepare and determine ordering. This also validates comp/ncomp
    call count_and_prepare(fields, .true., .false., order, tcomp, nucpot = nonz)
    ! early return if zero
    if (.not.nonz) return
    ! find number of GEO, and two, one or none external fields
    nf = size(fields)
    ngeo = count(fields%label == 'GEO ')
    ext_label(:) = (/'NONE','NONE'/)
    ext_ncomp = 1
    if (ngeo == nf-1) then
       ext_label(1) = fields(order(nf))%label
       ext_ncomp    = tcomp(nf)
    else if (ngeo == nf-2) then
       ext_label(1) = fields(order(nf-1))%label
       ext_label(2) = fields(order(nf))%label
       ext_ncomp    = tcomp(nf-1) * tcomp(nf)
    end if
    ! use inner to avoid allocate'ing 'nucpot' below
    ncor = 1
    if (ngeo > 0) ncor = tcomp(1)
    call inner
  contains
    subroutine inner
      real(8) nucpot(ncor**ngeo * ext_ncomp)
      call nuclear_potential(ngeo, ncor, ext_label, ext_ncomp, nucpot)
      ! add requested component ranges to rspfunc
      call permute_selcomp_add((1d0,0d0), nf, order, fields(:)%comp, &
                               tcomp, fields(:)%ncomp, nucpot, rspfunc)
    end subroutine
  end subroutine


  ! MaR: Seems to work properly, but memory usage is not tensor symmetry nonredundant
  !> Contribution from nuclear repulsion and nuclei--field interaction
  !> to response functions. Fields (type rsp_field) are here in arbitrary order.
  !> (in normal mode) Fields are sorted, component ranges extended to what
  !> rsp_backend expects (currently only full ranges). Then the call is then relayed
  !> to rsp_backend's nuclear_potential, which computes and returns the requested real
  !> tensor in standard order. The requested ranges of this tensor is then reordered
  !> and added to rspfunc
  subroutine rsp_nucpot_tr(fields, propsize, rspfunc_output)
    !> field descriptors (label freq comp ncomp)
    type(rsp_field), intent(in)    :: fields(:)
    !> output tensor, to which nuclear contribution is *ADDED*
    integer                        :: propsize
    complex(8),      intent(inout) :: rspfunc_output(propsize)
    !> tmp tensor, to which nuclear contribution is *ADDED*
    complex(8) :: rspfunc(product(fields%ncomp))
    !---------------------------------------------------------------
    integer      nf, ncor, ngeo, ext_ncomp, i
    integer      order(size(fields)), tcomp(size(fields))
    character(4) ext_label(2)
    logical      nonz
    rspfunc = 0.0
    ! prepare and determine ordering. This also validates comp/ncomp
    call count_and_prepare(fields, .true., .false., order, tcomp, nucpot = nonz)
    ! early return if zero
    if (.not.nonz) return
    ! find number of GEO, and two, one or none external fields
    nf = size(fields)
    ngeo = count(fields%label == 'GEO ')
    ext_label(:) = (/'NONE','NONE'/)
    ext_ncomp = 1
    if (ngeo == nf-1) then
       ext_label(1) = fields(order(nf))%label
       ext_ncomp    = tcomp(nf)
    else if (ngeo == nf-2) then
       ext_label(1) = fields(order(nf-1))%label
       ext_label(2) = fields(order(nf))%label
       ext_ncomp    = tcomp(nf-1) * tcomp(nf)
    end if
    ! use inner to avoid allocate'ing 'nucpot' below
    ncor = 1
    if (ngeo > 0) ncor = tcomp(1)
    call inner
  contains
    subroutine inner
      real(8) nucpot(ncor**ngeo * ext_ncomp)
      integer h, i, j, k, m, n, p
      call nuclear_potential(ngeo, ncor, ext_label, ext_ncomp, nucpot)
      ! add requested component ranges to rspfunc
      call permute_selcomp_add((1d0,0d0), nf, order, fields(:)%comp, &
                               tcomp, fields(:)%ncomp, nucpot, rspfunc)

! MaR: SIMPLE LOOPS TO ASSIGN VALUES WITH ONLY GEOMETRIC PERTURBATIONS

      if (ngeo == nf) then

         if (ngeo == 1) then

            rspfunc_output = rspfunc

         else if (ngeo == 2) then

            h = 0
            do i = 1, ncor
               do j = i, ncor
                  h = h + 1
                  rspfunc_output(h) = rspfunc((i - 1)*ncor + j)
               end do
            end do

         else if (ngeo == 3) then

            h = 0
            do i = 1, ncor
               do j = i, ncor
                  do k = j, ncor
                     h = h + 1
                     rspfunc_output(h) = rspfunc((i - 1)*ncor**2 + (j - 1)*ncor + k)
                  end do
               end do
            end do

         else if (ngeo == 4) then

            h = 0
            do i = 1, ncor
               do j = i, ncor
                  do k = j, ncor
                     do m = k, ncor
                        h = h + 1
                        rspfunc_output(h) = rspfunc((i - 1)*ncor**3 + (j - 1)*ncor**2 + &
                                            (k - 1)*ncor + m)
                     end do
                  end do
               end do
            end do

         else if (ngeo == 5) then

            h = 0
            do i = 1, ncor
               do j = i, ncor
                  do k = j, ncor
                     do m = k, ncor
                        do n = m, ncor
                           h = h + 1
                           rspfunc_output(h) = rspfunc((i - 1)*ncor**4 + (j - 1)*ncor**3 + &
                                                       (k - 1)*ncor**2 + (m - 1)*ncor + n)
                        end do
                     end do
                  end do
               end do
            end do

         else if (ngeo == 6) then

            h = 0
            do i = 1, ncor
               do j = i, ncor
                  do k = j, ncor
                     do m = k, ncor
                        do n = m, ncor
                           do p = n, ncor
                              h = h + 1
                              rspfunc_output(h) = rspfunc((i - 1)*ncor**5 + &
                              (j - 1)*ncor**4 + (k - 1)*ncor**3 + &
                              (m - 1)*ncor**2 + (n - 1)*ncor + p)
                           end do
                        end do
                     end do
                  end do
               end do
            end do


         else if (ngeo > 6) then

            write(*,*) 'rsp_nucpot_tr error: No support for ngeo > 6 yet'
            call quit('rsp_nucpot_tr error: No support for ngeo > 6 yet')

         else

            write(*,*) 'rsp_nucpot_tr error: Unknown field setup'
            call quit('rsp_nucpot_tr error: Unknown field setup')

         end if

         else if (ngeo == (nf - 1)) then

            if (nf == 2) then

               write(*,*) 'rsp_nucpot_tr warning: Generally untested support for one non-geo. field with nf = 2'

               h = 0
               do i = 1, ncor
                  do j = 1, 3
                     h = h + 1
                     rspfunc_output(h) = rspfunc(ncor * (j - 1) + i)
                  end do
               end do

            else if (nf == 1) then

               write(*,*) 'rsp_nucpot_tr warning: Generally untested support for one non-geo. field with nf = 1'

               rspfunc_output = rspfunc
                              
            else

               rspfunc_output = 0.0

            end if

         else

            rspfunc_output = 0.0

         end if

       end subroutine
  end subroutine


  !> average f-perturbed overlap integrals with perturbed density D
  !> and energy-weighted density DFD
  subroutine rsp_ovlave(nf, f, c, nc, DFD, ave, w, D)
    !> number of fields
    integer,       intent(in)  :: nf
    !> field labels in std order
    character(4),  intent(in)  :: f(nf)
    !> first and number of- components in each field
    integer,       intent(in)  :: c(nf), nc(nf)
    !> energy-weighted density matrix
    type(matrix),  intent(in)  :: DFD
    !> output average
    complex(8),    intent(out) :: ave(product(nc))
    !> field frequencies corresponding to each field
    complex(8),    intent(in), optional  :: w(nf)
    !> density matrix to contract half-differentiated overlap against
    type(matrix),  intent(in), optional  :: D
    call interface_1el_ovlave(nf, f, c, nc, DFD, ave, w, D)
  end subroutine


  !> average f-perturbed overlap integrals with perturbed density D
  !> and energy-weighted density DFD
  subroutine rsp_ovlave_tr(nf, f, c, nc, nblks, blk_info, & 
                                      blk_sizes, propsize, DFD, ave)
    !> number of fields
    integer,       intent(in)  :: nf, propsize
    integer :: nblks
    integer, dimension(nblks) :: blk_sizes
    integer, dimension(nblks, 3) :: blk_info
    !> field labels in std order
    character(4),  intent(in)  :: f(nf)
    !> first and number of- components in each field
    integer,       intent(in)  :: c(nf), nc(nf)
    !> energy-weighted density matrix
    type(matrix),  intent(in) :: DFD
    !> output average
    complex(8),    intent(inout) :: ave(propsize)

    !----------------------------------------------

    call interface_1el_ovlave_tr(nf, f, c, nc, nblks, blk_info, & 
                                 blk_sizes, propsize, ave = ave, DFD = DFD)

  end subroutine









!   ! MaR: This routine not tested - awaiting development in integral code
!   !> Compute half-differentiated overlap contribution to property
  recursive subroutine rsp_ovlave_t_matrix(num_fields, fields, bra, &
                                           ket, D, propsize, ave)

    implicit none

    integer :: num_fields, propsize, i, j, k, under_construction, merged_nblks
    integer :: ave_offset, tmp_result_offset, merged_triang_size, tmp_ave_size
    integer :: total_num_perturbations
    type(matrix) :: D
    type(p_tuple) :: fields, bra, ket, merged_p_tuple
    integer, dimension(bra%n_perturbations + ket%n_perturbations) :: pids_current_contribution
    complex(8), dimension(propsize) :: ave, tmp_ave
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size, &
                                          blk_sizes_merged, translated_index
    integer, allocatable, dimension(:,:) :: blk_sizes, merged_indices
    integer, allocatable, dimension(:,:,:) :: blks_tuple_info, merged_blk_info


under_construction = 1

if (under_construction == 1) then

! write(*,*) 'Called for T matrix contribution - currently unfinished so nothing was added'

else
    

    if (num_fields > 0) then

       call rsp_ovlave_t_matrix(num_fields - 1, p_tuple_remove_first(fields), &
            merge_p_tuple(bra, p_tuple_getone(fields, 1)), ket, D, propsize, ave)

       call rsp_ovlave_t_matrix(num_fields - 1, p_tuple_remove_first(fields), &
            bra, merge_p_tuple(ket, p_tuple_getone(fields, 1)), D, propsize, ave)

    else

       tmp_ave = 0.0

       merged_p_tuple = merge_p_tuple(bra, ket)

       ! Make frequency independent blocks for bra/ket

       bra%freq = 0.0
       ket%freq = 0.0

       allocate(nfields(2))
       allocate(nblks_tuple(2))
    
       nfields(1) = bra%n_perturbations
       nblks_tuple(1) = get_num_blks(bra)

       nfields(2) = ket%n_perturbations
       nblks_tuple(2) = get_num_blks(ket)
               
       total_num_perturbations = sum(nfields)

       allocate(blks_tuple_info(2, total_num_perturbations, 3))
       allocate(blks_tuple_triang_size(2))
       allocate(blk_sizes(2, total_num_perturbations))
       allocate(blk_sizes_merged(total_num_perturbations))
    
       blks_tuple_info(1, :, :) = get_blk_info(nblks_tuple(1), bra)
       blks_tuple_triang_size(1) = get_triangulated_size(nblks_tuple(1), &
                                   blks_tuple_info(1, 1:nblks_tuple(1), :))
       blk_sizes(1, 1:nblks_tuple(1)) = get_triangular_sizes(nblks_tuple(1), &
       blks_tuple_info(1,1:nblks_tuple(1),2), blks_tuple_info(1,1:nblks_tuple(1),3))

       blks_tuple_info(2, :, :) = get_blk_info(nblks_tuple(2), ket)
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

       call interface_1el_ovlave_half_diff(bra%n_perturbations, bra%plab, &
            (/(j/j, j = 1, bra%n_perturbations)/), bra%pdim, ket%n_perturbations, ket%plab, &
            (/(j/j, j = 1, ket%n_perturbations)/), ket%pdim, nblks_tuple, &
            blks_tuple_info, blk_sizes, D, tmp_ave_size, tmp_ave)

       k = 1
       do j = 1, bra%n_perturbations
          pids_current_contribution(k) = bra%pid(j)
          k = k + 1
       end do

       do j = 1, ket%n_perturbations
          pids_current_contribution(k) = ket%pid(j)
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
    
          if (bra%n_perturbations == 0) then
    
             tmp_result_offset = get_triang_blks_tuple_offset(1, &
                                 total_num_perturbations, nblks_tuple(2), &
                                 nfields(2), blks_tuple_info(2, :, :), &
                                 blk_sizes(2,:), blks_tuple_triang_size(2), & 
                                 (/ translated_index(:) /))
    
          elseif (ket%n_perturbations == 0) then
    
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











  !> Average 1-electron integrals perturbed by fields f
  !> with the (perturbed) density matrix D

! radovan: there is code repetition in ave and int setup

  subroutine rsp_oneave(nf, f, c, nc, D, ave)
    !> number of fields
    integer,       intent(in)  :: nf
    !> field labels in std order
    character(4),  intent(in)  :: f(nf)
    !> first and number of- components in each field
    integer,       intent(in)  :: c(nf), nc(nf)
    !> density matrix to average over
    type(matrix),  intent(in)  :: D
    !> output average
    complex(8),    intent(out) :: ave(product(nc))
    call interface_1el_oneave(nf, f, c, nc, D, ave)
  end subroutine



  !> Average 1-electron integrals perturbed by fields f
  !> with the (perturbed) density matrix D

! radovan: there is code repetition in ave and int setup

  subroutine rsp_oneave_tr(nf, f, c, nc, D, nblks, blk_info, blk_sizes, propsize, ave)
!     use dalton_ifc, only: SHELLS_NUCLEI_displace
    ! Gen1Int interface in Dalton
!     use gen1int_api
    !> structure containing integral program settings
!     type(rsp_cfg), intent(in)  :: mol
    !> number of fields
    integer,       intent(in)  :: nf
    integer :: nblks
    integer, dimension(nblks) :: blk_sizes
    integer, dimension(nblks, 3) :: blk_info
    !> field labels in std order
    character(4),  intent(in)  :: f(nf)
    !> first and number of- components in each field
    integer,       intent(in)  :: c(nf), nc(nf), propsize
    !> density matrix to average over
    type(matrix),  intent(in)  :: D
    !> output average
    complex(8),    intent(out) :: ave(propsize)
    call interface_1el_oneave_tr(nf, f, c, nc, D, nblks, blk_info, & 
                                 blk_sizes, propsize, ave)

  end subroutine

  function rank_one_pointer(siz, arr) result(ptr)
     integer,         intent(in) :: siz
     real(8), target, intent(in) :: arr(siz)
     real(8), pointer            :: ptr(:)
     ptr => arr
  end function


  !> Average 2-electron integrals perturbed by fields f over the
  !> product of density matrces D1 and D2
  subroutine rsp_twoave(nf, f, c, nc, D1, D2, ave)

    use eri_contractions, only: ctr_arg
    use eri_basis_loops,  only: unopt_geodiff_loop
    use interface_interest

    !> number of fields
    integer,              intent(in)  :: nf
    !> field labels in std order
    character(4),         intent(in)  :: f(nf)
    !> first and number of- components in each field
    integer,              intent(in)  :: c(nf), nc(nf)
    !> density matrix to average over
    type(matrix), target, intent(in)  :: D1, D2
    !> output average
    complex(8),           intent(out) :: ave(product(nc))
    !----------------------------------------------
    real(8), allocatable              :: real_ave(:)
    real(8), pointer :: tmp(:,:,:,:) !scratch
    type(matrix)  A(1) !scratch matrices
    type(ctr_arg) arg(1)
    real(8)       r
    integer       i, j, k, l, n, ncor

  if (any(f== 'EL  ')) then
     ave = 0.0
  else

    if (nf==0) then
       ! contract second density to Fock matrix, then trace with first
       A(1) = 0*D1
       call mat_ensure_alloc(A(1), only_alloc=.true.)
       call interface_scf_get_g(D2, A(1)) !Coulomb and exchange
       ave(1) = trace(A(1),D1)
    else if (nf==1 .and. f(1)=='GEO ') then

#ifdef PRG_DALTON
       ncor = 3 * get_nr_atoms()
       allocate(tmp(ncor,1,1,1))
#ifdef GRCONT_NOT_AVAILABLE
       arg(1) = ctr_arg(1, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor, tmp(:,1,1,1)))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       tmp = 2.0d0*tmp
#else
       n = D1%nrow
       f77_memory(     :n*n)   = reshape(D1%elms,(/n*n/))
       f77_memory(n*n+1:n*n*2) = reshape(D2%elms,(/n*n/))
       call GRCONT(f77_memory(n*n*2+1:), size(f77_memory)-n*n*2, &
                   tmp(:,1,1,1), ncor, .true., .false., &
                   1, 0, .true., .false., f77_memory(:n*n*2), 2)
#endif
       ave(:nc(1)) = tmp(c(1):c(1)+nc(1)-1,1,1,1)
       deallocate(tmp)
#endif /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
       allocate(real_ave(size(ave)))
       real_ave = 0.0
       call interest_mpi_wake_up()
       call interest_get_int(D1%nrow, D1%elms, D2%elms, 1, 0, size(real_ave), real_ave)
       do i = 1, size(ave)
          ave(i) = 2.0d0*real_ave(i)
       end do
       deallocate(real_ave)
#endif /* ifdef PRG_DIRAC */

    else if (nf==2 .and. all(f==(/'GEO ','GEO '/))) then

#ifdef PRG_DALTON
       ncor = 3 * get_nr_atoms()
       allocate(tmp(ncor,ncor,1,1))
#ifdef GRCONT_NOT_AVAILABLE
       arg(1) = ctr_arg(2, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor**2, tmp(:,:,1,1)))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       ! symmetrize
       do j = 1, ncor
          do i = 1, j
             r = tmp(i, j, 1, 1) + tmp(j, i, 1, 1)
             tmp(i, j, 1, 1) = 2.0d0*r
             tmp(j, i, 1, 1) = 2.0d0*r
          end do
       end do
#else
       n = D1%nrow
       f77_memory(     :n*n)   = reshape(D1%elms,(/n*n/))
       f77_memory(n*n+1:n*n*2) = reshape(D2%elms,(/n*n/))
       call GRCONT(f77_memory(n*n*2+1:), size(f77_memory)-n*n*2, &
                   tmp(:,:,1,1), ncor**2, .true., .false., &
                   2, 0, .true., .false., f77_memory(:n*n*2), 2)
#endif
       ave = reshape(tmp(c(1):c(1)+nc(1)-1, &
                         c(2):c(2)+nc(2)-1,1,1), shape(ave))
       deallocate(tmp)
#endif /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
       allocate(real_ave(size(ave)))
       real_ave = 0.0
       call interest_mpi_wake_up()
       call interest_get_int(D1%nrow, D1%elms, D2%elms, 2, 0, size(real_ave), real_ave)
       do i = 1, size(ave)
          ave(i) = 2.0d0*real_ave(i)
       end do
       deallocate(real_ave)
#endif /* ifdef PRG_DIRAC */

    else if (nf==3 .and. all(f==(/'GEO ','GEO ','GEO '/))) then

#ifdef PRG_DALTON
       ! contract FULL cubic in tmp, unsymmetrized divided by six
       ncor = 3 * get_nr_atoms()
       allocate(tmp(ncor,ncor,ncor,1))
       arg(1) = ctr_arg(3, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor**3, tmp(:,:,:,1)))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       ! symmetrize
       do k = 1, ncor
          do j = 1, k
             do i = 1, j
                r = tmp(i,j,k,1) + tmp(i,k,j,1) + tmp(k,i,j,1) &
                  + tmp(k,j,i,1) + tmp(j,k,i,1) + tmp(j,i,k,1)
                tmp(i,j,k,1) = r;  tmp(i,k,j,1) = r;  tmp(k,i,j,1) = r
                tmp(k,j,i,1) = r;  tmp(j,k,i,1) = r;  tmp(j,i,k,1) = r
             end do
          end do
       end do
       ! extract requested block
       ave = 2 * reshape(tmp(c(1):c(1)+nc(1)-1, &
                             c(2):c(2)+nc(2)-1, &
                             c(3):c(3)+nc(3)-1, 1), shape(ave))
       deallocate(tmp)
#endif /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
       allocate(real_ave(size(ave)))
       real_ave = 0.0
       call interest_mpi_wake_up()
       call interest_get_int(D1%nrow, D1%elms, D2%elms, 3, 0, size(real_ave), real_ave)
       do i = 1, size(ave)
          ave(i) = 2.0d0*real_ave(i)
       end do
       deallocate(real_ave)
#endif /* ifdef PRG_DIRAC */

    else if (nf==4 .and. all(f==(/'GEO ','GEO ','GEO ','GEO '/))) then
       ncor = 3 * get_nr_atoms()
       allocate(tmp(ncor,ncor,ncor,ncor))
       ! contract FULL quartic in tmp, unsymmetrized divided by 24
       arg(1) = ctr_arg(4, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor**4, tmp))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       ! symmetrize
       do l = 1, ncor
          do k = 1, l
             do j = 1, k
                do i = 1, j
                   r = tmp(i,j,k,l) + tmp(i,k,j,l) + tmp(k,i,j,l) &
                     + tmp(k,j,i,l) + tmp(j,k,i,l) + tmp(j,i,k,l) &
                     + tmp(i,j,l,k) + tmp(i,k,l,j) + tmp(k,i,l,j) &
                     + tmp(k,j,l,i) + tmp(j,k,l,i) + tmp(j,i,l,k) &
                     + tmp(i,l,j,k) + tmp(i,l,k,j) + tmp(k,l,i,j) &
                     + tmp(k,l,j,i) + tmp(j,l,k,i) + tmp(j,l,i,k) &
                     + tmp(l,i,j,k) + tmp(l,i,k,j) + tmp(l,k,i,j) &
                     + tmp(l,k,j,i) + tmp(l,j,k,i) + tmp(l,j,i,k)
                   tmp(i,j,k,l) = r;  tmp(i,k,j,l) = r;  tmp(k,i,j,l) = r
                   tmp(k,j,i,l) = r;  tmp(j,k,i,l) = r;  tmp(j,i,k,l) = r
                   tmp(i,j,l,k) = r;  tmp(i,k,l,j) = r;  tmp(k,i,l,j) = r
                   tmp(k,j,l,i) = r;  tmp(j,k,l,i) = r;  tmp(j,i,l,k) = r
                   tmp(i,l,j,k) = r;  tmp(i,l,k,j) = r;  tmp(k,l,i,j) = r
                   tmp(k,l,j,i) = r;  tmp(j,l,k,i) = r;  tmp(j,l,i,k) = r
                   tmp(l,i,j,k) = r;  tmp(l,i,k,j) = r;  tmp(l,k,i,j) = r
                   tmp(l,k,j,i) = r;  tmp(l,j,k,i) = r;  tmp(l,j,i,k) = r
                end do
             end do
          end do
       end do
       ! extract requested block
       ave = 2 * reshape(tmp(c(1):c(1)+nc(1)-1, &
                             c(2):c(2)+nc(2)-1, &
                             c(3):c(3)+nc(3)-1, &
                             c(4):c(4)+nc(4)-1), shape(ave))
       deallocate(tmp)
    else
       print *, 'rsp_twoave error: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('rsp_twoave error: not implented or in wrong order')
    end if

  end if

  end subroutine


! MR: NOT SURE IF WORKING PROPERLY
  !> Average 2-electron integrals perturbed by fields f over the
  !> product of density matrces D1 and D2
  subroutine rsp_twoave_tr(nf, f, c, nc, D1, D2, propsize, ave)

    use eri_contractions, only: ctr_arg
    use eri_basis_loops,  only: unopt_geodiff_loop
    use interface_interest

    !> number of fields
    integer,              intent(in)  :: nf, propsize
    !> field labels in std order
    character(4),         intent(in)  :: f(nf)
    !> first and number of- components in each field
    integer,              intent(in)  :: c(nf), nc(nf)
    !> density matrix to average over
    type(matrix), target, intent(in)  :: D1, D2
    !> output average
    complex(8),           intent(out) :: ave(propsize)
    !----------------------------------------------
    real(8), allocatable              :: real_ave(:)
    real(8), pointer :: tmp(:,:,:,:) !scratch
    real(8), pointer :: tmp_5(:,:,:,:,:) !scratch
    real(8), pointer :: tmp_6(:,:,:,:,:,:) !scratch
    type(matrix)  A(1) !scratch matrices
    type(ctr_arg) arg(1)
    real(8)       r
    integer       h, i, j, k, l, m, n, p, q, ncor

  if (any(f== 'EL  ')) then
     ave = 0.0
  else

    if (nf==0) then
       ! contract second density to Fock matrix, then trace with first
       A(1) = 0*D1
       call mat_ensure_alloc(A(1), only_alloc=.true.)
       call interface_scf_get_g(D2, A(1)) !Coulomb and exchange
       ave(1) = trace(A(1),D1)
    else if (nf==1 .and. f(1)=='GEO ') then

#ifdef PRG_DALTON
       ncor = 3 * get_nr_atoms()
       allocate(tmp(ncor,1,1,1))
       tmp = 0.0
#ifdef GRCONT_NOT_AVAILABLE
       arg(1) = ctr_arg(1, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor, tmp(:,1,1,1)))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       tmp = 2.0d0*tmp
#else
       n = D1%nrow
       f77_memory(     :n*n)   = reshape(D1%elms,(/n*n/))
       f77_memory(n*n+1:n*n*2) = reshape(D2%elms,(/n*n/))
       call GRCONT(f77_memory(n*n*2+1:), size(f77_memory)-n*n*2, &
                   tmp(:,1,1,1), ncor, .true., .false., &
                   1, 0, .true., .false., f77_memory(:n*n*2), 2)
#endif
       ave(:nc(1)) = tmp(c(1):c(1)+nc(1)-1,1,1,1)
       deallocate(tmp)
#endif /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
       allocate(real_ave(size(ave)))
       real_ave = 0.0
       call interest_mpi_wake_up()
       call interest_get_int(D1%nrow, D1%elms, D2%elms, 1, 0, size(real_ave), real_ave)
       do i = 1, size(ave)
          ave(i) = 2.0d0*real_ave(i)
       end do
       deallocate(real_ave)
#endif /* ifdef PRG_DIRAC */

    else if (nf==2 .and. all(f==(/'GEO ','GEO '/))) then
       ncor = 3 * get_nr_atoms()
       allocate(tmp(ncor,ncor,1,1))
       tmp = 0.0
#ifdef GRCONT_NOT_AVAILABLE
       arg(1) = ctr_arg(2, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor**2, tmp(:,:,1,1)))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       ! symmetrize
       do j = 1, ncor
          do i = 1, j
             r = tmp(i, j, 1, 1) + tmp(j, i, 1, 1)
             tmp(i, j, 1, 1) = 2.0d0*r
             tmp(j, i, 1, 1) = 2.0d0*r
          end do
       end do
#else
       n = D1%nrow
       f77_memory(     :n*n)   = reshape(D1%elms,(/n*n/))
       f77_memory(n*n+1:n*n*2) = reshape(D2%elms,(/n*n/))
       call GRCONT(f77_memory(n*n*2+1:), size(f77_memory)-n*n*2, &
                   tmp(:,:,1,1), ncor**2, .true., .false., &
                   2, 0, .true., .false., f77_memory(:n*n*2), 2)
#endif

       h = 0
       do i = 1, ncor
          do j = i, ncor
             h = h + 1
             ave(h) = tmp(i,j,1,1)
          end do
       end do

       deallocate(tmp)

    else if (nf==3 .and. all(f==(/'GEO ','GEO ','GEO '/))) then
       ! contract FULL cubic in tmp, unsymmetrized divided by six
       ncor = 3 * get_nr_atoms()
       allocate(tmp(ncor,ncor,ncor,1))
       tmp = 0.0
       arg(1) = ctr_arg(3, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor**3, tmp(:,:,:,1)))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       ! symmetrize
       do k = 1, ncor
          do j = 1, k
             do i = 1, j
                r = tmp(i,j,k,1) + tmp(i,k,j,1) + tmp(k,i,j,1) &
                  + tmp(k,j,i,1) + tmp(j,k,i,1) + tmp(j,i,k,1)
                tmp(i,j,k,1) = r
             end do
          end do
       end do

       h = 0
       do i = 1, ncor
          do j = i, ncor
             do k = j, ncor
                h = h + 1
                ! MR: UNSURE ABOUT ORIGINS OF FACTOR 2 IN NEXT LINE
                ave(h) = 2 * tmp(i,j,k,1)
             end do
          end do
       end do

       deallocate(tmp)

    else if (nf==4 .and. all(f==(/'GEO ','GEO ','GEO ','GEO '/))) then
       ncor = 3 * get_nr_atoms()
       allocate(tmp(ncor,ncor,ncor,ncor))
       tmp = 0.0
       ! contract FULL quartic in tmp, unsymmetrized divided by 24
       arg(1) = ctr_arg(4, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor**4, tmp))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       ! symmetrize
       do l = 1, ncor
          do k = 1, l
             do j = 1, k
                do i = 1, j
                   r = tmp(i,j,k,l) + tmp(i,k,j,l) + tmp(k,i,j,l) &
                     + tmp(k,j,i,l) + tmp(j,k,i,l) + tmp(j,i,k,l) &
                     + tmp(i,j,l,k) + tmp(i,k,l,j) + tmp(k,i,l,j) &
                     + tmp(k,j,l,i) + tmp(j,k,l,i) + tmp(j,i,l,k) &
                     + tmp(i,l,j,k) + tmp(i,l,k,j) + tmp(k,l,i,j) &
                     + tmp(k,l,j,i) + tmp(j,l,k,i) + tmp(j,l,i,k) &
                     + tmp(l,i,j,k) + tmp(l,i,k,j) + tmp(l,k,i,j) &
                     + tmp(l,k,j,i) + tmp(l,j,k,i) + tmp(l,j,i,k)
                   tmp(i,j,k,l) = r
                end do
             end do
          end do
       end do

       h = 0
       do i = 1, ncor
          do j = i, ncor
             do k = j, ncor
                do m = k, ncor
                   h = h + 1
                   ! MR: UNSURE ABOUT ORIGINS OF FACTOR 2 IN NEXT LINE
                   ave(h) = 2 * tmp(i,j,k,m)
                end do
             end do
          end do
       end do

       deallocate(tmp)

    else if (nf==5 .and. all(f==(/'GEO ','GEO ','GEO ','GEO ', 'GEO '/))) then
       ncor = 3 * get_nr_atoms()
       allocate(tmp_5(ncor,ncor,ncor,ncor,ncor))
       tmp_5 = 0.0
       ! contract FULL quartic in tmp, unsymmetrized divided by 24
       arg(1) = ctr_arg(5, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor**5, tmp_5))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       ! symmetrize
       do m = 1, ncor
          do l = 1, m
             do k = 1, l
                do j = 1, k
                   do i = 1, j

                      r = tmp_5(i,j,k,l,m) + tmp_5(i,j,k,m,l) + tmp_5(i,j,l,k,m) + &
                          tmp_5(i,j,l,m,k) + tmp_5(i,j,m,k,l) + tmp_5(i,j,m,l,k) + &
                          tmp_5(i,k,j,l,m) + tmp_5(i,k,j,m,l) + tmp_5(i,k,l,j,m) + &
                          tmp_5(i,k,l,m,j) + tmp_5(i,k,m,j,l) + tmp_5(i,k,m,l,j) + &
                          tmp_5(i,l,j,k,m) + tmp_5(i,l,j,m,k) + tmp_5(i,l,k,j,m) + &
                          tmp_5(i,l,k,m,j) + tmp_5(i,l,m,j,k) + tmp_5(i,l,m,k,j) + &
                          tmp_5(i,m,j,k,l) + tmp_5(i,m,j,l,k) + tmp_5(i,m,k,j,l) + &
                          tmp_5(i,m,k,l,j) + tmp_5(i,m,l,j,k) + tmp_5(i,m,l,k,j) + &
                          tmp_5(j,i,k,l,m) + tmp_5(j,i,k,m,l) + tmp_5(j,i,l,k,m) + &
                          tmp_5(j,i,l,m,k) + tmp_5(j,i,m,k,l) + tmp_5(j,i,m,l,k) + &
                          tmp_5(j,k,i,l,m) + tmp_5(j,k,i,m,l) + tmp_5(j,k,l,i,m) + &
                          tmp_5(j,k,l,m,i) + tmp_5(j,k,m,i,l) + tmp_5(j,k,m,l,i) + &
                          tmp_5(j,l,i,k,m) + tmp_5(j,l,i,m,k) + tmp_5(j,l,k,i,m) + &
                          tmp_5(j,l,k,m,i) + tmp_5(j,l,m,i,k) + tmp_5(j,l,m,k,i) + &
                          tmp_5(j,m,i,k,l) + tmp_5(j,m,i,l,k) + tmp_5(j,m,k,i,l) + &
                          tmp_5(j,m,k,l,i) + tmp_5(j,m,l,i,k) + tmp_5(j,m,l,k,i) + &
                          tmp_5(k,i,j,l,m) + tmp_5(k,i,j,m,l) + tmp_5(k,i,l,j,m) + &
                          tmp_5(k,i,l,m,j) + tmp_5(k,i,m,j,l) + tmp_5(k,i,m,l,j) + &
                          tmp_5(k,j,i,l,m) + tmp_5(k,j,i,m,l) + tmp_5(k,j,l,i,m) + &
                          tmp_5(k,j,l,m,i) + tmp_5(k,j,m,i,l) + tmp_5(k,j,m,l,i) + &
                          tmp_5(k,l,i,j,m) + tmp_5(k,l,i,m,j) + tmp_5(k,l,j,i,m) + &
                          tmp_5(k,l,j,m,i) + tmp_5(k,l,m,i,j) + tmp_5(k,l,m,j,i) + &
                          tmp_5(k,m,i,j,l) + tmp_5(k,m,i,l,j) + tmp_5(k,m,j,i,l) + &
                          tmp_5(k,m,j,l,i) + tmp_5(k,m,l,i,j) + tmp_5(k,m,l,j,i) + &
                          tmp_5(l,i,j,k,m) + tmp_5(l,i,j,m,k) + tmp_5(l,i,k,j,m) + &
                          tmp_5(l,i,k,m,j) + tmp_5(l,i,m,j,k) + tmp_5(l,i,m,k,j) + &
                          tmp_5(l,j,i,k,m) + tmp_5(l,j,i,m,k) + tmp_5(l,j,k,i,m) + &
                          tmp_5(l,j,k,m,i) + tmp_5(l,j,m,i,k) + tmp_5(l,j,m,k,i) + &
                          tmp_5(l,k,i,j,m) + tmp_5(l,k,i,m,j) + tmp_5(l,k,j,i,m) + &
                          tmp_5(l,k,j,m,i) + tmp_5(l,k,m,i,j) + tmp_5(l,k,m,j,i) + &
                          tmp_5(l,m,i,j,k) + tmp_5(l,m,i,k,j) + tmp_5(l,m,j,i,k) + &
                          tmp_5(l,m,j,k,i) + tmp_5(l,m,k,i,j) + tmp_5(l,m,k,j,i) + &
                          tmp_5(m,i,j,k,l) + tmp_5(m,i,j,l,k) + tmp_5(m,i,k,j,l) + &
                          tmp_5(m,i,k,l,j) + tmp_5(m,i,l,j,k) + tmp_5(m,i,l,k,j) + &
                          tmp_5(m,j,i,k,l) + tmp_5(m,j,i,l,k) + tmp_5(m,j,k,i,l) + &
                          tmp_5(m,j,k,l,i) + tmp_5(m,j,l,i,k) + tmp_5(m,j,l,k,i) + &
                          tmp_5(m,k,i,j,l) + tmp_5(m,k,i,l,j) + tmp_5(m,k,j,i,l) + &
                          tmp_5(m,k,j,l,i) + tmp_5(m,k,l,i,j) + tmp_5(m,k,l,j,i) + &
                          tmp_5(m,l,i,j,k) + tmp_5(m,l,i,k,j) + tmp_5(m,l,j,i,k) + &
                          tmp_5(m,l,j,k,i) + tmp_5(m,l,k,i,j) + tmp_5(m,l,k,j,i) 
                          tmp_5(i,j,k,l,m) = r

                   end do
                end do
             end do
          end do
       end do

       h = 0
       do i = 1, ncor
          do j = i, ncor
             do k = j, ncor
                do l = k, ncor
                   do m = l, ncor
                      h = h + 1
                      ! MR: UNSURE ABOUT ORIGINS OF FACTOR 2 IN NEXT LINE
                      ave(h) = 2 * tmp_5(i,j,k,l,m)
                   end do
                end do
             end do
          end do
       end do

       deallocate(tmp_5)

    else if (nf==6 .and. all(f==(/'GEO ','GEO ','GEO ','GEO ', 'GEO ', 'GEO '/))) then
       ncor = 3 * get_nr_atoms()
       allocate(tmp_6(ncor,ncor,ncor,ncor,ncor,ncor))
       tmp_6 = 0.0
       ! contract FULL quartic in tmp, unsymmetrized divided by 24
       arg(1) = ctr_arg(6, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor**6, tmp_6))
       call unopt_geodiff_loop(basis_large, &
                               basis_small, &
                               arg)
       ! symmetrize
       do n = 1, ncor
          do m = 1, n
             do l = 1, m
                do k = 1, l
                   do j = 1, k
                      do i = 1, j

                         r = tmp_6(i,j,k,l,m,n) + tmp_6(i,j,k,l,n,m) + tmp_6(i,j,k,m,l,n) + &
                             tmp_6(i,j,k,m,n,l) + tmp_6(i,j,k,n,l,m) + tmp_6(i,j,k,n,m,l) + &
                             tmp_6(i,j,l,k,m,n) + tmp_6(i,j,l,k,n,m) + tmp_6(i,j,l,m,k,n) + &
                             tmp_6(i,j,l,m,n,k) + tmp_6(i,j,l,n,k,m) + tmp_6(i,j,l,n,m,k) + &
                             tmp_6(i,j,m,k,l,n) + tmp_6(i,j,m,k,n,l) + tmp_6(i,j,m,l,k,n) + &
                             tmp_6(i,j,m,l,n,k) + tmp_6(i,j,m,n,k,l) + tmp_6(i,j,m,n,l,k) + &
                             tmp_6(i,j,n,k,l,m) + tmp_6(i,j,n,k,m,l) + tmp_6(i,j,n,l,k,m) + &
                             tmp_6(i,j,n,l,m,k) + tmp_6(i,j,n,m,k,l) + tmp_6(i,j,n,m,l,k) + &
                             tmp_6(i,k,j,l,m,n) + tmp_6(i,k,j,l,n,m) + tmp_6(i,k,j,m,l,n) + &
                             tmp_6(i,k,j,m,n,l) + tmp_6(i,k,j,n,l,m) + tmp_6(i,k,j,n,m,l) + &
                             tmp_6(i,k,l,j,m,n) + tmp_6(i,k,l,j,n,m) + tmp_6(i,k,l,m,j,n) + &
                             tmp_6(i,k,l,m,n,j) + tmp_6(i,k,l,n,j,m) + tmp_6(i,k,l,n,m,j) + &
                             tmp_6(i,k,m,j,l,n) + tmp_6(i,k,m,j,n,l) + tmp_6(i,k,m,l,j,n) + &
                             tmp_6(i,k,m,l,n,j) + tmp_6(i,k,m,n,j,l) + tmp_6(i,k,m,n,l,j) + &
                             tmp_6(i,k,n,j,l,m) + tmp_6(i,k,n,j,m,l) + tmp_6(i,k,n,l,j,m) + &
                             tmp_6(i,k,n,l,m,j) + tmp_6(i,k,n,m,j,l) + tmp_6(i,k,n,m,l,j) + &
                             tmp_6(i,l,j,k,m,n) + tmp_6(i,l,j,k,n,m) + tmp_6(i,l,j,m,k,n) + &
                             tmp_6(i,l,j,m,n,k) + tmp_6(i,l,j,n,k,m) + tmp_6(i,l,j,n,m,k) + &
                             tmp_6(i,l,k,j,m,n) + tmp_6(i,l,k,j,n,m) + tmp_6(i,l,k,m,j,n) + &
                             tmp_6(i,l,k,m,n,j) + tmp_6(i,l,k,n,j,m) + tmp_6(i,l,k,n,m,j) + &
                             tmp_6(i,l,m,j,k,n) + tmp_6(i,l,m,j,n,k) + tmp_6(i,l,m,k,j,n) + &
                             tmp_6(i,l,m,k,n,j) + tmp_6(i,l,m,n,j,k) + tmp_6(i,l,m,n,k,j) + &
                             tmp_6(i,l,n,j,k,m) + tmp_6(i,l,n,j,m,k) + tmp_6(i,l,n,k,j,m) + &
                             tmp_6(i,l,n,k,m,j) + tmp_6(i,l,n,m,j,k) + tmp_6(i,l,n,m,k,j) + &
                             tmp_6(i,m,j,k,l,n) + tmp_6(i,m,j,k,n,l) + tmp_6(i,m,j,l,k,n) + &
                             tmp_6(i,m,j,l,n,k) + tmp_6(i,m,j,n,k,l) + tmp_6(i,m,j,n,l,k) + &
                             tmp_6(i,m,k,j,l,n) + tmp_6(i,m,k,j,n,l) + tmp_6(i,m,k,l,j,n) + &
                             tmp_6(i,m,k,l,n,j) + tmp_6(i,m,k,n,j,l) + tmp_6(i,m,k,n,l,j) + &
                             tmp_6(i,m,l,j,k,n) + tmp_6(i,m,l,j,n,k) + tmp_6(i,m,l,k,j,n) + &
                             tmp_6(i,m,l,k,n,j) + tmp_6(i,m,l,n,j,k) + tmp_6(i,m,l,n,k,j) + &
                             tmp_6(i,m,n,j,k,l) + tmp_6(i,m,n,j,l,k) + tmp_6(i,m,n,k,j,l) + &
                             tmp_6(i,m,n,k,l,j) + tmp_6(i,m,n,l,j,k) + tmp_6(i,m,n,l,k,j) + &
                             tmp_6(i,n,j,k,l,m) + tmp_6(i,n,j,k,m,l) + tmp_6(i,n,j,l,k,m) + &
                             tmp_6(i,n,j,l,m,k) + tmp_6(i,n,j,m,k,l) + tmp_6(i,n,j,m,l,k) + &
                             tmp_6(i,n,k,j,l,m) + tmp_6(i,n,k,j,m,l) + tmp_6(i,n,k,l,j,m) + &
                             tmp_6(i,n,k,l,m,j) + tmp_6(i,n,k,m,j,l) + tmp_6(i,n,k,m,l,j) + &
                             tmp_6(i,n,l,j,k,m) + tmp_6(i,n,l,j,m,k) + tmp_6(i,n,l,k,j,m) + &
                             tmp_6(i,n,l,k,m,j) + tmp_6(i,n,l,m,j,k) + tmp_6(i,n,l,m,k,j) + &
                             tmp_6(i,n,m,j,k,l) + tmp_6(i,n,m,j,l,k) + tmp_6(i,n,m,k,j,l) + &
                             tmp_6(i,n,m,k,l,j) + tmp_6(i,n,m,l,j,k) + tmp_6(i,n,m,l,k,j) + &
                             tmp_6(j,i,k,l,m,n) + tmp_6(j,i,k,l,n,m) + tmp_6(j,i,k,m,l,n) + &
                             tmp_6(j,i,k,m,n,l) + tmp_6(j,i,k,n,l,m) + tmp_6(j,i,k,n,m,l) + &
                             tmp_6(j,i,l,k,m,n) + tmp_6(j,i,l,k,n,m) + tmp_6(j,i,l,m,k,n) + &
                             tmp_6(j,i,l,m,n,k) + tmp_6(j,i,l,n,k,m) + tmp_6(j,i,l,n,m,k) + &
                             tmp_6(j,i,m,k,l,n) + tmp_6(j,i,m,k,n,l) + tmp_6(j,i,m,l,k,n) + &
                             tmp_6(j,i,m,l,n,k) + tmp_6(j,i,m,n,k,l) + tmp_6(j,i,m,n,l,k) + &
                             tmp_6(j,i,n,k,l,m) + tmp_6(j,i,n,k,m,l) + tmp_6(j,i,n,l,k,m) + &
                             tmp_6(j,i,n,l,m,k) + tmp_6(j,i,n,m,k,l) + tmp_6(j,i,n,m,l,k) + &
                             tmp_6(j,k,i,l,m,n) + tmp_6(j,k,i,l,n,m) + tmp_6(j,k,i,m,l,n) + &
                             tmp_6(j,k,i,m,n,l) + tmp_6(j,k,i,n,l,m) + tmp_6(j,k,i,n,m,l) + &
                             tmp_6(j,k,l,i,m,n) + tmp_6(j,k,l,i,n,m) + tmp_6(j,k,l,m,i,n) + &
                             tmp_6(j,k,l,m,n,i) + tmp_6(j,k,l,n,i,m) + tmp_6(j,k,l,n,m,i) + &
                             tmp_6(j,k,m,i,l,n) + tmp_6(j,k,m,i,n,l) + tmp_6(j,k,m,l,i,n) + &
                             tmp_6(j,k,m,l,n,i) + tmp_6(j,k,m,n,i,l) + tmp_6(j,k,m,n,l,i) + &
                             tmp_6(j,k,n,i,l,m) + tmp_6(j,k,n,i,m,l) + tmp_6(j,k,n,l,i,m) + &
                             tmp_6(j,k,n,l,m,i) + tmp_6(j,k,n,m,i,l) + tmp_6(j,k,n,m,l,i) + &
                             tmp_6(j,l,i,k,m,n) + tmp_6(j,l,i,k,n,m) + tmp_6(j,l,i,m,k,n) + &
                             tmp_6(j,l,i,m,n,k) + tmp_6(j,l,i,n,k,m) + tmp_6(j,l,i,n,m,k) + &
                             tmp_6(j,l,k,i,m,n) + tmp_6(j,l,k,i,n,m) + tmp_6(j,l,k,m,i,n) + &
                             tmp_6(j,l,k,m,n,i) + tmp_6(j,l,k,n,i,m) + tmp_6(j,l,k,n,m,i) + &
                             tmp_6(j,l,m,i,k,n) + tmp_6(j,l,m,i,n,k) + tmp_6(j,l,m,k,i,n) + &
                             tmp_6(j,l,m,k,n,i) + tmp_6(j,l,m,n,i,k) + tmp_6(j,l,m,n,k,i) + &
                             tmp_6(j,l,n,i,k,m) + tmp_6(j,l,n,i,m,k) + tmp_6(j,l,n,k,i,m) + &
                             tmp_6(j,l,n,k,m,i) + tmp_6(j,l,n,m,i,k) + tmp_6(j,l,n,m,k,i) + &
                             tmp_6(j,m,i,k,l,n) + tmp_6(j,m,i,k,n,l) + tmp_6(j,m,i,l,k,n) + &
                             tmp_6(j,m,i,l,n,k) + tmp_6(j,m,i,n,k,l) + tmp_6(j,m,i,n,l,k) + &
                             tmp_6(j,m,k,i,l,n) + tmp_6(j,m,k,i,n,l) + tmp_6(j,m,k,l,i,n) + &
                             tmp_6(j,m,k,l,n,i) + tmp_6(j,m,k,n,i,l) + tmp_6(j,m,k,n,l,i) + &
                             tmp_6(j,m,l,i,k,n) + tmp_6(j,m,l,i,n,k) + tmp_6(j,m,l,k,i,n) + &
                             tmp_6(j,m,l,k,n,i) + tmp_6(j,m,l,n,i,k) + tmp_6(j,m,l,n,k,i) + &
                             tmp_6(j,m,n,i,k,l) + tmp_6(j,m,n,i,l,k) + tmp_6(j,m,n,k,i,l) + &
                             tmp_6(j,m,n,k,l,i) + tmp_6(j,m,n,l,i,k) + tmp_6(j,m,n,l,k,i) + &
                             tmp_6(j,n,i,k,l,m) + tmp_6(j,n,i,k,m,l) + tmp_6(j,n,i,l,k,m) + &
                             tmp_6(j,n,i,l,m,k) + tmp_6(j,n,i,m,k,l) + tmp_6(j,n,i,m,l,k) + &
                             tmp_6(j,n,k,i,l,m) + tmp_6(j,n,k,i,m,l) + tmp_6(j,n,k,l,i,m) + &
                             tmp_6(j,n,k,l,m,i) + tmp_6(j,n,k,m,i,l) + tmp_6(j,n,k,m,l,i) + &
                             tmp_6(j,n,l,i,k,m) + tmp_6(j,n,l,i,m,k) + tmp_6(j,n,l,k,i,m) + &
                             tmp_6(j,n,l,k,m,i) + tmp_6(j,n,l,m,i,k) + tmp_6(j,n,l,m,k,i) + &
                             tmp_6(j,n,m,i,k,l) + tmp_6(j,n,m,i,l,k) + tmp_6(j,n,m,k,i,l) + &
                             tmp_6(j,n,m,k,l,i) + tmp_6(j,n,m,l,i,k) + tmp_6(j,n,m,l,k,i) + &
                             tmp_6(k,i,j,l,m,n) + tmp_6(k,i,j,l,n,m) + tmp_6(k,i,j,m,l,n) + &
                             tmp_6(k,i,j,m,n,l) + tmp_6(k,i,j,n,l,m) + tmp_6(k,i,j,n,m,l) + &
                             tmp_6(k,i,l,j,m,n) + tmp_6(k,i,l,j,n,m) + tmp_6(k,i,l,m,j,n) + &
                             tmp_6(k,i,l,m,n,j) + tmp_6(k,i,l,n,j,m) + tmp_6(k,i,l,n,m,j) + &
                             tmp_6(k,i,m,j,l,n) + tmp_6(k,i,m,j,n,l) + tmp_6(k,i,m,l,j,n) + &
                             tmp_6(k,i,m,l,n,j) + tmp_6(k,i,m,n,j,l) + tmp_6(k,i,m,n,l,j) + &
                             tmp_6(k,i,n,j,l,m) + tmp_6(k,i,n,j,m,l) + tmp_6(k,i,n,l,j,m) + &
                             tmp_6(k,i,n,l,m,j) + tmp_6(k,i,n,m,j,l) + tmp_6(k,i,n,m,l,j) + &
                             tmp_6(k,j,i,l,m,n) + tmp_6(k,j,i,l,n,m) + tmp_6(k,j,i,m,l,n) + &
                             tmp_6(k,j,i,m,n,l) + tmp_6(k,j,i,n,l,m) + tmp_6(k,j,i,n,m,l) + &
                             tmp_6(k,j,l,i,m,n) + tmp_6(k,j,l,i,n,m) + tmp_6(k,j,l,m,i,n) + &
                             tmp_6(k,j,l,m,n,i) + tmp_6(k,j,l,n,i,m) + tmp_6(k,j,l,n,m,i) + &
                             tmp_6(k,j,m,i,l,n) + tmp_6(k,j,m,i,n,l) + tmp_6(k,j,m,l,i,n) + &
                             tmp_6(k,j,m,l,n,i) + tmp_6(k,j,m,n,i,l) + tmp_6(k,j,m,n,l,i) + &
                             tmp_6(k,j,n,i,l,m) + tmp_6(k,j,n,i,m,l) + tmp_6(k,j,n,l,i,m) + &
                             tmp_6(k,j,n,l,m,i) + tmp_6(k,j,n,m,i,l) + tmp_6(k,j,n,m,l,i) + &
                             tmp_6(k,l,i,j,m,n) + tmp_6(k,l,i,j,n,m) + tmp_6(k,l,i,m,j,n) + &
                             tmp_6(k,l,i,m,n,j) + tmp_6(k,l,i,n,j,m) + tmp_6(k,l,i,n,m,j) + &
                             tmp_6(k,l,j,i,m,n) + tmp_6(k,l,j,i,n,m) + tmp_6(k,l,j,m,i,n) + &
                             tmp_6(k,l,j,m,n,i) + tmp_6(k,l,j,n,i,m) + tmp_6(k,l,j,n,m,i) + &
                             tmp_6(k,l,m,i,j,n) + tmp_6(k,l,m,i,n,j) + tmp_6(k,l,m,j,i,n) + &
                             tmp_6(k,l,m,j,n,i) + tmp_6(k,l,m,n,i,j) + tmp_6(k,l,m,n,j,i) + &
                             tmp_6(k,l,n,i,j,m) + tmp_6(k,l,n,i,m,j) + tmp_6(k,l,n,j,i,m) + &
                             tmp_6(k,l,n,j,m,i) + tmp_6(k,l,n,m,i,j) + tmp_6(k,l,n,m,j,i) + &
                             tmp_6(k,m,i,j,l,n) + tmp_6(k,m,i,j,n,l) + tmp_6(k,m,i,l,j,n) + &
                             tmp_6(k,m,i,l,n,j) + tmp_6(k,m,i,n,j,l) + tmp_6(k,m,i,n,l,j) + &
                             tmp_6(k,m,j,i,l,n) + tmp_6(k,m,j,i,n,l) + tmp_6(k,m,j,l,i,n) + &
                             tmp_6(k,m,j,l,n,i) + tmp_6(k,m,j,n,i,l) + tmp_6(k,m,j,n,l,i) + &
                             tmp_6(k,m,l,i,j,n) + tmp_6(k,m,l,i,n,j) + tmp_6(k,m,l,j,i,n) + &
                             tmp_6(k,m,l,j,n,i) + tmp_6(k,m,l,n,i,j) + tmp_6(k,m,l,n,j,i) + &
                             tmp_6(k,m,n,i,j,l) + tmp_6(k,m,n,i,l,j) + tmp_6(k,m,n,j,i,l) + &
                             tmp_6(k,m,n,j,l,i) + tmp_6(k,m,n,l,i,j) + tmp_6(k,m,n,l,j,i) + &
                             tmp_6(k,n,i,j,l,m) + tmp_6(k,n,i,j,m,l) + tmp_6(k,n,i,l,j,m) + &
                             tmp_6(k,n,i,l,m,j) + tmp_6(k,n,i,m,j,l) + tmp_6(k,n,i,m,l,j) + &
                             tmp_6(k,n,j,i,l,m) + tmp_6(k,n,j,i,m,l) + tmp_6(k,n,j,l,i,m) + &
                             tmp_6(k,n,j,l,m,i) + tmp_6(k,n,j,m,i,l) + tmp_6(k,n,j,m,l,i) + &
                             tmp_6(k,n,l,i,j,m) + tmp_6(k,n,l,i,m,j) + tmp_6(k,n,l,j,i,m) + &
                             tmp_6(k,n,l,j,m,i) + tmp_6(k,n,l,m,i,j) + tmp_6(k,n,l,m,j,i) + &
                             tmp_6(k,n,m,i,j,l) + tmp_6(k,n,m,i,l,j) + tmp_6(k,n,m,j,i,l) + &
                             tmp_6(k,n,m,j,l,i) + tmp_6(k,n,m,l,i,j) + tmp_6(k,n,m,l,j,i) + &
                             tmp_6(l,i,j,k,m,n) + tmp_6(l,i,j,k,n,m) + tmp_6(l,i,j,m,k,n) + &
                             tmp_6(l,i,j,m,n,k) + tmp_6(l,i,j,n,k,m) + tmp_6(l,i,j,n,m,k) + &
                             tmp_6(l,i,k,j,m,n) + tmp_6(l,i,k,j,n,m) + tmp_6(l,i,k,m,j,n) + &
                             tmp_6(l,i,k,m,n,j) + tmp_6(l,i,k,n,j,m) + tmp_6(l,i,k,n,m,j) + &
                             tmp_6(l,i,m,j,k,n) + tmp_6(l,i,m,j,n,k) + tmp_6(l,i,m,k,j,n) + &
                             tmp_6(l,i,m,k,n,j) + tmp_6(l,i,m,n,j,k) + tmp_6(l,i,m,n,k,j) + &
                             tmp_6(l,i,n,j,k,m) + tmp_6(l,i,n,j,m,k) + tmp_6(l,i,n,k,j,m) + &
                             tmp_6(l,i,n,k,m,j) + tmp_6(l,i,n,m,j,k) + tmp_6(l,i,n,m,k,j) + &
                             tmp_6(l,j,i,k,m,n) + tmp_6(l,j,i,k,n,m) + tmp_6(l,j,i,m,k,n) + &
                             tmp_6(l,j,i,m,n,k) + tmp_6(l,j,i,n,k,m) + tmp_6(l,j,i,n,m,k) + &
                             tmp_6(l,j,k,i,m,n) + tmp_6(l,j,k,i,n,m) + tmp_6(l,j,k,m,i,n) + &
                             tmp_6(l,j,k,m,n,i) + tmp_6(l,j,k,n,i,m) + tmp_6(l,j,k,n,m,i) + &
                             tmp_6(l,j,m,i,k,n) + tmp_6(l,j,m,i,n,k) + tmp_6(l,j,m,k,i,n) + &
                             tmp_6(l,j,m,k,n,i) + tmp_6(l,j,m,n,i,k) + tmp_6(l,j,m,n,k,i) + &
                             tmp_6(l,j,n,i,k,m) + tmp_6(l,j,n,i,m,k) + tmp_6(l,j,n,k,i,m) + &
                             tmp_6(l,j,n,k,m,i) + tmp_6(l,j,n,m,i,k) + tmp_6(l,j,n,m,k,i) + &
                             tmp_6(l,k,i,j,m,n) + tmp_6(l,k,i,j,n,m) + tmp_6(l,k,i,m,j,n) + &
                             tmp_6(l,k,i,m,n,j) + tmp_6(l,k,i,n,j,m) + tmp_6(l,k,i,n,m,j) + &
                             tmp_6(l,k,j,i,m,n) + tmp_6(l,k,j,i,n,m) + tmp_6(l,k,j,m,i,n) + &
                             tmp_6(l,k,j,m,n,i) + tmp_6(l,k,j,n,i,m) + tmp_6(l,k,j,n,m,i) + &
                             tmp_6(l,k,m,i,j,n) + tmp_6(l,k,m,i,n,j) + tmp_6(l,k,m,j,i,n) + &
                             tmp_6(l,k,m,j,n,i) + tmp_6(l,k,m,n,i,j) + tmp_6(l,k,m,n,j,i) + &
                             tmp_6(l,k,n,i,j,m) + tmp_6(l,k,n,i,m,j) + tmp_6(l,k,n,j,i,m) + &
                             tmp_6(l,k,n,j,m,i) + tmp_6(l,k,n,m,i,j) + tmp_6(l,k,n,m,j,i) + &
                             tmp_6(l,m,i,j,k,n) + tmp_6(l,m,i,j,n,k) + tmp_6(l,m,i,k,j,n) + &
                             tmp_6(l,m,i,k,n,j) + tmp_6(l,m,i,n,j,k) + tmp_6(l,m,i,n,k,j) + &
                             tmp_6(l,m,j,i,k,n) + tmp_6(l,m,j,i,n,k) + tmp_6(l,m,j,k,i,n) + &
                             tmp_6(l,m,j,k,n,i) + tmp_6(l,m,j,n,i,k) + tmp_6(l,m,j,n,k,i) + &
                             tmp_6(l,m,k,i,j,n) + tmp_6(l,m,k,i,n,j) + tmp_6(l,m,k,j,i,n) + &
                             tmp_6(l,m,k,j,n,i) + tmp_6(l,m,k,n,i,j) + tmp_6(l,m,k,n,j,i) + &
                             tmp_6(l,m,n,i,j,k) + tmp_6(l,m,n,i,k,j) + tmp_6(l,m,n,j,i,k) + &
                             tmp_6(l,m,n,j,k,i) + tmp_6(l,m,n,k,i,j) + tmp_6(l,m,n,k,j,i) + &
                             tmp_6(l,n,i,j,k,m) + tmp_6(l,n,i,j,m,k) + tmp_6(l,n,i,k,j,m) + &
                             tmp_6(l,n,i,k,m,j) + tmp_6(l,n,i,m,j,k) + tmp_6(l,n,i,m,k,j) + &
                             tmp_6(l,n,j,i,k,m) + tmp_6(l,n,j,i,m,k) + tmp_6(l,n,j,k,i,m) + &
                             tmp_6(l,n,j,k,m,i) + tmp_6(l,n,j,m,i,k) + tmp_6(l,n,j,m,k,i) + &
                             tmp_6(l,n,k,i,j,m) + tmp_6(l,n,k,i,m,j) + tmp_6(l,n,k,j,i,m) + &
                             tmp_6(l,n,k,j,m,i) + tmp_6(l,n,k,m,i,j) + tmp_6(l,n,k,m,j,i) + &
                             tmp_6(l,n,m,i,j,k) + tmp_6(l,n,m,i,k,j) + tmp_6(l,n,m,j,i,k) + &
                             tmp_6(l,n,m,j,k,i) + tmp_6(l,n,m,k,i,j) + tmp_6(l,n,m,k,j,i) + &
                             tmp_6(m,i,j,k,l,n) + tmp_6(m,i,j,k,n,l) + tmp_6(m,i,j,l,k,n) + &
                             tmp_6(m,i,j,l,n,k) + tmp_6(m,i,j,n,k,l) + tmp_6(m,i,j,n,l,k) + &
                             tmp_6(m,i,k,j,l,n) + tmp_6(m,i,k,j,n,l) + tmp_6(m,i,k,l,j,n) + &
                             tmp_6(m,i,k,l,n,j) + tmp_6(m,i,k,n,j,l) + tmp_6(m,i,k,n,l,j) + &
                             tmp_6(m,i,l,j,k,n) + tmp_6(m,i,l,j,n,k) + tmp_6(m,i,l,k,j,n) + &
                             tmp_6(m,i,l,k,n,j) + tmp_6(m,i,l,n,j,k) + tmp_6(m,i,l,n,k,j) + &
                             tmp_6(m,i,n,j,k,l) + tmp_6(m,i,n,j,l,k) + tmp_6(m,i,n,k,j,l) + &
                             tmp_6(m,i,n,k,l,j) + tmp_6(m,i,n,l,j,k) + tmp_6(m,i,n,l,k,j) + &
                             tmp_6(m,j,i,k,l,n) + tmp_6(m,j,i,k,n,l) + tmp_6(m,j,i,l,k,n) + &
                             tmp_6(m,j,i,l,n,k) + tmp_6(m,j,i,n,k,l) + tmp_6(m,j,i,n,l,k) + &
                             tmp_6(m,j,k,i,l,n) + tmp_6(m,j,k,i,n,l) + tmp_6(m,j,k,l,i,n) + &
                             tmp_6(m,j,k,l,n,i) + tmp_6(m,j,k,n,i,l) + tmp_6(m,j,k,n,l,i) + &
                             tmp_6(m,j,l,i,k,n) + tmp_6(m,j,l,i,n,k) + tmp_6(m,j,l,k,i,n) + &
                             tmp_6(m,j,l,k,n,i) + tmp_6(m,j,l,n,i,k) + tmp_6(m,j,l,n,k,i) + &
                             tmp_6(m,j,n,i,k,l) + tmp_6(m,j,n,i,l,k) + tmp_6(m,j,n,k,i,l) + &
                             tmp_6(m,j,n,k,l,i) + tmp_6(m,j,n,l,i,k) + tmp_6(m,j,n,l,k,i) + &
                             tmp_6(m,k,i,j,l,n) + tmp_6(m,k,i,j,n,l) + tmp_6(m,k,i,l,j,n) + &
                             tmp_6(m,k,i,l,n,j) + tmp_6(m,k,i,n,j,l) + tmp_6(m,k,i,n,l,j) + &
                             tmp_6(m,k,j,i,l,n) + tmp_6(m,k,j,i,n,l) + tmp_6(m,k,j,l,i,n) + &
                             tmp_6(m,k,j,l,n,i) + tmp_6(m,k,j,n,i,l) + tmp_6(m,k,j,n,l,i) + &
                             tmp_6(m,k,l,i,j,n) + tmp_6(m,k,l,i,n,j) + tmp_6(m,k,l,j,i,n) + &
                             tmp_6(m,k,l,j,n,i) + tmp_6(m,k,l,n,i,j) + tmp_6(m,k,l,n,j,i) + &
                             tmp_6(m,k,n,i,j,l) + tmp_6(m,k,n,i,l,j) + tmp_6(m,k,n,j,i,l) + &
                             tmp_6(m,k,n,j,l,i) + tmp_6(m,k,n,l,i,j) + tmp_6(m,k,n,l,j,i) + &
                             tmp_6(m,l,i,j,k,n) + tmp_6(m,l,i,j,n,k) + tmp_6(m,l,i,k,j,n) + &
                             tmp_6(m,l,i,k,n,j) + tmp_6(m,l,i,n,j,k) + tmp_6(m,l,i,n,k,j) + &
                             tmp_6(m,l,j,i,k,n) + tmp_6(m,l,j,i,n,k) + tmp_6(m,l,j,k,i,n) + &
                             tmp_6(m,l,j,k,n,i) + tmp_6(m,l,j,n,i,k) + tmp_6(m,l,j,n,k,i) + &
                             tmp_6(m,l,k,i,j,n) + tmp_6(m,l,k,i,n,j) + tmp_6(m,l,k,j,i,n) + &
                             tmp_6(m,l,k,j,n,i) + tmp_6(m,l,k,n,i,j) + tmp_6(m,l,k,n,j,i) + &
                             tmp_6(m,l,n,i,j,k) + tmp_6(m,l,n,i,k,j) + tmp_6(m,l,n,j,i,k) + &
                             tmp_6(m,l,n,j,k,i) + tmp_6(m,l,n,k,i,j) + tmp_6(m,l,n,k,j,i) + &
                             tmp_6(m,n,i,j,k,l) + tmp_6(m,n,i,j,l,k) + tmp_6(m,n,i,k,j,l) + &
                             tmp_6(m,n,i,k,l,j) + tmp_6(m,n,i,l,j,k) + tmp_6(m,n,i,l,k,j) + &
                             tmp_6(m,n,j,i,k,l) + tmp_6(m,n,j,i,l,k) + tmp_6(m,n,j,k,i,l) + &
                             tmp_6(m,n,j,k,l,i) + tmp_6(m,n,j,l,i,k) + tmp_6(m,n,j,l,k,i) + &
                             tmp_6(m,n,k,i,j,l) + tmp_6(m,n,k,i,l,j) + tmp_6(m,n,k,j,i,l) + &
                             tmp_6(m,n,k,j,l,i) + tmp_6(m,n,k,l,i,j) + tmp_6(m,n,k,l,j,i) + &
                             tmp_6(m,n,l,i,j,k) + tmp_6(m,n,l,i,k,j) + tmp_6(m,n,l,j,i,k) + &
                             tmp_6(m,n,l,j,k,i) + tmp_6(m,n,l,k,i,j) + tmp_6(m,n,l,k,j,i) + &
                             tmp_6(n,i,j,k,l,m) + tmp_6(n,i,j,k,m,l) + tmp_6(n,i,j,l,k,m) + &
                             tmp_6(n,i,j,l,m,k) + tmp_6(n,i,j,m,k,l) + tmp_6(n,i,j,m,l,k) + &
                             tmp_6(n,i,k,j,l,m) + tmp_6(n,i,k,j,m,l) + tmp_6(n,i,k,l,j,m) + &
                             tmp_6(n,i,k,l,m,j) + tmp_6(n,i,k,m,j,l) + tmp_6(n,i,k,m,l,j) + &
                             tmp_6(n,i,l,j,k,m) + tmp_6(n,i,l,j,m,k) + tmp_6(n,i,l,k,j,m) + &
                             tmp_6(n,i,l,k,m,j) + tmp_6(n,i,l,m,j,k) + tmp_6(n,i,l,m,k,j) + &
                             tmp_6(n,i,m,j,k,l) + tmp_6(n,i,m,j,l,k) + tmp_6(n,i,m,k,j,l) + &
                             tmp_6(n,i,m,k,l,j) + tmp_6(n,i,m,l,j,k) + tmp_6(n,i,m,l,k,j) + &
                             tmp_6(n,j,i,k,l,m) + tmp_6(n,j,i,k,m,l) + tmp_6(n,j,i,l,k,m) + &
                             tmp_6(n,j,i,l,m,k) + tmp_6(n,j,i,m,k,l) + tmp_6(n,j,i,m,l,k) + &
                             tmp_6(n,j,k,i,l,m) + tmp_6(n,j,k,i,m,l) + tmp_6(n,j,k,l,i,m) + &
                             tmp_6(n,j,k,l,m,i) + tmp_6(n,j,k,m,i,l) + tmp_6(n,j,k,m,l,i) + &
                             tmp_6(n,j,l,i,k,m) + tmp_6(n,j,l,i,m,k) + tmp_6(n,j,l,k,i,m) + &
                             tmp_6(n,j,l,k,m,i) + tmp_6(n,j,l,m,i,k) + tmp_6(n,j,l,m,k,i) + &
                             tmp_6(n,j,m,i,k,l) + tmp_6(n,j,m,i,l,k) + tmp_6(n,j,m,k,i,l) + &
                             tmp_6(n,j,m,k,l,i) + tmp_6(n,j,m,l,i,k) + tmp_6(n,j,m,l,k,i) + &
                             tmp_6(n,k,i,j,l,m) + tmp_6(n,k,i,j,m,l) + tmp_6(n,k,i,l,j,m) + &
                             tmp_6(n,k,i,l,m,j) + tmp_6(n,k,i,m,j,l) + tmp_6(n,k,i,m,l,j) + &
                             tmp_6(n,k,j,i,l,m) + tmp_6(n,k,j,i,m,l) + tmp_6(n,k,j,l,i,m) + &
                             tmp_6(n,k,j,l,m,i) + tmp_6(n,k,j,m,i,l) + tmp_6(n,k,j,m,l,i) + &
                             tmp_6(n,k,l,i,j,m) + tmp_6(n,k,l,i,m,j) + tmp_6(n,k,l,j,i,m) + &
                             tmp_6(n,k,l,j,m,i) + tmp_6(n,k,l,m,i,j) + tmp_6(n,k,l,m,j,i) + &
                             tmp_6(n,k,m,i,j,l) + tmp_6(n,k,m,i,l,j) + tmp_6(n,k,m,j,i,l) + &
                             tmp_6(n,k,m,j,l,i) + tmp_6(n,k,m,l,i,j) + tmp_6(n,k,m,l,j,i) + &
                             tmp_6(n,l,i,j,k,m) + tmp_6(n,l,i,j,m,k) + tmp_6(n,l,i,k,j,m) + &
                             tmp_6(n,l,i,k,m,j) + tmp_6(n,l,i,m,j,k) + tmp_6(n,l,i,m,k,j) + &
                             tmp_6(n,l,j,i,k,m) + tmp_6(n,l,j,i,m,k) + tmp_6(n,l,j,k,i,m) + &
                             tmp_6(n,l,j,k,m,i) + tmp_6(n,l,j,m,i,k) + tmp_6(n,l,j,m,k,i) + &
                             tmp_6(n,l,k,i,j,m) + tmp_6(n,l,k,i,m,j) + tmp_6(n,l,k,j,i,m) + &
                             tmp_6(n,l,k,j,m,i) + tmp_6(n,l,k,m,i,j) + tmp_6(n,l,k,m,j,i) + &
                             tmp_6(n,l,m,i,j,k) + tmp_6(n,l,m,i,k,j) + tmp_6(n,l,m,j,i,k) + &
                             tmp_6(n,l,m,j,k,i) + tmp_6(n,l,m,k,i,j) + tmp_6(n,l,m,k,j,i) + &
                             tmp_6(n,m,i,j,k,l) + tmp_6(n,m,i,j,l,k) + tmp_6(n,m,i,k,j,l) + &
                             tmp_6(n,m,i,k,l,j) + tmp_6(n,m,i,l,j,k) + tmp_6(n,m,i,l,k,j) + &
                             tmp_6(n,m,j,i,k,l) + tmp_6(n,m,j,i,l,k) + tmp_6(n,m,j,k,i,l) + &
                             tmp_6(n,m,j,k,l,i) + tmp_6(n,m,j,l,i,k) + tmp_6(n,m,j,l,k,i) + &
                             tmp_6(n,m,k,i,j,l) + tmp_6(n,m,k,i,l,j) + tmp_6(n,m,k,j,i,l) + &
                             tmp_6(n,m,k,j,l,i) + tmp_6(n,m,k,l,i,j) + tmp_6(n,m,k,l,j,i) + &
                             tmp_6(n,m,l,i,j,k) + tmp_6(n,m,l,i,k,j) + tmp_6(n,m,l,j,i,k) + &
                             tmp_6(n,m,l,j,k,i) + tmp_6(n,m,l,k,i,j) + tmp_6(n,m,l,k,j,i)
                             
                         tmp_6(i,j,k,l,m,n) = r

                      end do
                   end do
                end do
             end do
          end do
       end do

       h = 0
       do i = 1, ncor
          do j = i, ncor
             do k = j, ncor
                do l = k, ncor
                   do m = l, ncor
                      do n = m, ncor
                         h = h + 1
                         ! MR: UNSURE ABOUT ORIGINS OF FACTOR 2 IN NEXT LINE
                         ave(h) = 2 * tmp_6(i,j,k,l,m,n)
                      end do
                   end do
                end do
             end do
          end do
       end do

       deallocate(tmp_6)


    else
       print *, 'rsp_twoave_tr error: Contribution not implemented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('rsp_twoave_tr error: Contribution not implemented or in wrong order')
    end if

  end if

  end subroutine



  !> Compute differentiated overlap matrices, and optionally
  !> add half-differentiated overlap contribution to Fock matrices
  subroutine rsp_ovlint(nr_ao, nf, f, c, nc, ovl, w, fock)
    !> structure containing integral program settings
    integer,       intent(in)    :: nr_ao
    !> number of fields
    integer,       intent(in)    :: nf
    !> field labels in std order
    character(4),  intent(in)    :: f(nf)
    !> first and number of- components in each field
    integer,       intent(in)    :: c(nf), nc(nf)
    !> resulting overlap integral matrices (incoming content deleted)
    type(matrix),  intent(inout) :: ovl(product(nc))
    !> frequencies of each field
    complex(8),    intent(in),    optional :: w(nf)
    !> Fock matrices to which the half-differentiated overlap
    !> contribution is ADDED
    type(matrix),  intent(inout), optional :: fock(product(nc))
    call interface_1el_ovlint(nr_ao, nf, f, c, nc, ovl, w, fock)
  end subroutine


  !> Compute differentiated overlap matrices,
  subroutine rsp_ovlint_tr(nr_ao, nf, f, c, nc, nblks, blk_info, & 
                                      blk_sizes, propsize, ovl)
    !> number of fields
    integer,       intent(in)    :: nf, propsize
    integer :: nblks
    integer, dimension(nblks) :: blk_sizes
    integer, dimension(nblks, 3) :: blk_info
    !> field labels in std order
    character(4),  intent(in)    :: f(nf)
    !> first and number of- components in each field
    integer,       intent(in)    :: c(nf), nc(nf)
    !> resulting overlap integral matrices (incoming content deleted)
    type(matrix),  intent(inout) :: ovl(propsize)
    !------------------------------------------------
    integer      i, nr_ao

    call interface_1el_ovlint_tr(nr_ao, nf, f, c, nc, nblks, blk_info, & 
                                 blk_sizes, propsize, ovl = ovl)

  end subroutine



!   ! MaR: This routine not tested - awaiting development in integral code
!   !> Compute half-differentiated overlap contribution to Fock matrices
  recursive subroutine rsp_ovlint_t_matrix(nr_ao, num_fields, fields, bra, &
                                           ket, propsize, fock)

    implicit none

    integer :: nr_ao, num_fields, propsize, i, j, k, under_construction, merged_nblks
    integer :: fock_offset, int_result_offset, merged_triang_size, tmp_fock_size
    integer :: total_num_perturbations
    type(p_tuple) :: fields, bra, ket, merged_p_tuple
    integer, dimension(bra%n_perturbations + ket%n_perturbations) :: pids_current_contribution
    type(matrix), dimension(propsize) :: fock, tmp_fock
    integer, allocatable, dimension(:) :: nfields, nblks_tuple, blks_tuple_triang_size, &
                                          blk_sizes_merged, translated_index
    integer, allocatable, dimension(:,:) :: blk_sizes, merged_indices
    integer, allocatable, dimension(:,:,:) :: blks_tuple_info, merged_blk_info


under_construction = 1

if (under_construction == 1) then

write(*,*) 'Called for T matrix contribution - currently unfinished so nothing was added'

else
    

    if (num_fields > 0) then

       call rsp_ovlint_t_matrix(nr_ao, num_fields - 1, p_tuple_remove_first(fields), &
            merge_p_tuple(bra, p_tuple_getone(fields, 1)), ket, propsize, fock)

       call rsp_ovlint_t_matrix(nr_ao, num_fields - 1, p_tuple_remove_first(fields), &
            bra, merge_p_tuple(ket, p_tuple_getone(fields, 1)), propsize, fock)

    else

       do i = 1, propsize
          
          tmp_fock(i) = 0 * fock(1)
          call mat_ensure_alloc(tmp_fock(i))

       end do

       merged_p_tuple = merge_p_tuple(bra, ket)

       ! Make frequency independent blocks for bra/ket

       bra%freq = 0.0
       ket%freq = 0.0

       allocate(nfields(2))
       allocate(nblks_tuple(2))
    
       nfields(1) = bra%n_perturbations
       nblks_tuple(1) = get_num_blks(bra)

       nfields(2) = ket%n_perturbations
       nblks_tuple(2) = get_num_blks(ket)
               
       total_num_perturbations = sum(nfields)

       allocate(blks_tuple_info(2, total_num_perturbations, 3))
       allocate(blks_tuple_triang_size(2))
       allocate(blk_sizes(2, total_num_perturbations))
       allocate(blk_sizes_merged(total_num_perturbations))
    
       blks_tuple_info(1, :, :) = get_blk_info(nblks_tuple(1), bra)
       blks_tuple_triang_size(1) = get_triangulated_size(nblks_tuple(1), &
                                   blks_tuple_info(1, 1:nblks_tuple(1), :))
       blk_sizes(1, 1:nblks_tuple(1)) = get_triangular_sizes(nblks_tuple(1), &
       blks_tuple_info(1,1:nblks_tuple(1),2), blks_tuple_info(1,1:nblks_tuple(1),3))

       blks_tuple_info(2, :, :) = get_blk_info(nblks_tuple(2), ket)
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

       call interface_1el_ovlint_half_diff(nr_ao, bra%n_perturbations, bra%plab, &
            (/(j/j, j = 1, bra%n_perturbations)/), bra%pdim, ket%n_perturbations, ket%plab, &
            (/(j/j, j = 1, ket%n_perturbations)/), ket%pdim, nblks_tuple, &
            blks_tuple_info, blk_sizes, tmp_fock_size, tmp_fock)

       k = 1
       do j = 1, bra%n_perturbations
          pids_current_contribution(k) = bra%pid(j)
          k = k + 1
       end do

       do j = 1, ket%n_perturbations
          pids_current_contribution(k) = ket%pid(j)
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
    
          if (bra%n_perturbations == 0) then
    
             int_result_offset = get_triang_blks_tuple_offset(1, &
                                 total_num_perturbations, nblks_tuple(2), &
                                 nfields(2), blks_tuple_info(2, :, :), &
                                 blk_sizes(2,:), blks_tuple_triang_size(2), & 
                                 (/ translated_index(:) /))
    
          elseif (ket%n_perturbations == 0) then
    
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
    
          fock(fock_offset) = fock(fock_offset) + &
          0.5 * (sum(bra%freq) - sum(ket%freq)) * tmp_fock(int_result_offset)

       end do

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




  subroutine rsp_oneint(nr_ao, nf, f, c, nc, oneint)
    !> structure containing integral program settings
    integer,       intent(in)    :: nr_ao
    !> number of fields
    integer,       intent(in)    :: nf
    !> field labels in std order
    character(4),  intent(in)    :: f(nf)
    !> first and number of- components in each field
    integer,       intent(in)    :: c(nf), nc(nf)
    !> output perturbed integrals
    type(matrix),  intent(inout) :: oneint(product(nc))
    call interface_1el_oneint(nr_ao, nf, f, c, nc, oneint)
  end subroutine


  subroutine rsp_oneint_tr(nr_ao, nf, f, c, nc, nblks, blk_info, & 
                                      blk_sizes, propsize, oneint)
    !> number of fields
    integer,       intent(in)    :: nf, propsize
    integer :: nblks
    integer, dimension(nblks) :: blk_sizes
    integer, dimension(nblks, 3) :: blk_info
    !> field labels in std order
    character(4),  intent(in)    :: f(nf)
    !> first and number of- components in each field
    integer,       intent(in)    :: c(nf), nc(nf)
    !> output perturbed integrals
    type(matrix),  intent(inout) :: oneint(propsize)
    !--------------------------------------------------
    integer order_mom  !order of Cartesian multipole moments
    integer num_mom    !number of Cartesian multipole moments
    integer order_geo  !order of total geometric derivatives
    integer num_atom   !number of atoms
    integer num_coord  !number of atomic coordinates
    integer num_geom   !number of total geometric derivatives
    integer num_ints   !number of all integral matrices
    integer imat       !incremental recorder over matrices
    integer :: i, nr_ao
    type(matrix) :: A

    call interface_1el_oneint_tr(nr_ao, nf, f, c, nc, nblks, blk_info, & 
                                      blk_sizes, propsize, oneint)

  end subroutine


  !> Contract 2-electron integrals perturbed by fields 'f' with density
  !> matrix 'dens', and add to Fock matrices 'fock' Average 2-electron integrals perturbed by fields f over the
  !> product of density matrces D1 and D2
  subroutine rsp_twoint(nr_ao, nf, f, c, nc, dens, fock)
    ! work array to be passed to GRCONT
    use eri_contractions, only: ctr_arg
    use eri_basis_loops,  only: unopt_geodiff_loop
    use interface_interest
    integer,              intent(in)    :: nr_ao
    !> number of fields
    integer,              intent(in)    :: nf
    !> field labels in std order
    character(4),         intent(in)    :: f(nf)
    !> first and number of- components in each field
    integer,              intent(in)    :: c(nf), nc(nf)
    !> density matrix to average over
    type(matrix), target, intent(in)    :: dens
    !> Fock matrix to which two-electron contribution is ADDED
    type(matrix), target, intent(inout) :: fock(product(nc))
    !--------------------------------------------------
    real(8), pointer :: null_ptr(:) !because null() isn't f90
    integer       i, j, n, ij, ncor
    type(ctr_arg) arg(1)
    type(matrix)  A !scratch
    real(8)          :: dummy(1)

    nullify(null_ptr) !because null() isn't f90

  if (any(f=='EL  ')) then
     call mat_init(A, nr_ao, nr_ao, &
                   .false., .false., .false., .false., .false.)
     do i = 1, product(nc)
        if (iszero(fock(i))) then
           call mat_ensure_alloc(fock(i))
           fock(i)%elms = fock(i)%elms + A%elms
        else
           fock(i)%elms = fock(i)%elms + A%elms
        end if
     end do
  else

    if (nf==0) then
       A = 0*dens
       call mat_ensure_alloc(A)
       call interface_scf_get_g(dens, A)
       fock(1) = fock(1) + A
       A = 0
    else if (nf==1 .and. f(1)=='GEO ') then

#ifdef PRG_DALTON
       n = nr_ao
       ncor = 3 * get_nr_atoms()
       do i = 0, nc(1)-1
          if (iszero(fock(i+1))) then
             call mat_ensure_alloc(fock(i+1))
          end if
#ifdef GRCONT_NOT_AVAILABLE
          arg(1) = ctr_arg(1, i+1, &
                           ncor, dens, fock(i+1), null_ptr)
          call unopt_geodiff_loop(basis_large, &
                                  basis_small, &
                                  arg)
#else
          ! if first or an x-coord, call GRCONT
          if (i==0 .or. mod(c(1)+i,3) == 1) then
             call GRCONT(f77_memory(n*n*3+1:), size(f77_memory)-n*n*3, &
                         f77_memory(:n*n*3), n*n*3, .true., .false., &
                         1, (c(1)+i+2)/3, .false., .true., dens%elms, 1)
          end if
          j = 1 + mod(c(1)+i-1,3) !x y z = 1 2 3
          if (iszero(fock(1+i))) then
             fock(1+i)%elms(:, :, 1) = reshape(f77_memory(n*n*(j-1)+1:n*n*j),(/n,n/))
          else
             fock(1+i)%elms(:, :, 1) = fock(1+i)%elms(:, :, 1) &
                            + reshape(f77_memory(n*n*(j-1)+1:n*n*j),(/n,n/))
          end if
#endif
       end do
#endif /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
       do i = 1, nc(1)
          if (iszero(fock(i))) then
             call mat_ensure_alloc(fock(i))
          end if
       end do
       do i = 1, nc(1)
          call interest_mpi_wake_up()
          call interest_get_int(dens%nrow, dens%elms, fock(i)%elms, 1, i, 0, dummy)
       end do
#endif /* ifdef PRG_DIRAC */

    else if (nf==2 .and. all(f==(/'GEO ','GEO '/))) then

#ifdef PRG_DALTON
       ncor = 3 * get_nr_atoms()
       do j = 0, nc(2)-1
          do i = 0, nc(1)-1
             ij = 1 + i + nc(1)*j
             if (iszero(fock(ij))) then
                call mat_ensure_alloc(fock(ij))
                fock(ij)%elms = 0 !ajt FIXME use mat_axpy
             end if
             arg(1) = ctr_arg(2, c(1)+i + ncor * (c(2)+j-1), &
                              ncor, dens, fock(ij), null_ptr)
             call unopt_geodiff_loop(basis_large, &
                                     basis_small, &
                                     arg)
          end do
       end do
#endif /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
       print *, 'error: the 2nd order twoint contribution is not available in DIRAC'
       stop 1
#endif /* ifdef PRG_DIRAC */

    else
       print *, 'error in rsp_oneave: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('error in rsp_oneave: not implented or in wrong order')
    end if

  end if

  end subroutine


! MR: SHELLS_NUCLEI_displace copied from old file. Just used temporarily 
! for testing 6th order force field (the (/'GEO ', 'GEO ', 'GEO '/) contribution
! in twoint is not implememted analytically yet).

  !> Move the nuclei and basis functions, as seen from the integral
  !> programs. For doing finite difference differentiation.
  subroutine SHELLS_NUCLEI_displace(ic, dc)
    implicit integer (i,m-n)
#include <implicit.h>
    integer,  intent(in) :: ic
    real(8), intent(in) :: dc !step
    ! need MXSHEL
#include <maxorb.h>
    ! need CENT
#include <shells.h>
    ! need MXCOOR
#include <mxcent.h>
    ! need CORD
#include <nuclei.h>
    integer a, r, i
    a = 1 + (ic-1)/3
    r = 1 + mod(ic-1,3)
    if (a <= 0 .or. a > NATOMS) &
       call quit("SHELLS_NUCLEI_displace error: arg. 'ic' out of range")
    ! move atom in CORD, then in CENT
    CORD(r,a)  = CORD(r,a)  + dc
    do i = 1, NLRGSH
       if (NCENT(i) == a) &
          CENT(i,r,:) = CENT(i,r,:) + dc
    end do
  end subroutine

  !> Contract 2-electron integrals perturbed by fields 'f' with density
  !> matrix 'dens', and add to Fock matrices 'fock' Average 2-electron integrals perturbed by fields f over the
  !> product of density matrces D1 and D2
  subroutine rsp_twoint_tr(nr_ao, nf, f, c, nc, dens, propsize, fock)
    use eri_contractions, only: ctr_arg
    use eri_basis_loops,  only: unopt_geodiff_loop
    use interface_interest
    use rsp_indices_and_addressing, only: make_triangulated_indices
    !> number of fields
    integer,              intent(in)    :: nf, propsize
    !> field labels in std order
    character(4),         intent(in)    :: f(nf)
    !> first and number of- components in each field
    integer,              intent(in)    :: c(nf), nc(nf)
    !> density matrix to average over
    type(matrix), target, intent(in)    :: dens
    !> Fock matrix to which two-electron contribution is ADDED
    type(matrix), target, intent(inout) :: fock(propsize)
    type(matrix), allocatable, target, dimension(:) :: tmpfock
    !--------------------------------------------------
    real(8), pointer :: null_ptr(:) !because null() isn't f90
    real(8)  :: fdistep = 2d0**(-25)
    integer       h, hincr, i, j, k, n, ij, ijk, incr, ncor, nr_ao
    type(ctr_arg) arg(1)
    type(matrix)  A !scratch
    real(8)          :: dummy(1)
    !--------------------------------------------------
    integer, allocatable, dimension(:,:) :: indices

    if (any(f=='EL  ')) then
       ! nothing to add
    else

#ifdef PRG_DALTON
! Begin new general geo code


       nullify(null_ptr)

       if (nf==0) then
       call mat_init(A, fock(1)%nrow, fock(1)%ncol, &
                     .false., .false., .false., .false., .false.)
          A = 0*dens
          call mat_ensure_alloc(A)
! write(*,*) 'A', A%elms
! write(*,*) 'dens', dens%elms

          call interface_scf_get_g(dens, A)

! write(*,*) 'A 2', A%elms
          fock(1) = fock(1) + A
          A = 0

       else

          if ((nf == count(f == 'GEO '))) then

             ncor = 3 * get_nr_atoms()
             allocate(indices(propsize, nf))
             call make_triangulated_indices(1, (/1, nf, nc(1)/), propsize, indices)

             do i = 1, propsize

                if (iszero(fock(i))) then
                   call mat_ensure_alloc(fock(i))
                   fock(i)%elms = 0 !ajt FIXME use mat_axpy
                end if
             
                k = 1

                do j = 1, nf

                   k = k + (indices(i,j) - 1) * ( nc(1)**(nf - j ) )

                end do

                   arg(1) = ctr_arg(nf, k, ncor, dens, fock(i), null_ptr)
                   call unopt_geodiff_loop(basis_large, basis_small, arg)

             end do

             deallocate(indices)

          end if

       end if


! End new general geo code

! MaR: I don't understand all the program-specific ifdefs below, so I kept this
! non-general code for "not-Dalton" runs

#else
       if (nf==0) then
          A = 0*dens
          call mat_ensure_alloc(A)
          call interface_scf_get_g(dens, A)
          fock(1) = fock(1) + A
          A = 0
       else if (nf==1 .and. f(1)=='GEO ') then
          n = nr_ao
          do i = 0, nc(1)-1
             if (iszero(fock(i+1))) then
                call mat_ensure_alloc(fock(i+1))
             end if
   
#ifdef PRG_DALTON
#ifdef GRCONT_NOT_AVAILABLE
             arg(1) = ctr_arg(1, i+1, &
                              ncor, dens, fock(i+1), null_ptr)
             call unopt_geodiff_loop(basis_large, &
                                     basis_small, &
                                     arg)
#else
             ! if first or an x-coord, call GRCONT
             if (i==0 .or. mod(c(1)+i,3) == 1) then
                call GRCONT(f77_memory(n*n*3+1:), size(f77_memory)-n*n*3, &
                            f77_memory(:n*n*3), n*n*3, .true., .false., &
                            1, (c(1)+i+2)/3, .false., .true., dens%elms, 1)
             end if
             j = 1 + mod(c(1)+i-1,3) !x y z = 1 2 3
             if (iszero(fock(1+i))) then
                fock(1+i)%elms(:, :, 1) = reshape(f77_memory(n*n*(j-1)+1:n*n*j),(/n,n/))
             else
                fock(1+i)%elms(:, :, 1) = fock(1+i)%elms(:, :, 1) &
                               + reshape(f77_memory(n*n*(j-1)+1:n*n*j),(/n,n/))
             end if
#endif
#endif /* ifdef PRG_DALTON */
   
#ifdef PRG_DIRAC
             call interest_mpi_wake_up()
             call interest_get_int(dens%nrow, dens%elms, fock(i+1)%elms, 1, i+1, 0, dummy)
#endif /* ifdef PRG_DIRAC */
   
          end do
       else if (nf==2 .and. all(f==(/'GEO ','GEO '/))) then
          ncor = 3 * get_nr_atoms()
          ij = 0
          do j = 0, nc(2)-1
             do i = j, nc(1)-1
                ij = ij + 1
                if (iszero(fock(ij))) then
                   call mat_ensure_alloc(fock(ij))
                   fock(ij)%elms = 0 !ajt FIXME use mat_axpy
                end if
                arg(1) = ctr_arg(2, c(1)+i + ncor * (c(2)+j-1), &
                                 ncor, dens, fock(ij), null_ptr)
                call unopt_geodiff_loop(basis_large, &
                                        basis_small, &
                                        arg)
             end do
          end do
   
       else if (nf==3 .and. all(f==(/'GEO ','GEO ', 'GEO '/))) then
          ncor = 3 * get_nr_atoms()
          ij = 0
          do k = 0, nc(3)-1
             do j = k, nc(2)-1
                do i = j, nc(1)-1
                   ij = ij + 1
                   if (iszero(fock(ij))) then
                      call mat_ensure_alloc(fock(ij))
                      fock(ij)%elms = 0 !ajt FIXME use mat_axpy
                   end if
                   arg(1) = ctr_arg(3, c(1)+i + ncor * (c(2)+j-1) + &
                                       ncor * ncor * (c(3)+k-1), &
                                       ncor, dens, fock(ij), null_ptr)
                   call unopt_geodiff_loop(basis_large, &
                                           basis_small, &
                                           arg)
                end do
             end do
          end do


   
       else
          print *, 'error in rsp_twoint_tr: not implemented or in wrong order - ', &
                  (' ' // f(i), i=1,nf)
          call quit('error in rsp_twoint_tr: not implemented or in wrong order')
       end if
#endif

    end if

  end subroutine


  subroutine rsp_xcint_tr_adapt(nr_ao, nf, f, c, nc, D, propsize, xcint)

    integer :: i, j, k, nr_ao, propsize
    !> number of fields
    integer,       intent(in)    :: nf
    !> field labels in std order
    character(4),  intent(in)    :: f(nf)
    !> first and number of- components in each field
    integer,       intent(in)    :: c(nf), nc(nf)
    !> output perturbed integrals
    type(matrix),  intent(inout) :: xcint(propsize)
type(matrix) :: Db, Fb
    type(matrix), allocatable, dimension(:,:) :: tmp_xcint
    type(matrix), dimension(:) :: D
    !--------------------------------------------------

    ! MaR: ONLY GEOMETRICAL PERTURBATIONS SUPPORTED SO FAR

! call mat_init(Db, nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
! call mat_init(Fb, nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)

    if (nf == 0) then

! write(*,*) 'D before', D(1)%elms
! write(*,*) 'F before', xcint(1)%elms
! 
! Db = D(1)
! Fb = xcint(1)

      call rsp_xcint(D, F=xcint(1))

! Db = D(1) - Db
! Fb = xcint(1) - Fb
! 
! write(*,*) 'D chg', Db%elms
! write(*,*) 'F after', Fb%elms

    else if (nf == 1) then

       if (all(f==(/'GEO '/))) then

          call rsp_xcint(D, Fg=xcint)

       else

!           write(*,*) 'WARNING (rsp_xcint_tr_adapt): UNSUPPORTED CONTRIBUTION:'
!           write(*,*) 'NO CONTRIBUTION WILL BE MADE'

       end if

    else if (nf == 2) then

       if (all(f==(/'GEO ','GEO '/))) then

          allocate(tmp_xcint(nc(1), nc(1)))

          do i = 1, nc(1)
             do j = 1, nc(1)

                ! MR: ASSUME CLOSED SHELL
                call mat_init(tmp_xcint(i,j), nr_ao, nr_ao, &
                              .false., .false., .false., .false., .false.)

             end do
          end do

          call rsp_xcint(D, Fgg=tmp_xcint)


          do i = 1, nc(1)
             do j = 1, i

                k = get_triang_blks_offset(1, 2, (/1, 2, nc(1)/), &
                                           (/propsize/), (/i, j/))

                xcint(k) = xcint(k) + tmp_xcint(i,j)
                tmp_xcint(i, j) = 0

             end do
          end do

          deallocate(tmp_xcint)

       else

!           write(*,*) 'WARNING: UNSUPPORTED CONTRIBUTION: NO CONTRIBUTION WILL BE MADE'

       end if

    else

!        write(*,*) 'WARNING: UNSUPPORTED NUMBER OF FIELDS: NO CONTRIBUTION WILL BE MADE'

    end if

  end subroutine



  function index_of_field(f) result(i)
    character(4) :: f
    integer      :: i
    do i = 1, size(all_known_fields)
        if (all_known_fields(i)%label == f) return
    end do
    call quit('Field not found: ' // f)
  end function


  function rsp_field_anti(f)
    character(4), intent(in) :: f(:)
    logical :: rsp_field_anti(size(f))
    integer :: i
    rsp_field_anti = (/(all_known_fields(index_of_field(f(i)))%anti, i=1,size(f))/)
  end function


  !> shape (dimensions) of property for fields f(:)
  function rsp_field_dim(f)
    character(4),  intent(in) :: f(:)
    integer :: rsp_field_dim(size(f)), i
    rsp_field_dim = (/(all_known_fields(index_of_field(f(i)))%ncomp, i=1,size(f))/)
    ! loop through mol-dependent
    do i = 1, size(f)
       if (rsp_field_dim(i) /= -1) then
          ! cycle
       else if (f(i) == 'GEO ') then
          rsp_field_dim(i) = 3 * get_nr_atoms()
       else
          call quit('rsp_field_dim error: Number of comp. unknown for ' // f(i))
       end if
    end do
  end function


  !> which fields are basis-perturbing
  function rsp_field_bas(f)
    character(4), intent(in) :: f(:)
    logical :: rsp_field_bas(size(f))
    integer :: i
    rsp_field_bas = (/(all_known_fields(index_of_field(f(i)))%bas, i=1,size(f))/)
  end function


  function prefix_zeros(n, l)
    integer, intent(in) :: n, l !number, length
    character(l)        :: prefix_zeros !resulting n in ascii
    character(1), parameter :: char0to9(0:9) &
          = (/'0','1','2','3','4','5','6','7','8','9'/)
    integer :: i, k
    k = n
    do i = l, 1, -1
       prefix_zeros(i:i) = char0to9(mod(k,10))
       k = k / 10
    end do
    if (k /= 0) call quit('prefix_zeros error: Argument integer does not fit ' &
                       // 'in the specified number of ASCII caracters')
  end function



  ! determine 'canonical' ordering 'o' of fields 'f'. Also, validate fields%comp/ncomp
  subroutine count_and_prepare(f, geo_first, bas_first, o, tot_ncomp, &
                               nucpot, overlap, hamilt)
    type(rsp_field),   intent(in)  :: f(:)
    logical,           intent(in)  :: geo_first
    logical,           intent(in)  :: bas_first
    integer,           intent(out) :: o(size(f))
    integer,           intent(out) :: tot_ncomp(size(f))
    logical, optional, intent(out) :: nucpot, overlap, hamilt
    integer i, j, k, jj, kk, nel, nmag, nbas, ngeo, fld_idx(size(f))
    ! look up field indices
    fld_idx(:) = (/(index_of_field(f(i)%label), i=1,size(f))/)
    ! initialize order to unity
    o(:) = (/(i,i=1,size(f))/)
    ! selection sort, select field to be in i'th place from i:size(f)
    do i = 1, size(f)
       jj = i
       j  = o(i)
       do kk = i+1, size(f)
          k = o(kk)
          ! firstly, if geo_first, select any 'GEO ' ahead of others
          if (geo_first .and. f(k)%label /= 'GEO ' &
                        .and. f(j)%label == 'GEO ') cycle
          if (.not.geo_first .or. (f(k)%label == 'GEO ' &
                             .eqv. f(j)%label == 'GEO ')) then
             ! secondly, if bas_first, those with field_stats%bas
             if (bas_first .and. .not.all_known_fields(fld_idx(k))%bas &
                           .and.      all_known_fields(fld_idx(j))%bas) cycle
             if (.not.bas_first .or. (all_known_fields(fld_idx(k))%bas &
                                .eqv. all_known_fields(fld_idx(j))%bas)) then
                ! thirdly, pick highest fld_idx
                if (fld_idx(k) <  fld_idx(j)) cycle
                if (fld_idx(k) == fld_idx(j)) then
                   ! fourthly, choose field with highest ncomp
                   if (f(k)%ncomp <  f(j)%ncomp) cycle
                   if (f(k)%ncomp == f(j)%ncomp) then
                      ! fifthly, lowest comp
                      if (f(k)%comp >  f(j)%comp) cycle
                      if (f(k)%comp == f(j)%comp) then
                         ! sixthly, highest abs(re(freq))
                         if (abs(real(f(k)%freq)) <  abs(real(f(j)%freq))) cycle
                         if (abs(real(f(k)%freq)) == abs(real(f(j)%freq))) then
                            ! seventhly, lowest sign(re(freq))
                            if (real(f(k)%freq) >  real(f(j)%freq)) cycle
                            if (real(f(k)%freq) == real(f(j)%freq)) then
                               ! eigthly, highest abs(im(freq))
                               if (abs(aimag(f(k)%freq)) <  abs(aimag(f(j)%freq))) cycle
                               if (abs(aimag(f(k)%freq)) == abs(aimag(f(j)%freq))) then
                                  ! ninethly, lowest sign(im(freq))
                                  if (aimag(f(k)%freq) >  aimag(f(j)%freq)) cycle
                                  if (aimag(f(k)%freq) == aimag(f(j)%freq)) then
                                     ! tenthly, if fields j and k are *identical*, go for the
                                     ! one with lowest input position
                                     if (k > j) cycle
                                  end if
                               end if
                            end if
                         end if
                      end if
                   end if
                end if
             end if
          end if
          jj = kk
          j  = k  !new minimum
       end do
       !swap entries i and j
       k     = o(i)
       o(i)  = j
       o(jj) = k
    end do
    ! place the full number of components in tot_ncomp
    do i = 1, size(f)
       tot_ncomp(i) = all_known_fields(fld_idx(o(i)))%ncomp
       ! ajt FIXME -1 means ncomp is uninitialized (molecule-dependent)
       !           This init should be done in a separate routine, together
       !           with all the other fields' ncomp
       if (tot_ncomp(i) == -1 .and. f(o(i))%label == 'GEO ') then
          tot_ncomp(i) = 3*get_nr_atoms()
          all_known_fields(fld_idx(o(i)))%ncomp = tot_ncomp(i)
       end if
    end do
    ! if requested, determine whether the nuclear potential, the overlap
    ! integrals, or the Hamiltonian, are preturbed by this tuple of fields
    if (present(nucpot)) &
       ngeo = count(f(:)%label == 'GEO ')
    if (present(overlap) .or. present(hamilt)) &
       nbas = count((/(all_known_fields(fld_idx(i))%bas, i=1,size(f))/))
    if (present(nucpot) .or. present(hamilt)) then
       nel  = count((/(all_known_fields(fld_idx(i))%lin,  i=1,size(f))/))
       nmag = count((/(all_known_fields(fld_idx(i))%quad, i=1,size(f))/))
    end if
    if (present(nucpot)) &
       nucpot = (((nel <= 1 .and. nmag == 0) &
             .or. (nel == 0 .and. nmag <= 2)) &
            .and. ngeo + nel + nmag == size(f))
    if (present(overlap)) &
       overlap = (nbas == size(f))
    if (present(hamilt)) &
       nucpot = (((nel <= 1 .and. nmag == 0) &
             .or. (nel == 0 .and. nmag <= 2)) &
            .and. nbas + nel + nmag == size(f))
  end subroutine



  function rsp_field_ordering(f) result(o)
    type(rsp_field), intent(in) :: f(:)
    integer                     :: o(size(f))
    integer tc(size(f))
    call count_and_prepare(f, .false., .false., o, tc)
  end function


  !> Add the average contribution 'ave', with dimensions 'dima' to response tensor
  !> 'rsp', with dimensions 'dimr', but where dimensions are permuted by 'perm'.
  !> 'idxr' selects (if >0) which starting components in 'ave' for each dim in 'rsp',
  !> or (if <0) index of a density matrix dimension.
  !> @prefac  sign or prefactor: rsp = rsp + prefac * ave
  !> @ndim    number of dimensions
  !> @perm    ordering of dimensions. For each dim in ave, the corresponding dim in rsp
  !> @idxr    starting component in ave for each dim in rsp, or if negative, the index
  !>          of a density dimension in rsp
  !> @dima    dimensions of ave. Density dimensions must have dim 1
  !> @dimr    dimensions of rsp
  !> @ave     real array of integral averages, as from integral program
  !> @rsp     complex respons function tensor
  subroutine permute_selcomp_add(prefac, ndim, perm, idxr, dima, dimr, ave, rsp)
    complex(8), intent(in)    :: prefac
    integer,    intent(in)    :: ndim, perm(ndim), idxr(ndim)
    integer,    intent(in)    :: dima(ndim), dimr(ndim)
    real(8),    intent(in)    :: ave(product(dima))
    complex(8), intent(inout) :: rsp(product(dimr))
    integer i, ia, ir, stpr(ndim), stpa(ndim), ii(ndim), dd(ndim)
    ! calculate dimension steps in ave (cumulative products of
    ! dimensions), as well as offset due to starting indices idxr (ia)
    ia = 0
    do i = 1, ndim
       if (i==1) stpa(i) = 1
       if (i/=1) stpa(i) = stpa(i-1) * dima(i-1)
       if (idxr(perm(i)) > 0) & !positive means starting comp
          ia = ia + stpa(i) * (idxr(perm(i)) - 1)
    end do
    ! calculate (permuted) dimension steps in rsp, and offset due to
    ! density indices (ir), and permuted dimensions (dd)
    ir = 0
    do i = 1, ndim
       if (i==1) stpr(perm(i)) = 1
       if (i/=1) stpr(perm(i)) = stpr(perm(i-1)) * dimr(i-1)
       if (idxr(i) <= 0) then !negative means density index
          ir = ir + stpr(perm(i)) * (-idxr(i) - 1)
          dd(perm(i)) = 1
       else
          dd(perm(i)) = dimr(i)
       end if
    end do
    ! loop over indices in ave and rsp
    ii = 0 !indices from zero to dd-1
    do
       rsp(ir+1) = rsp(ir+1) + prefac * ave(ia+1)
       ! increment indices
       do i = 1, ndim
          ii(i) = ii(i) + 1
          ia = ia + stpa(i)
          ir = ir + stpr(i)
          if (ii(i) /= dd(i)) exit
          ii(i) = 0
          ia = ia - stpa(i) * dd(i)
          ir = ir - stpr(i) * dd(i)
       end do
       if (i == ndim+1) exit
    end do
  end subroutine


end module
