#ifdef PRG_DIRAC
#define GRCONT_NOT_AVAILABLE
#endif

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

! MR: QUICK-FIX USE STATEMENT TO GET SUPPORT FOR DUMMY rsp_cfg TYPE
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

! MR: TEMPORARY ROUTINES RETURNING NON-REDUNDANT VALUES W.R.T TENSOR SYMMETRY
! MR: INCOMPLETE ADAPTATION TO TENSOR SYMMETRY NONREDUNDANCY IN SEVERAL OF THESE ROUTINES
  public rsp_nucpot_tr
  public rsp_ovlave_tr
  public rsp_oneave_tr
  public rsp_twoave_tr
public rsp_xcave_tr_adapt
  public rsp_ovlint_tr
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







! MR: NOT SURE IF WORKING PROPERLY
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

! MR: SIMPLE LOOPS TO ASSIGN VALUES WITH ONLY GEOMETRIC PERTURBATIONS

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
rspfunc_output(h) = rspfunc((i - 1)*ncor**3 + (j - 1)*ncor**2 + (k - 1)*ncor + m)
end do
end do
end do
end do


else if (ngeo > 4) then

write(*,*) 'rsp_nucpot_tr error: No support for ngeo > 4 yet'
call quit('rsp_nucpot_tr error: No support for ngeo > 4 yet')


else

write(*,*) 'rsp_nucpot_tr error: Unknown field setup'
call quit('rsp_nucpot_tr error: Unknown field setup')

end if

else if (ngeo == (nf - 1)) then

if (nf == 2) then

write(*,*) 'rsp_nucpot_tr error: No support for one non-geometrical field with nf = 2 yet'
call quit('rsp_nucpot_tr error: No support for one non-geometrical field  with nf = 2 yet')


else if (nf == 1) then

write(*,*) 'rsp_nucpot_tr error: No support for one non-geometrical field with nf = 1 yet'
call quit('rsp_nucpot_tr error: No support for one non-geometrical field  with nf = 1 yet')



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






! MR: NOT SURE IF WORKING PROPERLY
  !> average f-perturbed overlap integrals with perturbed density D
  !> and energy-weighted density DFD
  subroutine rsp_ovlave_tr(nf, f, c, nc, DFD, nblks, blk_info, & 
                                      blk_sizes, propsize, ave, w, D)
!     use dalton_ifc, only: SHELLS_NUCLEI_displace
    ! Gen1Int interface
!     use gen1int_api
    !> structure containing integral program settings
!     type(rsp_cfg), intent(in)  :: mol
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
    type(matrix),  intent(in)  :: DFD
    !> output average
    complex(8),    intent(out) :: ave(propsize)
    !> field frequencies corresponding to each field
    complex(8),    intent(in), optional  :: w(nf)
    !> density matrix to contract half-differentiated overlap against
    type(matrix),  intent(in), optional  :: D
    !----------------------------------------------
    call interface_1el_ovlave_tr(nf, f, c, nc, DFD, nblks, blk_info, & 
                                      blk_sizes, propsize, ave, w, D)

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



! MR: NOT SURE IF WORKING PROPERLY
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
       A(1) = mat_alloc_like(D1)
       call interface_scf_get_g(D2, A(1)) !Coulomb and exchange
       ave(1) = tr(A(1),D1)
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
       f77_memory(     :n*n)   = reshape(D1%elms_alpha,(/n*n/))
       f77_memory(n*n+1:n*n*2) = reshape(D2%elms_alpha,(/n*n/))
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
       call interest_get_int(D1%nrow, D1%elms_alpha, D2%elms_alpha, 1, 0, size(real_ave), real_ave)
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
       f77_memory(     :n*n)   = reshape(D1%elms_alpha,(/n*n/))
       f77_memory(n*n+1:n*n*2) = reshape(D2%elms_alpha,(/n*n/))
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
       call interest_get_int(D1%nrow, D1%elms_alpha, D2%elms_alpha, 2, 0, size(real_ave), real_ave)
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
       call interest_get_int(D1%nrow, D1%elms_alpha, D2%elms_alpha, 3, 0, size(real_ave), real_ave)
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
    type(matrix)  A(1) !scratch matrices
    type(ctr_arg) arg(1)
    real(8)       r
    integer       h, i, j, k, l, m, n, ncor

  if (any(f== 'EL  ')) then
     ave = 0.0
  else

    if (nf==0) then
       ! contract second density to Fock matrix, then trace with first
       A(1) = mat_alloc_like(D1)
       call interface_scf_get_g(D2, A(1)) !Coulomb and exchange
       ave(1) = tr(A(1),D1)
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
       f77_memory(     :n*n)   = reshape(D1%elms_alpha,(/n*n/))
       f77_memory(n*n+1:n*n*2) = reshape(D2%elms_alpha,(/n*n/))
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
       call interest_get_int(D1%nrow, D1%elms_alpha, D2%elms_alpha, 1, 0, size(real_ave), real_ave)
       do i = 1, size(ave)
          ave(i) = 2.0d0*real_ave(i)
       end do
       deallocate(real_ave)
#endif /* ifdef PRG_DIRAC */

    else if (nf==2 .and. all(f==(/'GEO ','GEO '/))) then
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
       f77_memory(     :n*n)   = reshape(D1%elms_alpha,(/n*n/))
       f77_memory(n*n+1:n*n*2) = reshape(D2%elms_alpha,(/n*n/))
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



!        ave = reshape(tmp(c(1):c(1)+nc(1)-1, &
!                          c(2):c(2)+nc(2)-1,1,1), shape(ave))
       deallocate(tmp)
    else if (nf==3 .and. all(f==(/'GEO ','GEO ','GEO '/))) then
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



       ! extract requested block
!        ave = 2 * reshape(tmp(c(1):c(1)+nc(1)-1, &
!                              c(2):c(2)+nc(2)-1, &
!                              c(3):c(3)+nc(3)-1, 1), shape(ave))
       deallocate(tmp)
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


       ! extract requested block
!        ave = 2 * reshape(tmp(c(1):c(1)+nc(1)-1, &
!                              c(2):c(2)+nc(2)-1, &
!                              c(3):c(3)+nc(3)-1, &
!                              c(4):c(4)+nc(4)-1), shape(ave))
       deallocate(tmp)
    else
       print *, 'rsp_twoave_tr error: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('rsp_twoave_tr error: not implented or in wrong order')
    end if

  end if

  end subroutine


! MR: TEMPORARY ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
! ADAPTED FROM rsp_oneave_tr

  subroutine rsp_xcave_tr_adapt(nr_ao, pert, D_sdf, propsize, ave)

    integer :: i, j, k, m, n, nr_ao, propsize
    type(p_tuple) :: pert, pg, pgg    
! integer,       intent(in)  :: nf
!     !> field labels in std order
!     character(4),  intent(in)  :: f(nf)
!     !> first and number of- components in each field
!     integer,       intent(in)  :: c(nf), nc(nf), propsize
    type(SDF)  :: D_sdf
    type(matrix) :: D
    type(matrix), allocatable, dimension(:) :: Dg
    type(matrix), allocatable, dimension(:,:) :: Dgg
    character :: xcave_pert_label
    !> output average
    complex(8),    intent(out) :: ave(propsize)
    complex(8), allocatable, dimension(:,:,:,:) :: tmp_ave


    allocate(tmp_ave(pert%pdim(1), pert%pdim(1), pert%pdim(1), pert%pdim(1)))
tmp_ave = 0.0

    if (pert%n_perturbations == 0) then
       
       write(*,*) 'rsp_xcave_tr_adapt CALLED WITH NO PERTURBATION'

    else if (pert%n_perturbations == 1) then

       if (all(pert%plab==(/'GEO '/))) then

!           allocate(xcave_pert_label(1))
          xcave_pert_label = 'g'

          ! MR: ASSUME CLOSED SHELL
          call mat_init(D, nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)

          call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), D)

          call rsp_xcave('g', tmp_ave(:, 1, 1, 1), D=D)

          do i = 1, pert%pdim(1)

             ave(i) = ave(i) + tmp_ave(i, 1, 1, 1)

          end do

!           deallocate(xcave_pert_label)

       else

          write(*,*) 'WARNING: UNSUPPORTED CONTRIBUTION: NO CONTRIBUTION WILL BE MADE'

       end if

    else if (pert%n_perturbations == 2) then

       if (all(pert%plab==(/'GEO ','GEO '/))) then

!           allocate(xcave_pert_label(2))
          xcave_pert_label = 'gg'
          allocate(Dg(pert%pdim(1)))

          ! MR: ASSUME CLOSED SHELL
          call mat_init(D, nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)

          call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), D)

          pg = p_tuple_getone(pert, 1)
          
          do i = 1, pert%pdim(1)

             ! MR: ASSUME CLOSED SHELL
             call mat_init(Dg(i), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
             call sdf_getdata_s(D_sdf, pg, (/i/), Dg(i))

          end do

          call rsp_xcave('gg', tmp_ave(:, :, 1, 1), D=D, Dg=Dg)

          do i = 1, pert%pdim(1)
             do j = 1, i

                n = get_triang_blks_offset(1, 2, (/1, 2, pert%pdim(1)/), &
                                           (/propsize/), (/i, j/))

                ave(n) = ave(n) + tmp_ave(i, j, 1, 1)

             end do
          end do

          deallocate(Dg)
!           deallocate(xcave_pert_label)

       else

          write(*,*) 'WARNING: UNSUPPORTED CONTRIBUTION: NO CONTRIBUTION WILL BE MADE'

       end if

    else if (pert%n_perturbations == 3) then

       if (all(pert%plab==(/'GEO ','GEO ','GEO '/))) then

!           allocate(xcave_pert_label(3))
          xcave_pert_label = 'ggg'
          allocate(Dg(pert%pdim(1)))

          ! MR: ASSUME CLOSED SHELL
          call mat_init(D, nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
          call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), D)

! write(*,*) 'elms of D', d%elms_alpha

          pg = p_tuple_getone(pert, 1)
! write(*,*) 'pg plab:', pg%plab
          do i = 1, pert%pdim(1)
! write(*,*) 'i', i
             ! MR: ASSUME CLOSED SHELL
             call mat_init(Dg(i), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
             call sdf_getdata_s(D_sdf, pg, (/i/), Dg(i))

! write(*,*) 'Dg(i)', Dg(i)%elms_alpha
          end do

          call rsp_xcave('ggg', tmp_ave(:, :, :, 1), D=D, Dg=Dg)

! write(*,*) 'ave size', size(ave)
! write(*,*) 'full ave', tmp_ave(:,:,:,1)
          do i = 1, pert%pdim(1)
             do j = 1, i
                do k = 1, j


                   n = get_triang_blks_offset(1, 3, (/1, 3, pert%pdim(1)/), &
                                              (/propsize/), (/i, j, k/))

write(*,*) 'n is', n

                   ave(n) = ave(n) + tmp_ave(i, j, k, 1)

                end do
             end do
          end do

          deallocate(Dg)
!           deallocate(xcave_pert_label)

       else

          write(*,*) 'WARNING: UNSUPPORTED CONTRIBUTION: NO CONTRIBUTION WILL BE MADE'

       end if

    else if (pert%n_perturbations == 4) then

       if (all(pert%plab==(/'GEO ','GEO ','GEO ','GEO '/))) then

!           allocate(xcave_pert_label(4))
          xcave_pert_label = 'gggg'
          allocate(Dg(pert%pdim(1)))
          allocate(Dgg(pert%pdim(1), pert%pdim(1)))

          ! MR: ASSUME CLOSED SHELL
          call mat_init(D, nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
          call sdf_getdata_s(D_sdf, get_emptypert(), (/1/), D)

          pg = p_tuple_getone(pert, 1)
          call p_tuple_p1_cloneto_p2(pert, pgg)
          pgg = p_tuple_remove_first(pgg)
          pgg = p_tuple_remove_first(pgg)


write(*,*) 'pgg plab:', pgg%plab

          do i = 1, pert%pdim(1)

             ! MR: ASSUME CLOSED SHELL
             call mat_init(Dg(i), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
             call sdf_getdata_s(D_sdf, pg, (/i/), Dg(i))

             do j = 1, pert%pdim(1)


                ! MR: ASSUME CLOSED SHELL
                call mat_init(Dgg(i,j), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
                call sdf_getdata_s(D_sdf, pgg, (/i,j/), Dgg(i,j))

             end do
          end do

          call rsp_xcave('gggg', tmp_ave, D=D, Dg=Dg, Dgg=Dgg)


          do i = 1, pert%pdim(1)
             do j = 1, i
                do k = 1, j
                   do m = 1, k

                       n = get_triang_blks_offset(1, 4, (/1, 4, pert%pdim(1)/), &
                                                  (/propsize/), (/i, j, k, m/))


                      ave(n) = ave(n) + tmp_ave(i, j, k, m)


                   end do
                end do
             end do
          end do

          deallocate(Dgg)
          deallocate(Dg)
!           deallocate(xcave_pert_label)

       else

          write(*,*) 'WARNING: UNSUPPORTED CONTRIBUTION: NO CONTRIBUTION WILL BE MADE'

       end if

    else

       write(*,*) 'WARNING: UNSUPPORTED NUMBER OF FIELDS: NO CONTRIBUTION WILL BE MADE'

    end if

    deallocate(tmp_ave)

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




! MR: TEMPORARY ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
  !> Compute differentiated overlap matrices, and optionally
  !> add half-differentiated overlap contribution to Fock matrices
  subroutine rsp_ovlint_tr(nr_ao, nf, f, c, nc, nblks, blk_info, & 
                                      blk_sizes, propsize, ovl, w, fock)
    ! Gen1Int interface in Dalton
!     use gen1int_api
    !> structure containing integral program settings
!     type(rsp_cfg), intent(in)    :: mol
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
    !> frequencies of each field
    complex(8),    intent(in),    optional :: w(nf)
    !> Fock matrices to which the half-differentiated overlap
    !> contribution is ADDED
    type(matrix),  intent(inout), optional :: fock(propsize)
    !------------------------------------------------
    integer      i, nr_ao

    call interface_1el_ovlint_tr(nr_ao, nf, f, c, nc, nblks, blk_info, & 
                                      blk_sizes, propsize, ovl, w, fock)

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





! MR: TEMPORARY ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN


  subroutine rsp_oneint_tr(nr_ao, nf, f, c, nc, nblks, blk_info, & 
                                      blk_sizes, propsize, oneint)
    ! Gen1Int interface in Dalton
!     use gen1int_api
    !> structure containing integral program settings
!     type(rsp_cfg), intent(in)    :: mol
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
     call mat_init(A, nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
     do i = 1, product(nc)
        if (iszero(fock(i))) then
           call mat_ensure_alloc(fock(i))
           fock(i)%elms_alpha = fock(i)%elms_alpha + A%elms_alpha
        else
           fock(i)%elms_alpha = fock(i)%elms_alpha + A%elms_alpha
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
                         1, (c(1)+i+2)/3, .false., .true., dens%elms_alpha, 1)
          end if
          j = 1 + mod(c(1)+i-1,3) !x y z = 1 2 3
          if (iszero(fock(1+i))) then
             fock(1+i)%elms_alpha(:, :, 1) = reshape(f77_memory(n*n*(j-1)+1:n*n*j),(/n,n/))
          else
             fock(1+i)%elms_alpha(:, :, 1) = fock(1+i)%elms_alpha(:, :, 1) &
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
          call interest_get_int(dens%nrow, dens%elms_alpha, fock(i)%elms_alpha, 1, i, 0, dummy)
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
                fock(ij)%elms_alpha = 0 !ajt FIXME use mat_axpy
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



! MR: TEMPORARY ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN

  !> Contract 2-electron integrals perturbed by fields 'f' with density
  !> matrix 'dens', and add to Fock matrices 'fock' Average 2-electron integrals perturbed by fields f over the
  !> product of density matrces D1 and D2
  subroutine rsp_twoint_tr(nr_ao, nf, f, c, nc, dens, propsize, fock)
    use eri_contractions, only: ctr_arg
    use eri_basis_loops,  only: unopt_geodiff_loop
    use interface_interest
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
    !--------------------------------------------------
    real(8), pointer :: null_ptr(:) !because null() isn't f90
    integer       i, j, n, ij, ncor, nr_ao
    type(ctr_arg) arg(1)
    type(matrix)  A !scratch
    real(8)          :: dummy(1)

! MR: DUMMY PLACEHOLDER MOL DECLARATION

! integer :: mol

  if (any(f=='EL  ')) then
     call mat_init(A, nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)
     do i = 1, propsize
        if (iszero(fock(i))) then
           call mat_ensure_alloc(fock(i))
           fock(i)%elms_alpha = fock(i)%elms_alpha + A%elms_alpha
        else
           fock(i)%elms_alpha = fock(i)%elms_alpha + A%elms_alpha
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
                         1, (c(1)+i+2)/3, .false., .true., dens%elms_alpha, 1)
          end if
          j = 1 + mod(c(1)+i-1,3) !x y z = 1 2 3
          if (iszero(fock(1+i))) then
             fock(1+i)%elms_alpha(:, :, 1) = reshape(f77_memory(n*n*(j-1)+1:n*n*j),(/n,n/))
          else
             fock(1+i)%elms_alpha(:, :, 1) = fock(1+i)%elms_alpha(:, :, 1) &
                            + reshape(f77_memory(n*n*(j-1)+1:n*n*j),(/n,n/))
          end if
#endif
#endif /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
          call interest_mpi_wake_up()
          call interest_get_int(dens%nrow, dens%elms_alpha, fock(i+1)%elms_alpha, 1, i+1, 0, dummy)
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
                fock(ij)%elms_alpha = 0 !ajt FIXME use mat_axpy
             end if
             arg(1) = ctr_arg(2, c(1)+i + ncor * (c(2)+j-1), &
                              ncor, dens, fock(ij), null_ptr)
             call unopt_geodiff_loop(basis_large, &
                                     basis_small, &
                                     arg)
          end do
       end do
    else
       print *, 'error in rsp_twoint_tr: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('error in rsp_twoint_tr: not implented or in wrong order')
    end if

  end if

  end subroutine






! MR: TEMPORARY ROUTINE FOR TENSOR SYMMETRY NONREDUNDANT DATA RETURN
! ADAPTED FROM rsp_oneint_tr


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
    type(matrix), allocatable, dimension(:,:) :: tmp_xcint
    type(matrix), dimension(:) :: D
    !--------------------------------------------------


    ! MR: ONLY GEOMETRICAL PERTURBATIONS SUPPORTED SO FAR


    if (nf == 0) then

      call rsp_xcint(D, F=xcint(1))

    else if (nf == 1) then

       if (all(f==(/'GEO '/))) then

          call rsp_xcint(D, Fg=xcint)

       else

          write(*,*) 'WARNING: UNSUPPORTED CONTRIBUTION: NO CONTRIBUTION WILL BE MADE'

       end if

    else if (nf == 2) then

       if (all(f==(/'GEO ','GEO '/))) then

          allocate(tmp_xcint(nc(1), nc(1)))

          do i = 1, nc(1)
             do j = 1, nc(1)

                ! MR: ASSUME CLOSED SHELL
                call mat_init(tmp_xcint(i,j), nrow=nr_ao, ncol=nr_ao, closed_shell=.true.)

             end do
          end do

          call rsp_xcint(D, Fgg=tmp_xcint)


          do i = 1, nc(1)
             do j = 1, i

                k = get_triang_blks_offset(1, 2, (/1, 2, nc(1)/), &
                                           (/propsize/), (/i, j/))

                xcint(k) = tmp_xcint(i,j)
                tmp_xcint(i, j) = 0

             end do
          end do

          deallocate(tmp_xcint)

       else

          write(*,*) 'WARNING: UNSUPPORTED CONTRIBUTION: NO CONTRIBUTION WILL BE MADE'

       end if

    else

       write(*,*) 'WARNING: UNSUPPORTED NUMBER OF FIELDS: NO CONTRIBUTION WILL BE MADE'

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
