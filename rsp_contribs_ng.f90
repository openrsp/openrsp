! Copyright 2009-2011 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License.

!> @file
!> Contains module prop_contribs

!> This module contains routines for calculating contributions
!> to molecular properties (1st order, linear response, etc.),
!> and perturbed Fock matrices. 
module rsp_contribs_ng

  use matrix_defop_ng

#ifdef BUILD_AORSP
  use dalton_ifc
  use basis_set, only: cgto
#endif

  implicit none
  public rsp_oneave
  public rsp_twoave
  public rsp_oneint
  public rsp_twoint
  public rsp_field
  public rsp_field_bas
  public rsp_field_dim
  public rsp_field_anti
  public rsp_cfg


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


  !> Molecule configuration type, abstracting data and settings
  !> to be passed to solver and integral routines.
  !> Should eventually be moved to separate program-specific
  !> interface modules
  type rsp_cfg
     !> prototype zero matrix, of same shape and multiplicity as
     !> then density matrix
     type(matrix) :: zeromat
     !> number of atoms, determining ncomp of GEO and NMAG
     integer      :: natom
     !> unit number for printing output
     integer      :: lupri
     !> whether energy/response is Kohn-Sham DFT
     logical      :: isdft
     !> nuclear charges (natom)
     real(8),    pointer :: charge(:)
     !> nuclear coordinates (3,natom)
     real(8),    pointer :: coord(:,:)
     !> basis
     type(cgto), pointer :: basis(:)
  end type


  !> private struct to collect properties of perturbing "fields"
  type fld_info
     !> four-letter abbreviation
     character(4)  :: code
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


  ! to compactify the table
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
  type(fld_info) :: field_list(13) = &                         !nc an ba ln qu
     (/fld_info('EXCI', 'Generalized "excitation" field'      , 1, F, F, T, T), &
       fld_info('FREQ', 'Generalized "freqency" field'        , 1, F, F, T, T), &
       fld_info('EL'  , 'Electric field'                      , 3, F, F, T, F), &
       fld_info('VEL' , 'Velocity'                            , 3, T, F, T, F), &
       fld_info('MAGO', 'Magnetic field w/o. London orbitals' , 3, T, F, F, T), &
       fld_info('MAG' , 'Magnetic field with London orbitals' , 3, T, T, F, F), &
       fld_info('ELGR', 'Electric field gradient'             , 6, F, F, T, F), &
       fld_info('VIBM', 'Displacement along vibrational modes',-1, F, T, F, F), &
       fld_info('GEO' , 'Nuclear coordinates'                 ,-1, F, T, F, F), & !-1=mol-dep
       fld_info('NUCM', 'Nuclear magnetic moment'             ,-1, F, T, F, T), & !-1=mol-dep
       fld_info('AOCC', 'AO contraction coefficients'         ,-1, F, T, F, F), & !-1=mol-dep
       fld_info('AOEX', 'AO exponents'                        ,-1, F, T, F, F)/)  !-1=mol-dep

  character(1), parameter :: xyz(3) = (/'X','Y','Z'/)

  private

contains


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



  !> nuclear repulsion and nuclei--field interaction
  subroutine rsp_nucpot(mol, nf, f, w, c, nc, pot)
    !> structure containing integral program settings
    type(rsp_cfg),   intent(in) :: mol
    !> number of fields
    integer,         intent(in) :: nf
    !> field lables in std order
    character(4),    intent(in) :: f(nf)
    !> field frequencies corresponding to each field
    complex(8),      intent(in) :: w(nf)
    !> first and number of- components in each field
    integer,         intent(in) :: c(nf), nc(nf)
    !> output average
    complex(8),     intent(out) :: pot(product(nc))
    !----------------------------------------------
    real(8) :: tmp(3*mol%natom)
    integer :: i
    if (nf==0) then
       call quit('rsp_nucpot error: unperturbed (nf=0) nuc.rep. not implemented')
    else if (nf==1 .and. f(1)=='GEO') then
       call GRADNN_ifc(mol%natom, tmp)
       pot(:nc(1)) = tmp(c(1):c(1)+nc(1)-1)
    else
       print *,'rsp_nucpot error: not implented or in wrong order - ', &
               (' ' // f(i), i=1,nf)
       call quit('rsp_nucpot error: not implented or in wrong order')
    end if
  end subroutine



  !> average f-perturbed overlap integrals with perturbed density D
  !> and energy-weighted density DFD
  subroutine rsp_ovlave(mol, nf, f, w, c, nc, D, DFD, ave)
    !> structure containing integral program settings
    type(rsp_cfg),   intent(in) :: mol
    !> number of fields
    integer,         intent(in) :: nf
    !> field lables in std order
    character(4),    intent(in) :: f(nf)
    !> field frequencies corresponding to each field
    complex(8),      intent(in) :: w(nf)
    !> first and number of- components in each field
    integer,         intent(in) :: c(nf), nc(nf)
    !> density and energy-weighted density matrix
    type(matrix),    intent(in) :: D, DFD
    !> output average
    complex(8),     intent(out) :: ave(product(nc))
    !----------------------------------------------
    type(matrix) :: A(2)
    integer      :: i
    if (nf==0) then
       call quit('rsp_ovlave error: unperturbed (nf=0) overlap not implemented')
    else if (nf==1 .and. f(1)=='GEO') then
       ! allocate matrices for integrals
       if (w(1)/=0) call mat_init(A(1), mol%zeromat)
       call mat_init(A(2), mol%zeromat)
       ! loop over nuclear coordinates
       do i = 0, nc(1)-1
          ! (half-) perturbed overlap -i/2 Tg into A(1), Sg in A(2)
          if (w(1)==0) then !w=0 means no -i/2 Tg contribution
             call di_read_operator_int('1DOVL' // prefix_zeros(c(1)+i,3), A(2))
             A(2) = -A(2) !1DOVL is really -dS/dg
          else
             call di_read_operator_int('SQHDR' // prefix_zeros(c(1)+i,3), A(2))
             A(2) = -A(2) !SQHDR is really -dS>/dg
             A(1) = A(1) - w(1)/2 * A(2)
             A(1) = A(1) + w(1)/2 * dag(A(2))
             A(2) = A(2) + dag(A(2))
          end if
          ave(1+i) = tr(A(1),D) - tr(A(2),DFD)
       end do
       A(1:2) = 0 !deallocate
    else
       print *, 'rsp_ovlave error: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('rsp_ovlave error: not implented or in wrong order')
    end if
  end subroutine



  !> Average 1-electron integrals perturbed by fields f
  !> with the (perturbed) density matrix D
  subroutine rsp_oneave(mol, nf, f, c, nc, D, ave)
    !> structure containing integral program settings
    type(rsp_cfg),   intent(in) :: mol
    !> number of fields
    integer,         intent(in) :: nf
    !> field lables in std order
    character(4),    intent(in) :: f(nf)
    !> first and number of- components in each field
    integer,         intent(in) :: c(nf), nc(nf)
    !> density matrix to average over
    type(matrix),    intent(in) :: D
    !> output average
    complex(8),     intent(out) :: ave(product(nc))
    !----------------------------------------------
    type(matrix) :: A(1)
    integer      :: i
    if (nf==0) then
       call quit('rsp_oneave error: unperturbed (nf=0) 1el.int. not implemented')
    else if (nf==1 .and. f(1)=='GEO') then
       call mat_init(A(1), mol%zeromat)
       do i = 0, nc(1)-1
          ! perturbed one-electron Hamiltonian integrals
          call di_read_operator_int('1DHAM' // prefix_zeros(c(1)+i,3), A(1))
          ave(1+i) = tr(A(1),D)
       end do
       A(1) = 0
    else
       print *, 'rsp_oneave error: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('rsp_oneave error: not implented or in wrong order')
    end if
  end subroutine



  !> Average 2-electron integrals perturbed by fields f over the
  !> product of density matrces D1 and D2
  subroutine rsp_twoave(mol, nf, f, c, nc, D1, D2, ave)
    !> structure containing integral program settings
    type(rsp_cfg),   intent(in) :: mol
    !> number of fields
    integer,         intent(in) :: nf
    !> field lables in std order
    character(4),    intent(in) :: f(nf)
    !> first and number of- components in each field
    integer,         intent(in) :: c(nf), nc(nf)
    !> density matrix to average over
    type(matrix),    intent(in) :: D1, D2
    !> output average
    complex(8),     intent(out) :: ave(product(nc))
    !----------------------------------------------
    real(8), allocatable :: work(:)
    real(8)      :: temp(3*mol%natoms) !scratch
    type(matrix) :: A(1) !scratch matrices
    integer      :: i, n
    if (nf==0) then
       ! contract second density to Fock matrix, then trace with first
       call mat_init(A(1), mol%zeromat)
       call di_get_gmat(D2, A(1)) !Coulomb and exchange
       ave(1) = tr(A(1),D1)
    else if (nf==1 .and. f(1)=='GEO') then
       n = mol%nbas
       allocate(work(5*n*n + 20*n + 1000))
       work(:n*n)        = reshape(D1%elms,(/n*n/))
       work(n*n+1:n*n*2) = reshape(D2%elms,(/n*n/))
       call GRCONT(work(n*n*2+1:), size(work)-n*n*2, temp, size(temp), &
                   .true., .false., 1, 0, .true., .false., work(:n*n*2), 2)
       ave(:nc(1)) = temp(c(1):c(1)+nc(1)-1)
       deallocate(work)
    else
       print *, 'rsp_twoave error: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('rsp_twoave error: not implented or in wrong order')
    end if
  end subroutine



  !> Exchange-correlation perturbed by fields f, averaged over densities D
  subroutine rsp_ksmave(mol, nf, f, c, nc, nd, D, ave)
    !> structure containing integral program settings
    type(rsp_cfg),   intent(in) :: mol
    !> number of fields
    integer,         intent(in) :: nf
    !> field lables in std order
    character(4),    intent(in) :: f(nf)
    !> first and number of- components in each field
    integer,         intent(in) :: c(nf), nc(nf)
    !> number of density matrices
    integer,         intent(in) :: nd
    !> density matrices to average over. D(1) must be unperturbed
    type(matrix),    intent(in) :: D(nd)
    !> output average
    complex(8),     intent(out) :: ave(product(nc))
    !----------------------------------------------
    ave = 0
  end subroutine



  function idx(f)
    character(4) :: f
    integer      :: idx
    do idx = 1, size(field_list)
        if (field_list(idx)%label == f) return
    end do
    call quit('Field not found: ' // f)
  end function


  function rsp_field_anti(f)
    character(4), intent(in) :: f(:)
    logical :: rsp_field_anti(size(f))
    integer :: i
    rsp_field_anti = (/(field_list(idx(f(i)))%anti, i=1,size(f))/)
  end function


  !> shape (dimensions) of property for fields f(:)
  function rsp_field_dim(mol, f)
    !> structure containing the integral program settings
    type(rsp_cfg), intent(in) :: mol
    !> field lables
    character(4),  intent(in) :: f(:)
    integer :: rsp_field_dim(size(f)), i
    rsp_field_dim = (/(field_list(idx(f(i)))%ncomp, i=1,size(f))/)
    ! loop through mol-dependent
    do i = 1, size(f)
       if (rsp_field_dim(i) /= -1) then
          ! cycle
       else if (f(i)=='GEO') then
          rsp_field_dim(i) = 3 * mol%natoms
       else
          call quit('rsp_field_dim error: Number of comp. unknown for ' // f(i))
       end if
    end do
  end function


  function rsp_field_bas(f)
    character(4), intent(in) :: f(:)
    logical :: rsp_field_bas(size(f))
    integer :: i
    rsp_field_bas = (/(field_list(idx(f(i)))%bas, i=1,size(f))/)
  end function


  !> Same as set_dsofso in fock-eval.f90, but pertrubed
  !> D and DFD is input instead of unperturbed D and F
  subroutine save_D_and_DFD_for_ABACUS(mol, anti, D, DFD)
    !> whether the integrals D and DFD are to be contracted with are
    !> symmetric or anti-symmetric
    logical,      intent(in) :: anti 
    type(matrix), intent(in) :: D, DFD
       ! perturbed (or un-) density and
       ! 'generalized Fock matrix' (energy-weighted density matrix)
    real(8) :: Dtri(mol%nbas*(mol%nbas+1)/2), &
             DFDtri(mol%nbas*(mol%nbas+1)/2)
    if (iszero(D)) then
       Dtri = 0
    else
       if (.not.anti) call DGEFSP(D%nrow, D%elms, Dtri)
       if (     anti) call DGETAP(D%nrow, D%elms, Dtri)
    end if
    if (iszero(DFD)) then
       DFDtri = 0
    else
       if (.not.anti) call DGEFSP(D%nrow, D%elms, DFDtri)
       if (     anti) call DGETAP(D%nrow, D%elms, DFDtri)
    end if
    ! write fo files
    call write_dsofso(Dtri,DFDtri)
  end subroutine


end module
