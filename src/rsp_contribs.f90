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
  use dalton_ifc
  use basis_set,  only: cgto

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
     !> basis
     type(cgto), pointer :: basis(:)
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
    integer      ncor, ngeo, last_ncomp
    character(4) last_field
    ncor = 3 * get_nr_atoms()
    ngeo = 0
    last_field = 'NONE'
    last_ncomp = 1
    !ajt FIXME validate comp/ncomp ranges
    !ajt FIXME determine sorting
    if (any(fields(:)%label == 'EL  ')) then
       if (size(fields) > 2) then
          rspfunc = 0.0
       end if
    else
       ! count the number of GEO
       if (all(fields(:size(fields)-1)%label == 'GEO ')) then
          ngeo = size(fields) - 1
          if (fields(size(fields))%label == 'GEO ') then
             ngeo = ngeo + 1
          else
             last_field = fields(size(fields))%label
             last_ncomp = 3 !ajt FIXME
          end if
       else
          call quit('rsp_nucpot error: failed to parse fields')
       end if
       call inner
    end if
  contains
    subroutine inner
      real(8) tmp(ncor**ngeo * last_ncomp)
      call nuclear_potential(ngeo, ncor, last_field, last_ncomp, tmp)
      ! add selected rectangle to rspfunc
      !ajt FIXME reordering and range selection
      if (ncor**ngeo * last_ncomp /= size(rspfunc)) &
         call quit('rsp_nucpot error: only full ranges implemented')
      rspfunc = rspfunc + tmp
    end subroutine
  end subroutine



  !> average f-perturbed overlap integrals with perturbed density D
  !> and energy-weighted density DFD
  subroutine rsp_ovlave(mol, nf, f, c, nc, DFD, ave, w, D)
    !> structure containing integral program settings
    type(rsp_cfg), intent(in)  :: mol
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



  !> Average 1-electron integrals perturbed by fields f
  !> with the (perturbed) density matrix D

! radovan: there is code repetition in ave and int setup

  subroutine rsp_oneave(mol, nf, f, c, nc, D, ave)
    !> structure containing integral program settings
    type(rsp_cfg), intent(in)  :: mol
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


  function rank_one_pointer(siz, arr) result(ptr)
     integer,         intent(in) :: siz
     real(8), target, intent(in) :: arr(siz)
     real(8), pointer            :: ptr(:)
     ptr => arr
  end function


  !> Average 2-electron integrals perturbed by fields f over the
  !> product of density matrces D1 and D2
  subroutine rsp_twoave(mol, nf, f, c, nc, D1, D2, ave)
    use eri_contractions, only: ctr_arg
    use eri_basis_loops,  only: unopt_geodiff_loop
    !> structure containing integral program settings
    type(rsp_cfg),        intent(in)  :: mol
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
       A(1) = mol%zeromat
       call mat_ensure_alloc(A(1))
       call di_get_gmat(D2, A(1)) !Coulomb and exchange
       ave(1) = tr(A(1),D1)
    else if (nf==1 .and. f(1)=='GEO ') then
       ncor = 3 * get_nr_atoms()
       allocate(tmp(ncor,1,1,1))
       n = mol%zeromat%nrow
       f77_memory(     :n*n)   = reshape(D1%elms,(/n*n/))
       f77_memory(n*n+1:n*n*2) = reshape(D2%elms,(/n*n/))
       call GRCONT(f77_memory(n*n*2+1:), size(f77_memory)-n*n*2, &
                   tmp(:,1,1,1), ncor, .true., .false., &
                   1, 0, .true., .false., f77_memory(:n*n*2), 2)
       ave(:nc(1)) = tmp(c(1):c(1)+nc(1)-1,1,1,1)
       deallocate(tmp)
    else if (nf==2 .and. all(f==(/'GEO ','GEO '/))) then
       ncor = 3 * get_nr_atoms()
       allocate(tmp(ncor,ncor,1,1))
       n = mol%zeromat%nrow
       f77_memory(     :n*n)   = reshape(D1%elms,(/n*n/))
       f77_memory(n*n+1:n*n*2) = reshape(D2%elms,(/n*n/))
       call GRCONT(f77_memory(n*n*2+1:), size(f77_memory)-n*n*2, &
                   tmp(:,:,1,1), ncor**2, .true., .false., &
                   2, 0, .true., .false., f77_memory(:n*n*2), 2)
       ave = reshape(tmp(c(1):c(1)+nc(1)-1, &
                         c(2):c(2)+nc(2)-1,1,1), shape(ave))
       deallocate(tmp)
    else if (nf==3 .and. all(f==(/'GEO ','GEO ','GEO '/))) then
       ! contract FULL cubic in tmp, unsymmetrized divided by six
       ncor = 3 * get_nr_atoms()
       allocate(tmp(ncor,ncor,ncor,1))
       arg(1) = ctr_arg(3, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor**3, tmp(:,:,:,1)))
       call unopt_geodiff_loop(mol%basis, arg)
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
    else if (nf==4 .and. all(f==(/'GEO ','GEO ','GEO ','GEO '/))) then
       ncor = 3 * get_nr_atoms()
       allocate(tmp(ncor,ncor,ncor,ncor))
       ! contract FULL quartic in tmp, unsymmetrized divided by 24
       arg(1) = ctr_arg(4, -huge(1), ncor, D1, D2, &
                        rank_one_pointer(ncor**4, tmp))
       call unopt_geodiff_loop(mol%basis, arg)
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



  !> Compute differentiated overlap matrices, and optionally
  !> add half-differentiated overlap contribution to Fock matrices
  subroutine rsp_ovlint(mol, nf, f, c, nc, ovl, w, fock)
    !> structure containing integral program settings
    type(rsp_cfg), intent(in)    :: mol
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
    call interface_1el_ovlint(mol%zeromat, nf, f, c, nc, ovl, w, fock)
  end subroutine



  subroutine rsp_oneint(mol, nf, f, c, nc, oneint)
    !> structure containing integral program settings
    type(rsp_cfg), intent(in)    :: mol
    !> number of fields
    integer,       intent(in)    :: nf
    !> field labels in std order
    character(4),  intent(in)    :: f(nf)
    !> first and number of- components in each field
    integer,       intent(in)    :: c(nf), nc(nf)
    !> output perturbed integrals
    type(matrix),  intent(inout) :: oneint(product(nc))
    call interface_1el_oneint(mol%zeromat, nf, f, c, nc, oneint)
  end subroutine



  !> Contract 2-electron integrals perturbed by fields 'f' with density
  !> matrix 'dens', and add to Fock matrices 'fock' Average 2-electron integrals perturbed by fields f over the
  !> product of density matrces D1 and D2
  subroutine rsp_twoint(mol, nf, f, c, nc, dens, fock)
    ! work array to be passed to GRCONT
    use eri_contractions, only: ctr_arg
    use eri_basis_loops,  only: unopt_geodiff_loop
    !> structure containing integral program settings
    type(rsp_cfg),        intent(in)    :: mol
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
  if (any(f=='EL  ')) then
        A = mol%zeromat
        call mat_ensure_alloc(A)
        A%elms = 0.0
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
       call di_get_gmat(dens, A)
       fock(1) = fock(1) + A
       A = 0
    else if (nf==1 .and. f(1)=='GEO ') then
       n = mol%zeromat%nrow
       do i = 0, nc(1)-1
          ! if first or an x-coord, call GRCONT
          if (i==0 .or. mod(c(1)+i,3) == 1) &
             call GRCONT(f77_memory(n*n*3+1:), size(f77_memory)-n*n*3, &
                         f77_memory(:n*n*3), n*n*3, .true., .false., &
                         1, (c(1)+i+2)/3, .false., .true., dens%elms, 1)
          j = 1 + mod(c(1)+i-1,3) !x y z = 1 2 3
          if (iszero(fock(1+i))) then
             call mat_ensure_alloc(fock(1+i))
             fock(1+i)%elms = reshape(f77_memory(n*n*(j-1)+1:n*n*j),(/n,n/))
          else
             fock(1+i)%elms = fock(1+i)%elms &
                            + reshape(f77_memory(n*n*(j-1)+1:n*n*j),(/n,n/))
          end if
       end do
    else if (nf==2 .and. all(f==(/'GEO ','GEO '/))) then
       ncor = 3 * get_nr_atoms()
       nullify(null_ptr) !because null() isn't f90
       do j = 0, nc(2)-1
          do i = 0, nc(1)-1
             ij = 1 + i + nc(1)*j
             if (iszero(fock(ij))) then
                call mat_ensure_alloc(fock(ij))
                fock(ij)%elms = 0 !ajt FIXME use mat_axpy
             end if
             arg(1) = ctr_arg(2, c(1)+i + ncor * (c(2)+j-1), &
                              ncor, dens, fock(ij), null_ptr)
             call unopt_geodiff_loop(mol%basis, arg)
          end do
       end do
    else
       print *, 'error in rsp_oneave: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('error in rsp_oneave: not implented or in wrong order')
    end if

  end if

  end subroutine





  function idx(f)
    character(4) :: f
    integer      :: idx
    do idx = 1, size(all_known_fields)
        if (all_known_fields(idx)%label == f) return
    end do
    call quit('Field not found: ' // f)
  end function


  function rsp_field_anti(f)
    character(4), intent(in) :: f(:)
    logical :: rsp_field_anti(size(f))
    integer :: i
    rsp_field_anti = (/(all_known_fields(idx(f(i)))%anti, i=1,size(f))/)
  end function


  !> shape (dimensions) of property for fields f(:)
  function rsp_field_dim(f)
    character(4),  intent(in) :: f(:)
    integer :: rsp_field_dim(size(f)), i
    rsp_field_dim = (/(all_known_fields(idx(f(i)))%ncomp, i=1,size(f))/)
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
    rsp_field_bas = (/(all_known_fields(idx(f(i)))%bas, i=1,size(f))/)
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



  !> Same as set_dsofso in fock-eval.f90, but pertrubed
  !> Call ONEDRV in ABACUS



  !> Find the reordering of type(rsp_field)s f(:) that puts them in "canonical"
  !> order, sorted by:
  !>    1) decreasing %label's index in all_known_fields (GEO before EL, etc.)
  !>    2) decreasing number of components %ncomp, (=1 for equations)
  !>    3) increasing starting component index %comp,
  !>    4) decreasing absolute value of real part of %freq (- before +)
  !>    5) decreasing absolute value of imaginary part of %freq (- before +)
  !> Canonical order should be used in response equation solution caching,
  !> and is due to 1) also the order delivered by the integral backends.
  function rsp_field_ordering(f) result(o)
    type(rsp_field), intent(in) :: f(:)
    integer                     :: o(size(f))
    integer i, j, k
    o = (/(i, i=1, size(f))/)
    do i = 1, size(f)
       !find perturbation after i, with highest idx (secondly highest %ncomp)
       j = i
       do k = j+1, size(f)
          if (idx(f(o(k))%label) <  idx(f(o(j))%label)) cycle
          if (idx(f(o(k))%label) == idx(f(o(j))%label)) then
             if (f(o(k))%ncomp <  f(o(j))%ncomp) cycle
             if (f(o(k))%ncomp == f(o(j))%ncomp) then
                if (f(o(k))%comp >  f(o(j))%comp) cycle
                if (f(o(k))%comp == f(o(j))%comp) then
                   if (   abs(real(f(o(k))%freq)) &
                       <  abs(real(f(o(j))%freq))) cycle
                   if (   abs(real(f(o(k))%freq)) &
                       == abs(real(f(o(j))%freq))) then
                      if (   abs(aimag(f(o(k))%freq)) &
                          <= abs(aimag(f(o(j))%freq))) cycle
                   end if
                end if
             end if
          end if
          j = k  !new minimum
       end do
       !swap entries i and j
       k    = o(i)
       o(i) = o(j)
       o(j) = k
    end do
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


  !> (derivatives of) nuclear repulsion and nuclei--field interaction
  subroutine nuclear_potential(derv, ncor, field, ncomp, nucpot)
    !> order of geometry derivative requested
    integer,      intent(in)  :: derv
    !> number of geometry coordinates
    integer,      intent(in)  :: ncor
    !> one external potential label, or 'NONE'
    character(4), intent(in)  :: field
    !> number of components of the field (must be all, ie. 1/3/6)
    integer,      intent(in)  :: ncomp
    !> output tensor
    real(8),      intent(out) :: nucpot(ncor**derv * ncomp)
    !-------------------------------------------------------
    integer i, na
    na = get_nr_atoms()
    if (ncor /= 3*na) &
       call quit('rsp_backend nuclear_potential error: expected ncor = 3*natom')
    if ((field == 'NONE' .and. ncomp /= 1) .or. &
        (field == 'EL  ' .and. ncomp /= 3)) &
       call quit('rsp_backend nuclear_potential error: expected ncomp = 1 or 3')
    ! switch
    if (derv == 0 .and. field == 'NONE') then
       call quit('rsp_ error: unperturbed nuc.rep. not implemented')
    else if (derv == 0 .and. field == 'EL  ') then
       call DIPNUC_ifc(nucpot)
       nucpot = -nucpot !dip = -dE/dF
    else if (derv == 1 .and. field == 'NONE') then
       call GRADNN_ifc(na, nucpot)
    else if (derv == 2 .and. field == 'NONE') then
       call HESSNN_ifc(na, nucpot)
    else if (derv == 3 .and. field == 'NONE') then
       call cubicff_nuc(na, nucpot)
    else if (derv == 4 .and. field == 'NONE') then
       call quarticff_nuc(na, nucpot)
    else
       print *,'rsp_backend nuclear_potential error: not implented, derv =', &
               derv, '  field = ', field
       call quit('rsp_backend nuclear_potential error: not implented')
    end if
  end subroutine


  !> Nuclear repulsion contribution to cubic force field, where \param na
  !> is the number of atoms, \param chg the charges of the nuclei, and
  !> \param cor the coordinates of the nuclei
  subroutine cubicff_nuc(na, cub)
     integer, intent(in)  :: na
     real(8), intent(out) :: cub(3,na,3,na,3,na)
     real(8) r(3), rrr(3,3,3), trace(3)
     integer i, j, k, l
     cub = 0 !start from zero
     ! loop over pairs i<j of nuclei
     do j = 2, na
        do i = 1, j-1
           ! construct triple tensor prod with traces removed
           r(1) = get_nuc_xyz(1, j) - get_nuc_xyz(1, i)
           r(2) = get_nuc_xyz(2, j) - get_nuc_xyz(2, i)
           r(3) = get_nuc_xyz(3, j) - get_nuc_xyz(3, i)
           rrr = reshape((/((r(:)*r(k)*r(l),k=1,3),l=1,3)/),(/3,3,3/))
           trace = rrr(:,1,1) + rrr(:,2,2) + rrr(:,3,3)
           do k = 1, 3
              rrr(:,k,k) = rrr(:,k,k) - trace/5
              rrr(k,:,k) = rrr(k,:,k) - trace/5
              rrr(k,k,:) = rrr(k,k,:) - trace/5
           end do
           ! apply scale factor 15 Qi Qj / r^7
           rrr = rrr * (15 * get_nuc_charge(i) * get_nuc_charge(j) / sum(r**2)**(7/2d0))
           cub(:,i,:,i,:,i) = cub(:,i,:,i,:,i) + rrr !iii
           cub(:,i,:,i,:,j) = -rrr !iij
           cub(:,i,:,j,:,i) = -rrr !iji
           cub(:,j,:,i,:,i) = -rrr !jii
           cub(:,i,:,j,:,j) =  rrr !ijj
           cub(:,j,:,i,:,j) =  rrr !jij
           cub(:,j,:,j,:,i) =  rrr !jji
           cub(:,j,:,j,:,j) = cub(:,j,:,j,:,j) - rrr !jjj
        end do
     end do
  end subroutine


  !> Nuclear repulsion contribution to quartic force field, where \param na
  !> is the number of atoms, \param chg the charges of the nuclei, and
  !> \param cor the coordinates of the nuclei
  !> Note: Dimensions 7 and 8 of 'qua' have been joined, to comply with
  !> fortran standard (max 7 dims)
  subroutine quarticff_nuc(na, qua)
     integer, intent(in)  :: na
     real(8), intent(out) :: qua(3,na,3,na,3,na,3*na)
     real(8) r(3), rrrr(3,3,3,3), trace(3,3), trace2
     integer i, j, k, l, m
     qua = 0 !start from zero
     ! loop over pairs i<j of nuclei
     do j = 2, na
        do i = 1, j-1
           ! construct quadruple tensor product with traces removed
           r(1) = get_nuc_xyz(1, j) - get_nuc_xyz(1, i)
           r(2) = get_nuc_xyz(2, j) - get_nuc_xyz(2, i)
           r(3) = get_nuc_xyz(3, j) - get_nuc_xyz(3, i)
           rrrr = reshape((/(((r(:)*r(k)*r(l)*r(m), &
                          k=1,3),l=1,3),m=1,3)/),(/3,3,3,3/))
           trace  = rrrr(:,:,1,1) + rrrr(:,:,2,2) + rrrr(:,:,3,3)
           trace2 = trace(1,1) + trace(2,2) + trace(3,3)
           do k = 1, 3
              trace(k,k) = trace(k,k) - trace2/10
           end do
           do k = 1, 3
              rrrr(:,:,k,k) = rrrr(:,:,k,k) - trace/7
              rrrr(:,k,:,k) = rrrr(:,k,:,k) - trace/7
              rrrr(k,:,:,k) = rrrr(k,:,:,k) - trace/7
              rrrr(:,k,k,:) = rrrr(:,k,k,:) - trace/7
              rrrr(k,:,k,:) = rrrr(k,:,k,:) - trace/7
              rrrr(k,k,:,:) = rrrr(k,k,:,:) - trace/7
           end do
           ! apply scale factor 105 Qi Qj / r^9
           rrrr = rrrr * (105 * get_nuc_charge(i) * get_nuc_charge(j) / sum(r**2)**(9/2d0))
           qua(:,i,:,i,:,i,3*i-2:3*i) = qua(:,i,:,i,:,i,3*i-2:3*i) + rrrr !iiii
           qua(:,i,:,i,:,i,3*j-2:3*j) = -rrrr !iiij
           qua(:,i,:,i,:,j,3*i-2:3*i) = -rrrr !iiji
           qua(:,i,:,j,:,i,3*i-2:3*i) = -rrrr !ijii
           qua(:,j,:,i,:,i,3*i-2:3*i) = -rrrr !jiii
           qua(:,i,:,i,:,j,3*j-2:3*j) =  rrrr !iijj
           qua(:,i,:,j,:,i,3*j-2:3*j) =  rrrr !ijij
           qua(:,j,:,i,:,i,3*j-2:3*j) =  rrrr !jiij
           qua(:,i,:,j,:,j,3*i-2:3*i) =  rrrr !ijji
           qua(:,j,:,i,:,j,3*i-2:3*i) =  rrrr !jiji
           qua(:,j,:,j,:,i,3*i-2:3*i) =  rrrr !jjii
           qua(:,i,:,j,:,j,3*j-2:3*j) = -rrrr !ijjj
           qua(:,j,:,i,:,j,3*j-2:3*j) = -rrrr !jijj
           qua(:,j,:,j,:,i,3*j-2:3*j) = -rrrr !jjij
           qua(:,j,:,j,:,j,3*i-2:3*i) = -rrrr !jjji
           qua(:,j,:,j,:,j,3*j-2:3*j) = qua(:,j,:,j,:,j,3*j-2:3*j) + rrrr !jjjj
        end do
     end do
  end subroutine



  ! derivatives of nuclear repulsion energy
  subroutine nucrep_deriv(na, charges, coords, deriv, nucrep)
    integer, intent(in)  :: na, deriv, charges(na)
    real(8), intent(in)  :: coords(3,na)
    real(8), intent(out) :: nucrep((3*na)**deriv)
    !--------------------------------------------
    real(8) displ(3), power(3**deriv)
    integer nuci, nucj, k, l, m, ncor_pow(deriv)
    ! powers of ncor=3*na, for navigating through nucrep
    ncor_pow(:) = (/((3*na)**k, k=0,deriv-1)/)
    ! start from zero 
    nucrep(:) = 0
    ! loop over pairs i<j of nuclei
    do nucj = 2, na
       do nuci = 1, nucj-1
          ! internuclear displacement
          displ(:) = coords(:,nucj) - coords(:,nuci)
          ! first unperturbed repulsion energy
          power(1) = charges(nuci) * charges(nucj) / sqrt(sum(displ**2))
          ! construct deriv'th tensor power of displacement vector,
          ! carrying the prefactor QiQj (2k-1)!! / r^(2k+1)
          do k = 1, deriv
             l = 3**(k-1)
             power(2*l+1:3*l) = power(:l) * displ(3) * (2*k-1) / sum(displ**2)
             power(  l+1:2*l) = power(:l) * displ(2) * (2*k-1) / sum(displ**2)
             power(     :  l) = power(:l) * displ(1) * (2*k-1) / sum(displ**2)
          end do
          ! remove traces from displacement power
          if (deriv >= 2) call remove_traces(deriv, power)
          ! add traceless displacement power at the 2^deriv different
          ! places in the nucrep derivative tensor
          call place_power_tensor
       end do
    end do

  contains

    recursive subroutine remove_traces(m,a)
      integer, intent(in)    :: m
      real(8), intent(inout) :: a(3**(m-2),3,3)
      real(8) b(3**(m-2))
      integer i, j
      b(:) = a(:,1,1) + a(:,2,2) + a(:,3,3)
      if (m >= 4) call remove_traces(m-2,b)
      do j = 2, m
         do i = 1, j-1
            call subtract_from_diag( &
                         (deriv+m-1) * (1+(deriv-m)/2), &
                         3**(i-1), 3**(j-i-1), 3**(m-j), b, a)
         end do
      end do
    end subroutine

    subroutine subtract_from_diag(denom, l, m, t, b, a)
      integer, intent(in)    :: denom, l, m, t
      real(8), intent(in)    :: b(l,m,t)
      real(8), intent(inout) :: a(l,3,m,3,t)
      a(:,1,:,1,:) = a(:,1,:,1,:) - b(:,:,:)/denom
      a(:,2,:,2,:) = a(:,2,:,2,:) - b(:,:,:)/denom
      a(:,3,:,3,:) = a(:,3,:,3,:) - b(:,:,:)/denom
    end subroutine
    
    subroutine place_power_tensor
      integer offset, i, j, y, z, k, s
      logical neg
      offset = (3*(nuci-1)) * sum(ncor_pow)
      neg = .false.
      ! loop over all index permutations iiiii jiiii ... ijjjj jjjjj
      do k = 0, 2**deriv - 1
         ! loop over Cartesian indices xxxxx yxxxx ... yzzzz zzzzz
         i = offset + 1
         y = 0
         z = 0
         do j = 1, 3**deriv
            ! add or set
            if (k==0 .or. k == 2**deriv-1) then
               nucrep(i) = nucrep(i) + merge(-1,1,neg) * power(j)
            else
               nucrep(i) = merge(-1,1,neg) * power(j)
            end if
            ! increment index i
            do s = 0, deriv-1
               ! x -> y transition
               if (.not.btest(y,s) .and. .not.btest(z,s)) then
                  y = ibset(y,s)
                  i = i + ncor_pow(s+1)
                  exit
               ! y -> z transition
               else if (btest(y,s)) then
                  y = ibclr(y,s)
                  z = ibset(z,s)
                  i = i + ncor_pow(s+1)
                  exit
               ! z -> x transition (followed by 'carry')
               else
                  z = ibclr(z,s)
                  i = i - 2*ncor_pow(s+1)
               end if
            end do
         end do
         ! increment offset ofs for next index permutation
         do s = 0, deriv-1
            neg = .not.neg
            if (.not.btest(k,s)) then
               offset = offset + 3*(nucj-nuci) * ncor_pow(s+1)
               exit
            else
               offset = offset - 3*(nucj-nuci) * ncor_pow(s+1)
            end if
         end do
      end do
    end subroutine

  end subroutine


end module
