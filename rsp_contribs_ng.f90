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
  use matrix_genop_ng, only: mat_alloc

#ifdef BUILD_AORSP
  use dalton_ifc_ng, only: di_read_operator_int, &
                           di_get_gmat, &
                           GRADNN_ifc, &
                           HESSNN_ifc
  use basis_set,  only: cgto
#endif

  implicit none
  public rsp_nucpot
  public rsp_ovlave
  public rsp_oneave
  public rsp_twoave
  public rsp_excave
  public rsp_ovlint
  public rsp_oneint
  public rsp_twoint
  public rsp_excint
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
     !> whether energy/response has exchange-correlation contributin
     logical      :: hasxc
#ifdef BUILD_AORSP
     !> nuclear charges (natom)
     real(8),    pointer :: charge(:)
     !> nuclear coordinates (3,natom)
     real(8),    pointer :: coord(:,:)
     !> basis
     type(cgto), pointer :: basis(:)
#endif
  end type


  !> private struct to collect properties of perturbing "fields"
  type fld_info
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
  type(fld_info) :: field_list(12) = &                         !nc an ba ln qu
     (/fld_info('EXCI', 'Generalized "excitation" field'      , 1, F, F, T, T), &
       fld_info('FREQ', 'Generalized "freqency" field'        , 1, F, F, T, T), &
       fld_info('EL  ', 'Electric field'                      , 3, F, F, T, F), &
       fld_info('VEL ', 'Velocity'                            , 3, T, F, T, F), &
       fld_info('MAGO', 'Magnetic field w/o. London orbitals' , 3, T, F, F, T), &
       fld_info('MAG ', 'Magnetic field with London orbitals' , 3, T, T, F, F), &
       fld_info('ELGR', 'Electric field gradient'             , 6, F, F, T, F), &
       fld_info('VIBM', 'Displacement along vibrational modes',-1, F, T, F, F), &
       fld_info('GEO ', 'Nuclear coordinates'                 ,-1, F, T, F, F), & !-1=mol-dep
       fld_info('NUCM', 'Nuclear magnetic moment'             ,-1, F, T, F, T), & !-1=mol-dep
       fld_info('AOCC', 'AO contraction coefficients'         ,-1, F, T, F, F), & !-1=mol-dep
       fld_info('AOEX', 'AO exponents'                        ,-1, F, T, F, F)/)  !-1=mol-dep

  character(1), parameter :: xyz(3) = (/'X','Y','Z'/)

  private

contains


  !> nuclear repulsion and nuclei--field interaction
  subroutine rsp_nucpot(mol, nf, f, w, c, nc, pot)
    !> structure containing integral program settings
    type(rsp_cfg), intent(in)  :: mol
    !> number of fields
    integer,       intent(in)  :: nf
    !> field labels in std order
    character(4),  intent(in)  :: f(nf)
    !> field frequencies corresponding to each field
    complex(8),    intent(in)  :: w(nf)
    !> first and number of- components in each field
    integer,       intent(in)  :: c(nf), nc(nf)
    !> output average
    complex(8),    intent(out) :: pot(product(nc))
    !----------------------------------------------
    real(8) :: tmp(3*mol%natom, 3*mol%natom, 3*mol%natom)
    integer :: i
    if (nf==0) then
       call quit('rsp_nucpot error: unperturbed (nf=0) nuc.rep. not implemented')
    else if (nf==1 .and. f(1)=='GEO ') then
       call GRADNN_ifc(mol%natom, tmp(:,1,1))
       pot(:nc(1)) = tmp(c(1):c(1)+nc(1)-1,1,1)
    else if (nf==2 .and. all(f==(/'GEO ','GEO '/))) then
       call HESSNN_ifc(mol%natom, tmp(:,:,1))
       pot = reshape(tmp(c(1):c(1)+nc(1)-1, &
                         c(2):c(2)+nc(2)-1,1), shape(pot))
    else if (nf==3 .and. all(f==(/'GEO ','GEO ','GEO '/))) then
       call cubicff_nuc(mol%natom, mol%charge, mol%coord, tmp)
       pot = reshape(tmp(c(1):c(1)+nc(1)-1, &
                         c(2):c(2)+nc(2)-1, &
                         c(3):c(3)+nc(3)-1), shape(pot))
    else
       print *,'rsp_nucpot error: not implented or in wrong order - ', &
               (' ' // f(i), i=1,nf)
       call quit('rsp_nucpot error: not implented or in wrong order')
    end if
  end subroutine



  !> average f-perturbed overlap integrals with perturbed density D
  !> and energy-weighted density DFD
  subroutine rsp_ovlave(mol, nf, f, c, nc, DFD, ave, w, D)
    use dalton_ifc_ng, only: SHELLS_NUCLEI_displace
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
    !----------------------------------------------
    type(matrix) A(2)
    real(8), parameter :: fdistep = 2d0**(-25)
    real(8) tmp(3*mol%natom, 3*mol%natom, 3*mol%natom)
    real(8) fdi(3*mol%natom, 3*mol%natom)
    integer i
    if (present(w) .and. .not.present(D)) &
       call quit("error in rsp_ovlave: frequencies 'w' and density 'D' " &
              // 'must both be present or both absent')
    if (nf==0) then
       call quit('rsp_ovlave error: unperturbed (nf=0) overlap not implemented')
    else if (nf==1 .and. f(1)=='GEO ') then
       ! allocate matrices for integrals
       A(1) = mol%zeromat
       call mat_alloc(A(1))
       if (present(w)) then
          A(2) = mol%zeromat
          call mat_alloc(A(2))
       end if
       ! loop over nuclear coordinates
       do i = 0, nc(1)-1
          ! (half-) perturbed overlap -i/2 Tg into A(1), Sg in A(2)
          if (present(w)) then !w=0 means no -i/2 Tg contribution
             call di_read_operator_int('SQHDR' // prefix_zeros(c(1)+i,3), A(1))
             A(1) = -A(1) !SQHDR is really -dS>/dg
             A(2) = (-w(1)/2) * (A(1) + trps(A(1)))
             A(1) = A(1) + trps(A(1)) !=1DOVL
          else
             call di_read_operator_int('1DOVL' // prefix_zeros(c(1)+i,3), A(1))
             A(1) = -A(1) !1DOVL is really -dS/dg
          end if
          ave(1+i) = -tr(A(1),DFD)
          if (present(w)) ave(1+i) = ave(1+i) + tr(A(2),D)
       end do
       A(1:2) = 0 !deallocate
    else if (nf==2 .and. all(f==(/'GEO ','GEO '/))) then
       if (present(w)) then
          if (.not.all(w==0)) &
             call quit('rsp_ovlave error: GEO-GEO with freqencies not implemented')
       end if
       call ONEDRV_ave_ifc(mol, f, size(tmp(:,:,1)), tmp(:,:,1), DFD=DFD)
       ave = reshape(tmp(c(1):c(1)+nc(1)-1, &
                         c(2):c(2)+nc(2)-1,1), shape(ave))
    else if (nf==3 .and. all(f==(/'GEO ','GEO ','GEO '/))) then
       if (present(w)) then
          if (.not.all(w==0)) &
             call quit('rsp_ovlave error: GEO-GEO-GEO with freqencies not implemented')
       end if
       ! gen1int_driver(prop_name, order_geo, order_mag,
       !                is_lao, get_int, wrt_int, vals_int,
       !                do_exp, ndens, ao_dens,
       !                get_exp, wrt_exp, vals_expect,
       !                len_work, dal_work, level_print)
       !ajt FIXME use finite difference intil Gao has committed update
       !call gen1int_driver('OVERLAP', 3, 0,
       !                    .false., .false., .false., tmp,
       !                    .true., 1, D%elms,
       !                    .true., .false., tmp,
       !                    size(work), work, 5)
       do i = 0, nc(3)-1
          call SHELLS_NUCLEI_displace(c(3)+i, fdistep)
          call ONEDRV_ave_ifc(mol, f(1:2), size(fdi), fdi, DFD=DFD)
          tmp(:,:,c(3)+i) = fdi / (2*fdistep)
          call SHELLS_NUCLEI_displace(c(3)+i, -2*fdistep)
          call ONEDRV_ave_ifc(mol, f(1:2), size(fdi), fdi, DFD=DFD)
          tmp(:,:,c(3)+i) = tmp(:,:,c(3)+i) - fdi / (2*fdistep)
          call SHELLS_NUCLEI_displace(c(3)+i, fdistep)
       end do
       ave = reshape(tmp(c(1):c(1)+nc(1)-1, &
                         c(2):c(2)+nc(2)-1, &
                         c(3):c(3)+nc(3)-1), shape(ave))
    else
       print *, 'rsp_ovlave error: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('rsp_ovlave error: not implented or in wrong order')
    end if
  end subroutine



  !> Average 1-electron integrals perturbed by fields f
  !> with the (perturbed) density matrix D
  subroutine rsp_oneave(mol, nf, f, c, nc, D, ave)
    use dalton_ifc_ng, only: SHELLS_NUCLEI_displace, &
                             work => dal_work
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
    !----------------------------------------------
    ! step length for finite difference. Power of 2 to ensure
    ! addition and subsequent subtraction yields original
    type(matrix) A(1)
    real(8), parameter :: fdistep = 2d0**(-25)
    real(8) tmpggg(3*mol%natom, 3*mol%natom, 3*mol%natom)
    real(8) fdigg(3*mol%natom, 3*mol%natom)
    real(8) tmpggf(3*mol%natom, 3*mol%natom, 3)
    real(8) fdigf(3*mol%natom, 3)
    integer i, j, k, n, indx
    if (nf==0) then
       call quit('rsp_oneave error: unperturbed (nf=0) 1el.int. not implemented')
    else if (nf==1 .and. f(1)=='GEO ') then
       A(1) = mol%zeromat
       call mat_alloc(A(1))
       do i = 0, nc(1)-1
          ! perturbed one-electron Hamiltonian integrals
          call di_read_operator_int('1DHAM' // prefix_zeros(c(1)+i,3), A(1))
          ave(1+i) = tr(A(1),D)
       end do
       A(1) = 0
    else if (nf==1 .and. f(1)=='EL  ') then
       A(1) = mol%zeromat
       call mat_alloc(A(1))
       do i = 0, nc(1)-1
          ! dipole integrals
          call di_read_operator_int(xyz(c(1)+i) // 'DIPLEN ', A(1))
          ave(1+i) = tr(A(1),D)
       end do
       A(1) = 0
    else if (nf==2 .and. all(f==(/'GEO ','EL  '/))) then
       ! put averages into tmp(:,:3,1)
       ! gen1int_driver('CARMOM', 1, 0,
       !                .false., .false., .false., tmp(:1,1,1),
       !                .true., 1, D%elms,
       !                .true., .false., tmp(:,:3,1),
       !                size(work), work, 5)
       !ave = reshape(tmp(c(1):c(1)+nc(1)-1, &
       !                  c(2):c(2)+nc(2)-1,1), shape(ave))
       ! read integrals from AOPROPER until gen1int works
       A(1) = mol%zeromat
       call mat_alloc(A(1))
       do j = 0, nc(2)-1 ! EL
          do i = 0, nc(1)-1 ! GEO
             call di_read_operator_int( &
                  prefix_zeros(c(1)+i,3) // 'DPG ' // xyz(c(2)+j), A(1))
                  ave(1+i+nc(1)*j) = -tr(A(1),D)
          end do
       end do
       A(1) = 0
    else if (nf==2 .and. all(f==(/'GEO ','GEO '/))) then
       call ONEDRV_ave_ifc(mol, f, size(tmpggg(:,:,1)), tmpggg(:,:,1), D=D)
       ave = reshape(tmpggg(c(1):c(1)+nc(1)-1, &
                            c(2):c(2)+nc(2)-1,1), shape(ave))
    else if (nf==3 .and. all(f==(/'GEO ','GEO ','EL  '/))) then
       !dj FIXME use finite difference until Gao has committed update
       fdigf = 0
       do i = 0, nc(1)-1
          call SHELLS_NUCLEI_displace(c(1)+i, fdistep)
          fdigf = 0
          call GET1IN_ave_ifc(mol, f(2:3), size(fdigf), fdigf, D)
          tmpggf(c(1)+i,:,:) = fdigf / (2*fdistep)
          call SHELLS_NUCLEI_displace(c(1)+i, -2*fdistep)
          fdigf = 0
          call GET1IN_ave_ifc(mol, f(2:3), size(fdigf), fdigf, D)
          tmpggf(c(1)+i,:,:) = tmpggf(c(1)+i,:,:) - fdigf / (2*fdistep)
          call SHELLS_NUCLEI_displace(c(1)+i, fdistep)
          !tmpggf(c(1)+i,:,:) = transpose(tmpggf(c(1)+i,:,:))
       end do
       indx = 1
       do k = 1, nc(3)
          do i = 1, nc(3)
             do j = k, nc(2), nc(3)
                ave(indx:indx+nc(1)) = -2d0*tmpggf(:,j,i)
                indx = indx + nc(1)
             end do
          end do
       end do
       !ave = reshape(-2d0*tmpggf(c(1):c(1)+nc(1)-1, &
       !                          c(2):c(2)+nc(2)-1, &
       !                          c(3):c(3)+nc(3)-1), shape(ave))
    else if (nf==3 .and. all(f==(/'GEO ','GEO ','GEO '/))) then
       ! gen1int_driver(prop_name, order_geo, order_mag,
       !                is_lao, get_int, wrt_int, vals_int,
       !                do_exp, ndens, ao_dens,
       !                get_exp, wrt_exp, vals_expect,
       !                len_work, dal_work, level_print)
       !ajt FIXME use finite difference intil Gao has committed update
       !call gen1int_driver('ONEHAMIL', 3, 0,
       !                    .false., .false., .false., tmp,
       !                    .true., 1, D%elms,
       !                    .true., .false., tmp,
       !                    size(work), work, 5)
       do i = 0, nc(3)-1
          call SHELLS_NUCLEI_displace(c(3)+i, fdistep)
          call ONEDRV_ave_ifc(mol, f(1:2), size(fdigg), fdigg, D=D)
          tmpggg(:,:,c(3)+i) = fdigg / (2*fdistep)
          call SHELLS_NUCLEI_displace(c(3)+i, -2*fdistep)
          call ONEDRV_ave_ifc(mol, f(1:2), size(fdigg), fdigg, D=D)
          tmpggg(:,:,c(3)+i) = tmpggg(:,:,c(3)+i) - fdigg / (2*fdistep)
          call SHELLS_NUCLEI_displace(c(3)+i, fdistep)
       end do
       ave = reshape(tmpggg(c(1):c(1)+nc(1)-1, &
                            c(2):c(2)+nc(2)-1, &
                            c(3):c(3)+nc(3)-1), shape(ave))
    else
       print *, 'rsp_oneave error: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('rsp_oneave error: not implented or in wrong order')
    end if
  end subroutine



  !> Average 2-electron integrals perturbed by fields f over the
  !> product of density matrces D1 and D2
  subroutine rsp_twoave(mol, nf, f, c, nc, D1, D2, ave)
    use dalton_ifc_ng, only: work => dal_work
    use cgto_diff_eri, only: geo_eri
    use contract_eri,  only: coul_exch_ave, geo_eri_loop_unopt
    !> structure containing integral program settings
    type(rsp_cfg),   intent(in) :: mol
    !> number of fields
    integer,         intent(in) :: nf
    !> field labels in std order
    character(4),    intent(in) :: f(nf)
    !> first and number of- components in each field
    integer,         intent(in) :: c(nf), nc(nf)
    !> density matrix to average over
    type(matrix),    intent(in) :: D1, D2
    !> output average
    complex(8),     intent(out) :: ave(product(nc))
    !----------------------------------------------
    real(8)      tmp(3*mol%natom, 3*mol%natom, 3*mol%natom), r !scratch
    type(matrix) A(1) !scratch matrices
    integer      i, j, k, n
    if (nf==0) then
       ! contract second density to Fock matrix, then trace with first
       A(1) = mol%zeromat
       call mat_alloc(A(1))
       call di_get_gmat(D2, A(1)) !Coulomb and exchange
       ave(1) = tr(A(1),D1)
    else if (nf==1 .and. f(1)=='GEO ') then
       n = mol%zeromat%nrow
       work(     :n*n)   = reshape(D1%elms,(/n*n/))
       work(n*n+1:n*n*2) = reshape(D2%elms,(/n*n/))
       call GRCONT(work(n*n*2+1:), size(work)-n*n*2, &
                   tmp(:,1,1), 3*mol%natom, .true., .false., &
                   1, 0, .true., .false., work(:n*n*2), 2)
       ave(:nc(1)) = tmp(c(1):c(1)+nc(1)-1,1,1)
    else if (nf==2 .and. all(f==(/'GEO ','GEO '/))) then
       n = mol%zeromat%nrow
       work(     :n*n)   = reshape(D1%elms,(/n*n/))
       work(n*n+1:n*n*2) = reshape(D2%elms,(/n*n/))
       call GRCONT(work(n*n*2+1:), size(work)-n*n*2, &
                   tmp(:,:,1), (3*mol%natom)**2, .true., .false., &
                   2, 0, .true., .false., work(:n*n*2), 2)
       ave = reshape(tmp(c(1):c(1)+nc(1)-1, &
                         c(2):c(2)+nc(2)-1,1), shape(ave))
    else if (nf==3 .and. all(f==(/'GEO ','GEO ','GEO '/))) then
       if (.not.associated(D1%elms,D2%elms)) &
          call quit('rsp_twoave error: FIXME')
       ! contract FULL cubic
       call geo_eri_loop_unopt(mol%natom, size(mol%basis), mol%basis, 3, &
                               coul_exch_ave, D1%elms, mol%natom, tmp)
       ! symmetrize
       do k = 1, 3*mol%natom
          do j = 1,k
             do i = 1,j
                r = tmp(i,j,k) + tmp(i,k,j) + tmp(k,i,j) &
                  + tmp(k,j,i) + tmp(j,k,i) + tmp(j,i,k)
                tmp(i,j,k) = r;  tmp(i,k,j) = r;  tmp(k,i,j) = r
                tmp(k,j,i) = r;  tmp(j,k,i) = r;  tmp(j,i,k) = r
             end do
          end do
       end do
       ! extract designated block
       ave = 2 * reshape(tmp(c(1):c(1)+nc(1)-1, &
                             c(2):c(2)+nc(2)-1, &
                             c(3):c(3)+nc(3)-1), shape(ave))
    else
       print *, 'rsp_twoave error: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('rsp_twoave error: not implented or in wrong order')
    end if
  end subroutine



  !> Exchange-correlation perturbed by fields f, averaged over densities D
  subroutine rsp_excave(mol, nf, f, c, nc, nd, D, ave)
    !> structure containing integral program settings
    type(rsp_cfg),   intent(in) :: mol
    !> number of fields
    integer,         intent(in) :: nf
    !> field labels in std order
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
    !ajt FIXME currently unimplemented
    if (mol%hasxc) call quit('error in rsp_excave: not implemented')
    ave = 0
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
    !------------------------------------------------
    integer      i
    if (present(w) .and. .not.present(fock)) &
       call quit("error in rsp_ovlint: frequencies 'w' and Fock matrix 'fock' " &
              // 'must both be present or both absent')
    if (nf==0) then
       call quit('rsp_ovlint error: unperturbed (nf=0) overlap not implemented')
    else if (nf==1 .and. f(1)=='GEO ') then
       ! loop over nuclear coordinates
       do i = 0, nc(1)-1
          ! allocate, if needed
          if (.not.isdef(ovl(1+i))) then
             ovl(1+i) = mol%zeromat
             call mat_alloc(ovl(1+i))
          end if
          ! overlap into ovl, half-perturbed overlap -i/2 Tg added to fock
          if (present(w)) then
             call di_read_operator_int('SQHDR' // prefix_zeros(c(1)+i,3), ovl(1+i))
             ovl(1+i)  = -ovl(1+i) !SQHDR is really -dS>/dg
             fock(1+i) = fock(1+i) - w(1)/2 * ovl(1+i)
             fock(1+i) = fock(1+i) + w(1)/2 * trps(ovl(1+i))
             ovl(1+i)  = ovl(1+i)  + trps(ovl(1+i)) !=dS/dg=-1DOVL
          else
             call di_read_operator_int('1DOVL' // prefix_zeros(c(1)+i,3), ovl(1+i))
             ovl(1+i) = -ovl(1+i) !1DOVL is really -dS/dg
          end if
       end do
    else
       print *, 'rsp_ovlint error: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('rsp_ovlint error: not implented or in wrong order')
    end if
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
    !--------------------------------------------------
    integer      i
    if (nf==0) then
       call quit('error in rsp_oneint: unperturbed (nf=0) 1el.int. not implemented')
    else if (nf==1 .and. f(1)=='GEO ') then
       do i = 0, nc(1)-1
          ! allocate, if needed
          if (.not.isdef(oneint(1+i))) then
             oneint(1+i) = mol%zeromat
             call mat_alloc(oneint(1+i))
          end if
          ! perturbed one-electron Hamiltonian integrals
          call di_read_operator_int('1DHAM' // prefix_zeros(c(1)+i,3), oneint(1+i))
       end do
    else if (nf==1 .and. f(1)=='EL  ') then
       do i = 0, nc(1)-1
          ! allocate, if needed
          if (.not.isdef(oneint(1+i))) then
             oneint(1+i) = mol%zeromat
             call mat_alloc(oneint(1+i))
          end if
          ! dipole integrals
          call di_read_operator_int(xyz(c(1)+i) // 'DIPLEN ', oneint(1+i))
       end do
    else
       print *, 'error in rsp_oneave: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('error in rsp_oneave: not implented or in wrong order')
    end if
  end subroutine



  !> Contract 2-electron integrals perturbed by fields 'f' with density
  !> matrix 'dens', and add to Fock matrices 'fock' Average 2-electron integrals perturbed by fields f over the
  !> product of density matrces D1 and D2
  subroutine rsp_twoint(mol, nf, f, c, nc, dens, fock)
    ! work array to be passed to GRCONT
    use dalton_ifc_ng, only: work => dal_work
    !> structure containing integral program settings
    type(rsp_cfg),   intent(in)    :: mol
    !> number of fields
    integer,         intent(in)    :: nf
    !> field labels in std order
    character(4),    intent(in)    :: f(nf)
    !> first and number of- components in each field
    integer,         intent(in)    :: c(nf), nc(nf)
    !> density matrix to average over
    type(matrix),    intent(in)    :: dens
    !> Fock matrix to which two-electron contribution is ADDED
    type(matrix),    intent(inout) :: fock(product(nc))
    !--------------------------------------------------
    integer      i, j, n
    type(matrix) A !scratch
    if (nf==0) then
       A = 0*dens
       call mat_alloc(A)
       call di_get_gmat(dens, A)
       fock(1) = fock(1) + A
       A = 0
    else if (nf==1 .and. f(1)=='GEO ') then
       n = mol%zeromat%nrow
       do i = 0, nc(1)-1
          ! if first or an x-coord, call GRCONT
          if (i==0 .or. mod(c(1)+i,3) == 1) &
             call GRCONT(work(n*n*3+1:), size(work)-n*n*3, &
                         work(:n*n*3), n*n*3, .true., .false., &
                         1, (c(1)+i+2)/3, .false., .true., dens%elms, 1)
          j = 1 + mod(c(1)+i-1,3) !x y z = 1 2 3
          if (iszero(fock(1+i))) then
             call mat_alloc(fock(1+i))
             fock(1+i)%elms = reshape(work(n*n*(j-1)+1:n*n*j),(/n,n/))
          else
             fock(1+i)%elms = fock(1+i)%elms &
                            + reshape(work(n*n*(j-1)+1:n*n*j),(/n,n/))
          end if
       end do
    else
       print *, 'error in rsp_oneave: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('error in rsp_oneave: not implented or in wrong order')
    end if
  end subroutine



  !> Exchange-correlation perturbed by fields 'f', contracted over
  !> densities 'D', added to Fock matrices 'fock'
  subroutine rsp_excint(mol, nf, f, c, nc, nd, D, fock)
    !> structure containing integral program settings
    type(rsp_cfg),  intent(in)    :: mol
    !> number of fields
    integer,        intent(in)    :: nf
    !> field labels in std order
    character(4),   intent(in)    :: f(nf)
    !> first and number of- components in each field
    integer,        intent(in)    :: c(nf), nc(nf)
    !> number of density matrices
    integer,        intent(in)    :: nd
    !> density matrices to average over. D(1) must be unperturbed
    type(matrix),   intent(in)    :: D(nd)
    !> output average
    type(matrix),   intent(inout) :: fock(product(nc))
    !----------------------------------------------
    !ajt FIXME currently unimplemented
    if (mol%hasxc) call quit('error in rsp_excave: not implemented')
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
    !> field labels
    character(4),  intent(in) :: f(:)
    integer :: rsp_field_dim(size(f)), i
    rsp_field_dim = (/(field_list(idx(f(i)))%ncomp, i=1,size(f))/)
    ! loop through mol-dependent
    do i = 1, size(f)
       if (rsp_field_dim(i) /= -1) then
          ! cycle
       else if (f(i)=='GEO ') then
          rsp_field_dim(i) = 3 * mol%natom
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
  subroutine ONEDRV_ave_ifc(mol, fld, siz, ave, D, DFD)
    use dalton_ifc_ng, only: dal_work
    type(rsp_cfg),          intent(in)  :: mol
    character(4),           intent(in)  :: fld(:)
    integer,                intent(in)  :: siz
    real(8),                intent(out) :: ave(siz)
    type(matrix), optional, intent(in)  :: D, DFD
    !--------------------------------------------
#include <mxcent.h>
#include <taymol.h>
    real(8) Dtri(  (mol%zeromat%nrow)*(mol%zeromat%nrow+1)/2)
    real(8) DFDtri((mol%zeromat%nrow)*(mol%zeromat%nrow+1)/2)
    logical anti
    integer nc
    ! create triangularly packed matrices from D, DFD
    anti = (mod(count(fld=='MAG '),2) == 1)
    if (.not.present(D)) then
       Dtri = 0
    else if (anti) then
       call DGETAP(D%nrow, D%elms, Dtri)
    else !symm
       call DGEFSP(D%nrow, D%elms, Dtri)
    end if
    if (.not.present(DFD)) then
       DFDtri = 0
    else if (anti) then
       call DGETAP(DFD%nrow, DFD%elms, DFDtri)
    else !symm
       call DGEFSP(DFD%nrow, DFD%elms, DFDtri)
    end if
    ! write to files
    call WRITE_DSOFSO(Dtri, DFDtri)
    nc = 3 * mol%natom
    HESMOL(:nc,:nc) = 0
    !  SUBROUTINE ONEDRV(WORK,LWORK,IPRINT,PROPTY,MAXDIF,
    ! &                  DIFINT,NODC,NODV,DIFDIP,DIFQDP,
    ! &                  HFONLY,NCLONE,PCM)
    call ONEDRV(dal_work, size(dal_work), 5, .true., size(fld), &
                .true., .true., .true., .false., .false., &
                .true., .false., .false.)
    ! HESMOL will contain either HESSKE+HESSNA or HESFS2 or their sum
    if (size(fld)==1) &
       call quit('error in ONEDRV_ave_ifc: not implemented')
    if (size(fld)==2) &
       ave = reshape(2*HESMOL(:nc,:nc), (/nc*nc/)) !factor 2 for total dens
  end subroutine


  !> Call GET1IN in ABACUS
  subroutine GET1IN_ave_ifc(mol, fld, siz, ave, D)
    use dalton_ifc_ng, only: dal_work
    type(rsp_cfg),          intent(in)  :: mol
    character(4),           intent(in)  :: fld(:)
    integer,                intent(in)  :: siz
    real(8),                intent(out) :: ave(siz)
    type(matrix),           intent(in)  :: D
    !--------------------------------------------
#include <dummy.h>
#include <mxcent.h>
    real(8) Dtri(  (mol%zeromat%nrow)*(mol%zeromat%nrow+1)/2)
    character(8), dimension(9*MXCENT) :: labint
    integer, dimension(9*MXCENT) :: intrep, intadr
    integer ncomp
    !ajt Dzero is dangerous! Should have been izero. =0 safe
    intrep = 0 !call dzero(intrep,9*MXCENT)
    intadr = 0 !call dzero(intadr,9*MXCENT)
    ! create triangularly packed matrix from D
    call DGEFSP(D%nrow, D%elms, Dtri)
    ncomp = 0
!      SUBROUTINE GET1IN(SINTMA,WORD,NCOMP,WORK,LWORK,LABINT,INTREP,
!     &                  INTADR,MPQUAD,TOFILE,KPATOM,TRIMAT,EXPVAL,
!     &                  EXP1VL,DENMAT,NPRINT)
    call GET1IN(dummy,'DPLGRA ',ncomp,dal_work,size(dal_work),labint,intrep, &
                       intadr,0,.false.,0,.false.,ave, &
                       .true.,Dtri,0)
  end subroutine


  !> Nuclear repulsion contribution to cubic force field, where \param na
  !> is the number of atoms, \param chg the charges of the nuclei, and
  !> \param cor the coordinates of the nuclei
  subroutine cubicff_nuc(na, chg, cor, cub)
     integer, intent(in)  :: na
     real(8), intent(in)  :: chg(na)
     real(8), intent(in)  :: cor(3,na)
     real(8), intent(out) :: cub(3,na,3,na,3,na)
     real(8) r(3), rrr(3,3,3), tr(3)
     integer i, j, k, l
     cub = 0 !start from zero
     ! loop over pairs i<j of nuclei
     do j = 2,na
	do i = 1,j-1
           ! contrcuct 3rd tensor prod with traces removed
           r = cor(:,j) - cor(:,i)
           rrr = reshape((/((r(:)*r(k)*r(l),k=1,3),l=1,3)/),(/3,3,3/))
           tr = rrr(:,1,1) + rrr(:,2,2) + rrr(:,3,3)
           do k = 1,3
              rrr(:,k,k) = rrr(:,k,k) - tr/5
              rrr(k,:,k) = rrr(k,:,k) - tr/5
              rrr(k,k,:) = rrr(k,k,:) - tr/5
           end do
           ! apply scale factor 15 Qi Qj / r^7
           rrr = rrr * (15 * chg(i) * chg(j) / sum(r**2)**(7/2d0))
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

end module
