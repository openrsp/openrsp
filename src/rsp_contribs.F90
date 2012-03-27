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
  use matrix_backend, only: mat_alloc
  use rsp_backend
  use dalton_ifc, only: di_read_operator_int, &
                        di_get_gmat
  use basis_set,  only: cgto

  implicit none
  public rsp_nucpot
  public rsp_ovlave
  public rsp_oneave
  public rsp_twoave
  public rsp_xcave
  public rsp_xcave_new
  public rsp_ovlint
  public rsp_oneint
  public rsp_twoint
  public rsp_xcint
  public rsp_field
  public rsp_field_bas
  public rsp_field_dim
  public rsp_field_anti
  public rsp_field_ordering
  public rsp_cfg

  !ajt Public for now, for use also in rsp_functions, until interfaces
  !    to integrals are all moved to rsp_backend, and transition of arguments
  !    of type(rsp_field) is complete.
  public perm_comp_add

!FIXME Gao: all Gen1Int routines will be moved to rsp_backend.F90
#ifdef BUILD_GEN1INT
  private rsp_set_geom
#endif

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
    ncor = 3 * get_natom()
    ngeo = 0
    last_field = 'NONE'
    last_ncomp = 1
    !ajt FIXME validate comp/ncomp ranges
    !ajt FIXME determine sorting
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
    use dalton_ifc, only: SHELLS_NUCLEI_displace
#ifdef BUILD_GEN1INT
    ! Gen1Int interface in Dalton
    use gen1int_interface
#endif
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
    real(8), parameter   :: fdistep = 2d0**(-25)
    real(8), allocatable :: tmp(:,:,:)
    real(8), allocatable :: fdi(:,:)
    integer i, ncor
#ifdef BUILD_GEN1INT
    ! origins in Dalton
#include "orgcom.h"
    integer order_geo                    !order of geometric derivatives
    integer num_coord                    !number of atomic coordinates
    integer num_expectations             !number of all expectation values
    real(8), allocatable :: val_expt(:)  !expectation values, real numbers
    type(one_prop_t) overlap             !operator of overlap integrals
    integer ierr                         !error information
#endif
    if (present(w) .and. .not.present(D)) &
       call quit("error in rsp_ovlave: frequencies 'w' and density 'D' " &
              // 'must both be present or both absent')
#ifdef BUILD_GEN1INT
    ! gets the order of geometric derivatives
    order_geo = count(f=='GEO ')
    if (order_geo/=nf) &
      call quit("rsp_ovlave>> only geometric derivatives implemented!")
    ! sets the number of atomic coordinates
    num_coord = 3*get_natom()
    ! allocates memory for expectation values
    num_expectations = num_coord**order_geo
    allocate(val_expt(num_expectations), stat=ierr)
    if (ierr/=0) call quit("rsp_ovlave>> failed to allocate val_expt!")
    ! creates operator for overlap integrals
    call OnePropCreate(prop_name=INT_OVERLAP, one_prop=overlap, &
                       info_prop=ierr, dipole_origin=DIPORG)
    if (ierr/=0) &
      call quit("rsp_ovlave>> failed to create operator of overlap integrals!")
!FIXME changes to call Gen1Int
!FIXME \sum_{j+k=n} (-\sum_{j}w_{j}+\sum_{k}w_{k})/2 S(>)^{j}(<)^{k}
! When it comes to w, I think it's a better idea to ask the integral program for S>> and S<>, and multiply the resulting average or integral by w afterwards, rather than to send w into the integral program.
    ! with field frequencies
    if (present(w)) then
      if (order_geo==1) then
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
      else
        call quit('rsp_ovlave>> GEO(>1) with freqencies not implemented!')
      end if
    else
      ! calculates the expectaion values of derivatives of overlap matrix
!FIXME: rewrites \fn(gen1int_ifc_main) so that we could pass information of basis sets
      val_expt = 0.0
      call gen1int_ifc_main(one_prop=overlap,                          &
                            order_geo_total=order_geo,                 &
                            max_num_cent=min(2,order_geo,get_natom()), &
                            num_ints=num_expectations,                 &
                            num_dens=1, ao_dens=(/DFD/),               &
                            val_expt=val_expt,                         &
                            redunt_expt=order_geo>1,                   &
                            io_viewer=get_lupri(),                     &
                            level_print=5)
      val_expt = -val_expt
      ! assigns the output average
      call rsp_set_geom(mol, nf, c, nc, val_expt, ave)
    end if
    ! frees space
    deallocate(val_expt)
    call OnePropDestroy(one_prop=overlap)
#else
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
       ncor = 3 * get_natom()
       allocate(tmp(ncor,ncor,1))
       call ONEDRV_ave_ifc(mol, f, size(tmp(:,:,1)), tmp(:,:,1), DFD=DFD)
       ave = reshape(tmp(c(1):c(1)+nc(1)-1, &
                         c(2):c(2)+nc(2)-1,1), shape(ave))
       deallocate(tmp)
    else if (nf==3 .and. all(f==(/'GEO ','GEO ','GEO '/))) then
       if (present(w)) then
          if (.not.all(w==0)) &
             call quit('rsp_ovlave error: GEO-GEO-GEO with freqencies not implemented')
       end if
       ncor = 3 * get_natom()
       allocate(tmp(ncor,ncor,ncor))
       allocate(fdi(ncor,ncor))
       do i = 0, nc(3)-1
          call SHELLS_NUCLEI_displace(c(3)+i, fdistep)
          call ONEDRV_ave_ifc(mol, f(1:2), size(fdi), fdi, DFD=DFD)
          tmp(:,:,c(3)+i) = fdi / (2*fdistep)
          call SHELLS_NUCLEI_displace(c(3)+i, -2*fdistep)
          call ONEDRV_ave_ifc(mol, f(1:2), size(fdi), fdi, DFD=DFD)
          tmp(:,:,c(3)+i) = tmp(:,:,c(3)+i) - fdi / (2*fdistep)
          call SHELLS_NUCLEI_displace(c(3)+i, fdistep)
       end do
       deallocate(fdi)
       ave = reshape(tmp(c(1):c(1)+nc(1)-1, &
                         c(2):c(2)+nc(2)-1, &
                         c(3):c(3)+nc(3)-1), shape(ave))
       deallocate(tmp)
    else
       print *, 'rsp_ovlave error: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('rsp_ovlave error: not implented or in wrong order')
    end if
#endif
  end subroutine



  !> Average 1-electron integrals perturbed by fields f
  !> with the (perturbed) density matrix D
  subroutine rsp_oneave(mol, nf, f, c, nc, D, ave)
    use dalton_ifc, only: SHELLS_NUCLEI_displace, dal_work
#ifdef BUILD_GEN1INT
    ! Gen1Int interface in Dalton
    use gen1int_interface
#endif
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
    real(8), parameter   :: fdistep = 2d0**(-25)
    real(8), allocatable :: tmpggg(:,:,:)
    real(8), allocatable :: fdigg(:,:)
    real(8), allocatable :: tmpfg(:,:)
    real(8), allocatable :: tmpfgg(:,:,:)
    real(8), allocatable :: tmpggf(:,:,:)
    real(8), allocatable :: fdigf(:,:)
    integer i, j, k, n, indx, ncor
#ifdef BUILD_GEN1INT
    ! origins in Dalton
#include "orgcom.h"
    ! uses \var(MXCENT)
#include "mxcent.h"
    ! coordinates and charges of atoms
#include "nuclei.h"
    integer order_geo                    !order of geometric derivatives
    integer num_coord                    !number of atomic coordinates
    integer num_expectations             !number of all expectation values
    real(8), allocatable :: val_expt(:)  !expectation values, real numbers
    type(one_prop_t) one_hamil           !operator of one-electron Hamiltonian
    integer ierr                         !error information
!FIXME: implements dipole length, etc.
    ! gets the order of geometric derivatives
    order_geo = count(f=='GEO ')
    if (order_geo/=nf) &
      call quit("rsp_oneave>> only geometric derivatives implemented!")
    ! sets the number of atomic coordinates
    num_coord = 3*get_natom()
    ! allocates memory for expectation values
    num_expectations = num_coord**order_geo
    allocate(val_expt(num_expectations), stat=ierr)
    if (ierr/=0) call quit("rsp_oneave>> failed to allocate val_expt!")
    ! creates operator of one-electron Hamiltonian
    call OnePropCreate(prop_name=INT_ONE_HAMIL, one_prop=one_hamil,   &
                       info_prop=ierr, coord_nuclei=CORD(:,1:NUCDEP), &
                       charge_nuclei=-CHARGE(1:NUCDEP))
    if (ierr/=0) &
      call quit("rsp_oneave>> failed to create operator of one-electron Hamiltonian!")
    ! calculates the expectaion values of derivatives of one-electron Hamiltonian
!FIXME: rewrites \fn(gen1int_ifc_main) so that we could pass information of basis sets
    val_expt = 0.0
    call gen1int_ifc_main(one_prop=one_hamil,                        &
                          order_geo_total=order_geo,                 &
                          max_num_cent=min(3,order_geo,get_natom()), &
                          num_ints=num_expectations,                 &
                          num_dens=1, ao_dens=(/D/),                 &
                          val_expt=val_expt,                         &
                          redunt_expt=order_geo>1,                   &
                          io_viewer=get_lupri(),                     &
                          level_print=5)
    ! assigns the output average
    call rsp_set_geom(mol, nf, c, nc, val_expt, ave)
    ! frees space
    deallocate(val_expt)
    call OnePropDestroy(one_prop=one_hamil)
#else
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
!-temporarily disabled by Bin Gao before implementing operator for dipole length integrals
!-#ifdef BUILD_GEN1INT
!-       ! calculates the first order total geometric derivatives of dipole length integrals
!-       tmpfg = 0.0
!-       call gen1int_ifc_main("DIPLEN", .false., 0, 3, 1,           &
!-                             .false., .false., .false.,            &
!-                             size(tmpfg), tmpfg, 0, D%nrow*D%ncol, &
!-                             1, D%elms, .true., .false., .false.,  &
!-                             tmpfg, get_lupri(), 5)
!-       tmpfg = 2.0*tmpfg
!-       do j = 0, nc(2)-1 ! EL
!-         do i = 0, nc(1)-1 ! GEO
!-           ave(1+i+nc(1)*j) = tmpfg(c(2)+j,c(1)+i)
!-         end do
!-       end do
!-#else
       ! read integrals from AOPROPER
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
!-#endif
    else if (nf==2 .and. all(f==(/'GEO ','GEO '/))) then
       ncor = 3 * get_natom()
       allocate(tmpggg(ncor,ncor,1))
       call ONEDRV_ave_ifc(mol, f, size(tmpggg(:,:,1)), tmpggg(:,:,1), D=D)
       ave = reshape(tmpggg(c(1):c(1)+nc(1)-1, &
                            c(2):c(2)+nc(2)-1,1), shape(ave))
       deallocate(tmpggg)
    else if (nf==3 .and. all(f==(/'GEO ','GEO ','EL  '/))) then
!-temporarily disabled by Bin Gao before implementing operator for dipole length integrals
!-#ifdef BUILD_GEN1INT
!-!FIXME: not so sure if \var(tmpfgg) is correct!! needs to check!!
!-       ! calculates the second order total geometric derivatives of dipole length integrals
!-       tmpfgg = 0.0
!-       call gen1int_ifc_main("DIPLEN", .false., 0, 3, 2,             &
!-                             .false., .false., .false.,              &
!-                             size(tmpfgg), tmpfgg, 0, D%nrow*D%ncol, &
!-                             1, D%elms, .true., .true., .false.,     &
!-                             tmpfgg, get_lupri(), 5)
!-       tmpfgg = 2.0*tmpfgg
!-       indx = 0
!-       do i = 1, nc(3)
!-         do j = 1, nc(2)
!-           do k = 1, nc(1)
!-             indx = indx+1
!-             ave(indx) = tmpfgg(i,k,j)
!-           end do
!-         end do
!-       end do
!-#else
       !dj use finite difference to calculate the second order total geometric
       !derivatives of dipole length integrals
       ncor = 3 * get_natom()
       allocate(fdigf(ncor,3))
       allocate(tmpggf(ncor,ncor,3))
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
       deallocate(fdigf)
       deallocate(tmpggf)
!-#endif
    else if (nf==3 .and. all(f==(/'GEO ','GEO ','GEO '/))) then
       ncor = 3 * get_natom()
       allocate(tmpggg(ncor,ncor,ncor))
       allocate(fdigg(ncor,ncor))
       do i = 0, nc(3)-1
          call SHELLS_NUCLEI_displace(c(3)+i, fdistep)
          call ONEDRV_ave_ifc(mol, f(1:2), size(fdigg), fdigg, D=D)
          tmpggg(:,:,c(3)+i) = fdigg / (2*fdistep)
          call SHELLS_NUCLEI_displace(c(3)+i, -2*fdistep)
          call ONEDRV_ave_ifc(mol, f(1:2), size(fdigg), fdigg, D=D)
          tmpggg(:,:,c(3)+i) = tmpggg(:,:,c(3)+i) - fdigg / (2*fdistep)
          call SHELLS_NUCLEI_displace(c(3)+i, fdistep)
       end do
       deallocate(fdigg)
       ave = reshape(tmpggg(c(1):c(1)+nc(1)-1, &
                            c(2):c(2)+nc(2)-1, &
                            c(3):c(3)+nc(3)-1), shape(ave))
       deallocate(tmpggg)
    else
       print *, 'rsp_oneave error: not implented or in wrong order - ', &
                (' ' // f(i), i=1,nf)
       call quit('rsp_oneave error: not implented or in wrong order')
    end if
#endif
  end subroutine



  !> Average 2-electron integrals perturbed by fields f over the
  !> product of density matrces D1 and D2
  subroutine rsp_twoave(mol, nf, f, c, nc, D1, D2, ave)
    use dalton_ifc,       only: dal_work
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
    if (nf==0) then
       ! contract second density to Fock matrix, then trace with first
       A(1) = mol%zeromat
       call mat_alloc(A(1))
       call di_get_gmat(D2, A(1)) !Coulomb and exchange
       ave(1) = tr(A(1),D1)
    else if (nf==1 .and. f(1)=='GEO ') then
       ncor = 3 * get_natom()
       allocate(tmp(ncor,1,1,1))
       n = mol%zeromat%nrow
       dal_work(     :n*n)   = reshape(D1%elms,(/n*n/))
       dal_work(n*n+1:n*n*2) = reshape(D2%elms,(/n*n/))
       call GRCONT(dal_work(n*n*2+1:), size(dal_work)-n*n*2, &
                   tmp(:,1,1,1), ncor, .true., .false., &
                   1, 0, .true., .false., dal_work(:n*n*2), 2)
       ave(:nc(1)) = tmp(c(1):c(1)+nc(1)-1,1,1,1)
       deallocate(tmp)
    else if (nf==2 .and. all(f==(/'GEO ','GEO '/))) then
       ncor = 3 * get_natom()
       allocate(tmp(ncor,ncor,1,1))
       n = mol%zeromat%nrow
       dal_work(     :n*n)   = reshape(D1%elms,(/n*n/))
       dal_work(n*n+1:n*n*2) = reshape(D2%elms,(/n*n/))
       call GRCONT(dal_work(n*n*2+1:), size(dal_work)-n*n*2, &
                   tmp(:,:,1,1), ncor**2, .true., .false., &
                   2, 0, .true., .false., dal_work(:n*n*2), 2)
       ave = reshape(tmp(c(1):c(1)+nc(1)-1, &
                         c(2):c(2)+nc(2)-1,1,1), shape(ave))
       deallocate(tmp)
    else if (nf==3 .and. all(f==(/'GEO ','GEO ','GEO '/))) then
       ! contract FULL cubic in tmp, unsymmetrized divided by six
       ncor = 3 * get_natom()
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
       ncor = 3 * get_natom()
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
  end subroutine

  !> Exchange-correlation perturbed by fields f, averaged over densities D
  subroutine rsp_xcave(geo_order, nr_dmat, D, res)

     use xcint_main
     
!    ---------------------------------------------------------------------------
     integer,      intent(in)  :: geo_order
     integer,      intent(in)  :: nr_dmat
     type(matrix), intent(in)  :: D(*)
     complex(8),   intent(out) :: res(*)
!    ---------------------------------------------------------------------------
     integer                   :: i
     integer                   :: imat
     integer                   :: mat_dim
     integer                   :: nr_atoms
     real(8),      allocatable :: xc_dmat(:)
     real(8),      allocatable :: res_real(:)
!    ---------------------------------------------------------------------------
    
     nr_atoms = get_natom()
 
     if (.not. is_ks_calculation()) then
        res(1:(nr_atoms*3)**geo_order) = 0.0d0
        return
     end if
     
     allocate(res_real((nr_atoms*3)**geo_order))
     res_real = 0.0d0
     
     mat_dim = D(1)%nrow
     allocate(xc_dmat(mat_dim*mat_dim*nr_dmat))
     xc_dmat = 0.0d0
     do imat = 1, nr_dmat
        call daxpy(mat_dim*mat_dim, 1.0d0, D(imat)%elms, 1, xc_dmat((imat-1)*mat_dim*mat_dim + 1), 1)
     end do
     call xc_integrate(                        &
                       xc_mat_dim=mat_dim,     &
                       xc_dmat=xc_dmat,        &
                       xc_res=res_real,        &
                       xc_nr_dmat=nr_dmat,     &
                       xc_get_ave=.true.,      &
                       xc_geo_order=geo_order, &
                       xc_nr_atoms=nr_atoms    &
                      )
     deallocate(xc_dmat)
     
     do i = 1, size(res_real)
        res(i) = cmplx(res_real(i), 0.0d0)
     end do
     deallocate(res_real)

  end subroutine

  !> Exchange-correlation perturbed by fields f, averaged over densities D
  subroutine rsp_xcave_new(geo_order, nr_dmat, D, res)

     use xcint_main
     
!    ---------------------------------------------------------------------------
     integer,      intent(in)  :: geo_order
     integer,      intent(in)  :: nr_dmat
     type(matrix), intent(in)  :: D(*)
     complex(8),   intent(out) :: res(*)
!    ---------------------------------------------------------------------------
     integer                   :: i, j, k, l
     integer                   :: imat
     integer                   :: mat_dim
     integer                   :: nr_atoms
     real(8)                   :: res_real(1) !fixme
!    ---------------------------------------------------------------------------
    
     nr_atoms = get_natom()
 
     if (.not. is_ks_calculation()) then
        res(1:(nr_atoms*3)**geo_order) = 0.0d0
        return
     end if
    
     mat_dim = D(1)%nrow

     select case (geo_order)
        case (1)
           do i = 1, nr_atoms*3
              call xc_integrate(                       &
                                xc_mat_dim=mat_dim,    &
                                xc_nr_dmat=1,          &
                                xc_dmat=(/D(1)%elms/), &
                                xc_res=res_real,       &
                                xc_get_ave=.true.,     &
                                xc_geo_coor=(/i/)      &
                               )
              res(i) = cmplx(res_real(1), 0.0d0)
           end do
        case default
           print *, 'error: order too hight in xcave_new'
           stop 1
     end select

  end subroutine


  !> Compute differentiated overlap matrices, and optionally
  !> add half-differentiated overlap contribution to Fock matrices
  subroutine rsp_ovlint(mol, nf, f, c, nc, ovl, w, fock)
#ifdef BUILD_GEN1INT
    ! Gen1Int interface in Dalton
    use gen1int_interface
#endif
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
#ifdef BUILD_GEN1INT
    ! origins in Dalton
#include "orgcom.h"
    integer order_geo         !order of geometric derivatives
    integer num_coord         !number of atomic coordinates
    integer dim_redunt_geo    !size of redundant total geometric derivatives
    type(one_prop_t) overlap  !operator of overlap integrals
    integer ierr              !error information
#endif
    if (present(w) .and. .not.present(fock)) &
       call quit("error in rsp_ovlint: frequencies 'w' and Fock matrix 'fock' " &
              // 'must both be present or both absent')
#ifdef BUILD_GEN1INT
    ! gets the order of geometric derivatives
    order_geo = count(f=='GEO ')
    if (order_geo/=nf) &
      call quit("rsp_ovlint>> only geometric derivatives implemented!")
    ! sets the number of atomic coordinates
    num_coord = 3*get_natom()
    ! sets the size of redundant total geometric derivatives
    dim_redunt_geo = num_coord**order_geo
    if (dim_redunt_geo/=size(ovl)) &
      call quit("rsp_ovlint>> returning specific components not implemented!")
    ! creates operator for overlap integrals
    call OnePropCreate(prop_name=INT_OVERLAP, one_prop=overlap, &
                       info_prop=ierr, dipole_origin=DIPORG)
    if (ierr/=0) &
      call quit("rsp_ovlint>> failed to create operator of overlap integrals!")
!FIXME changes to call Gen1Int
    ! with field frequencies
    if (present(w)) then
      if (order_geo==1) then
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
!FIXME to Andreas: do we need higher order geometric derivatives of overlap integrals with frequencies?
      else
        call quit('rsp_ovlint>> GEO(>1) with freqencies not implemented!')
      end if
    else
      ! allocates matrices
      do i = 1, dim_redunt_geo
        if (.not.isdef(ovl(i))) then
          ovl(i) = mol%zeromat
          call mat_alloc(ovl(i))
        end if
      end do
      ! calculates the geometric derivatives of overlap matrix
!FIXME: rewrites \fn(gen1int_ifc_main) so that we could pass information of basis sets
      call gen1int_ifc_main(one_prop=overlap,                          &
                            order_geo_total=order_geo,                 &
                            max_num_cent=min(2,order_geo,get_natom()), &
                            num_ints=dim_redunt_geo,                   &
                            val_ints=ovl,                              &
                            redunt_ints=order_geo>1,                   &
                            num_dens=1,                                &
                            io_viewer=get_lupri(),                     &
                            level_print=5)
    end if
    ! frees space
    call OnePropDestroy(one_prop=overlap)
#else
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
#endif
  end subroutine



  subroutine rsp_oneint(mol, nf, f, c, nc, oneint)
#ifdef BUILD_GEN1INT
    ! Gen1Int interface in Dalton
    use gen1int_interface
#endif
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
#ifdef BUILD_GEN1INT
    ! origins in Dalton
#include "orgcom.h"
    ! uses \var(MXCENT)
#include "mxcent.h"
    ! coordinates and charges of atoms
#include "nuclei.h"
    integer order_geo           !order of geometric derivatives
    integer num_coord           !number of atomic coordinates
    integer dim_redunt_geo      !size of redundant total geometric derivatives
    type(one_prop_t) one_hamil  !operator of one-electron Hamiltonian
    integer ierr                !error information
!FIXME: implements dipole length, etc.
    ! gets the order of geometric derivatives
    order_geo = count(f=='GEO ')
    if (order_geo/=nf) &
      call quit("rsp_oneint>> only geometric derivatives implemented!")
    ! sets the number of atomic coordinates
    num_coord = 3*get_natom()
    ! sets the size of redundant total geometric derivatives
    dim_redunt_geo = num_coord**order_geo
    if (dim_redunt_geo/=size(oneint)) &
      call quit("rsp_oneint>> returning specific components not implemented!")
    ! creates operator of one-electron Hamiltonian
    call OnePropCreate(prop_name=INT_ONE_HAMIL, one_prop=one_hamil,   &
                       info_prop=ierr, coord_nuclei=CORD(:,1:NUCDEP), &
                       charge_nuclei=-CHARGE(1:NUCDEP))
    if (ierr/=0) &
      call quit("rsp_oneint>> failed to create operator of one-electron Hamiltonian!")
    ! allocates matrices
    do i = 1, dim_redunt_geo
      if (.not.isdef(oneint(i))) then
        oneint(i) = mol%zeromat
        call mat_alloc(oneint(i))
      end if
    end do
    ! calculates the geometric derivatives of one-electron Hamiltonian
!FIXME: rewrites \fn(gen1int_ifc_main) so that we could pass information of basis sets
    call gen1int_ifc_main(one_prop=one_hamil,                        &
                          order_geo_total=order_geo,                 &
                          max_num_cent=min(3,order_geo,get_natom()), &
                          num_ints=dim_redunt_geo,                   &
                          val_ints=oneint,                           &
                          redunt_ints=order_geo>1,                   &
                          num_dens=1,                                &
                          io_viewer=get_lupri(),                     &
                          level_print=5)
    ! frees space
    call OnePropDestroy(one_prop=one_hamil)
#else
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
#endif
  end subroutine



  !> Contract 2-electron integrals perturbed by fields 'f' with density
  !> matrix 'dens', and add to Fock matrices 'fock' Average 2-electron integrals perturbed by fields f over the
  !> product of density matrces D1 and D2
  subroutine rsp_twoint(mol, nf, f, c, nc, dens, fock)
    ! work array to be passed to GRCONT
    use dalton_ifc,       only: dal_work
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
             call GRCONT(dal_work(n*n*3+1:), size(dal_work)-n*n*3, &
                         dal_work(:n*n*3), n*n*3, .true., .false., &
                         1, (c(1)+i+2)/3, .false., .true., dens%elms, 1)
          j = 1 + mod(c(1)+i-1,3) !x y z = 1 2 3
          if (iszero(fock(1+i))) then
             call mat_alloc(fock(1+i))
             fock(1+i)%elms = reshape(dal_work(n*n*(j-1)+1:n*n*j),(/n,n/))
          else
             fock(1+i)%elms = fock(1+i)%elms &
                            + reshape(dal_work(n*n*(j-1)+1:n*n*j),(/n,n/))
          end if
       end do
    else if (nf==2 .and. all(f==(/'GEO ','GEO '/))) then
       ncor = 3 * get_natom()
       nullify(null_ptr) !because null() isn't f90
       do j = 0, nc(2)-1
          do i = 0, nc(1)-1
             ij = 1 + i + nc(1)*j
             if (iszero(fock(ij))) then
                call mat_alloc(fock(ij))
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
  end subroutine



  !> Exchange-correlation perturbed by fields 'field', contracted over
  !> densities 'D', added to Fock matrices 'F'
  subroutine rsp_xcint(geo_order, nr_dmat, D, res)

     use xcint_main
     
!    ---------------------------------------------------------------------------
     integer,      intent(in)    :: geo_order
     integer,      intent(in)    :: nr_dmat
     type(matrix), intent(in)    :: D(*)
     type(matrix), intent(inout) :: res(*)
!    ---------------------------------------------------------------------------
     integer                     :: icenter
     integer                     :: ixyz
     integer                     :: ioff
     integer                     :: imat
     integer                     :: mat_dim
     integer                     :: nr_atoms
     integer                     :: nr_fmat
     real(8),      allocatable   :: xc_dmat(:)
     real(8),      allocatable   :: xc_res(:)
!    ---------------------------------------------------------------------------
     
     if (.not. is_ks_calculation()) then
        return
     end if
     
     mat_dim = D(1)%nrow
     allocate(xc_dmat(mat_dim*mat_dim*nr_dmat))
     xc_dmat = 0.0d0
     do imat = 1, nr_dmat
        call daxpy(mat_dim*mat_dim, 1.0d0, D(imat)%elms, 1, xc_dmat((imat-1)*mat_dim*mat_dim + 1), 1)
     end do

     nr_fmat = 3**geo_order
     allocate(xc_res(mat_dim*mat_dim*nr_fmat))

     select case (geo_order)
        case (0)
           call xc_integrate(                       &
                             xc_mat_dim=mat_dim,    &
                             xc_dmat=xc_dmat,       &
                             xc_res=xc_res,         &
                             xc_nr_dmat=nr_dmat,    &
                             xc_nr_fmat=nr_fmat,    &
                             xc_geo_order=geo_order &
                            )
           
           do imat = 1, nr_fmat
              call daxpy(mat_dim*mat_dim, 1.0d0, xc_res((imat-1)*mat_dim*mat_dim + 1), 1, res(imat)%elms, 1)
           end do
        case (1)
           do icenter = 1, get_natom()
              call xc_integrate(                         &
                                xc_mat_dim=mat_dim,      &
                                xc_dmat=xc_dmat,         &
                                xc_res=xc_res,           &
                                xc_nr_dmat=nr_dmat,      &
                                xc_nr_fmat=nr_fmat,      &
                                xc_geo_order=geo_order,  &
                                xc_nr_atoms=get_natom(), &
                                xc_cent=icenter          &
                               )
           
              ioff = (icenter-1)*3
              do imat = 1, nr_fmat
                 call daxpy(mat_dim*mat_dim, 1.0d0, xc_res((imat-1)*mat_dim*mat_dim + 1), 1, res(ioff+imat)%elms, 1)
              end do
           end do
        case default
           print *, 'error: order not implemented in rsp_xcint'
           stop 1
     end select
     
     deallocate(xc_dmat)
     deallocate(xc_res)

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
          rsp_field_dim(i) = 3 * get_natom()
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
  subroutine ONEDRV_ave_ifc(mol, fld, siz, ave, D, DFD)
    use dalton_ifc, only: dal_work
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
    nc = 3 * get_natom()
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
    use dalton_ifc, only: dal_work
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
  
  
  function rank_one_pointer(siz, arr) result(ptr)
     integer,         intent(in) :: siz
     real(8), target, intent(in) :: arr(siz)
     real(8), pointer            :: ptr(:)
     ptr => arr
  end function


#ifdef BUILD_GEN1INT
  !> gets required expectation values of geometric derivatives
  subroutine rsp_set_geom(mol, nf, c, nc, val_expt, ave)
    !> structure containing integral program settings
    type(rsp_cfg), intent(in)  :: mol
    !> number of fields
    integer,       intent(in)  :: nf
    !> first and number of- components in each field
    integer,       intent(in)  :: c(nf), nc(nf)
    !> all expectaion values
    real(8),       intent(in)  :: val_expt(:)
    !> output average
    complex(8),    intent(out) :: ave(product(nc))
    !----------------------------------------------
    integer num_coord                    !number of atomic coordinates
    integer, allocatable :: powers_nc(:) !powers of \var(num_coord)
    integer, allocatable :: max_geo(:)   !maximum indices of output average
    integer, allocatable :: idx_geo(:)   !indices of a specific expectaion value
    integer addr_geo                     !address of a specific expectation value
    logical do_geo                       !indicates assigning output average
    integer ierr                         !error information
    integer ixyz, igeo                   !incremental recorders
    ! sets the number of atomic coordinates
    num_coord = 3*get_natom()
    ! sets the powers of \var(num_coord)
    allocate(powers_nc(nf), stat=ierr)
    if (ierr/=0) call quit("rsp_set_geom>> failed to allocate powers_nc!")
    powers_nc(1) = 1
    do ixyz = 2, nf
      powers_nc(ixyz) = num_coord*powers_nc(ixyz-1)
    end do
    ! maximum indices for output average
    allocate(max_geo(nf), stat=ierr)
    if (ierr/=0) call quit("rsp_set_geom>> failed to allocate max_geo!")
    do ixyz = 1, nf
      max_geo(ixyz) = c(ixyz)+nc(ixyz)-1
    end do
    ! indices of a specific expectaion value
    allocate(idx_geo(nf), stat=ierr)
    if (ierr/=0) call quit("rsp_set_geom>> failed to allocate idx_geo!")
    idx_geo = c
    do_geo = .true.
    igeo = 0
    do while (do_geo)
      ! gets the address of returned expectaion value
      addr_geo = idx_geo(1)
      do ixyz = 2, nf
        addr_geo = addr_geo+(idx_geo(ixyz)-1)*powers_nc(ixyz)
      end do
      ! assigns the output average
      igeo = igeo+1
      ave(igeo) = val_expt(addr_geo)
      ! generates new indices of the specific expectation value
      do_geo = .false.
      do ixyz = 1, nf-1
        if (idx_geo(ixyz)<max_geo(ixyz)) then
          idx_geo(ixyz) = idx_geo(ixyz)+1
          do_geo = .true.
          exit
        else
          idx_geo(ixyz) = c(ixyz)
        end if
      end do
      if ((.not.do_geo) .and. (idx_geo(nf)<max_geo(nf))) then
        idx_geo(nf) = idx_geo(nf)+1
        do_geo = .true.
      end if
    end do
    ! frees space
    deallocate(powers_nc)
    deallocate(max_geo)
    deallocate(idx_geo)
  end subroutine


#endif


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


  !> Add the average contribution 'ave', with 'canonical' (integral program-)
  !> ordering of dimensions, and full/total component ranges 'tc', to response
  !> tensor 'rsp', in which the dimensions are ordered arbitrarily, related to
  !> those in 'ave' by ordering 'perm', and with component ranges
  !> 'comp : comp + nc-1'. 'dcomp', either zero or a number within each range,
  !> specify the indices belong to the density matrix (!=0), whose dimensions
  !> are not present in ave. After each call to an integral program, this subroutine
  !> will be used to place the contribution in the resulting response tensor.
  !> 'neg' and 'imag' determine the phase factor of addition.
  !> @neg   minus sign, rsp = rsp - ave
  !> @imag  imaginary,  rsp = rsp + i*ave
  !> @ndim  number of dimensions
  !> @perm  ordering of dimensions. For each dim in rsp, the corresponding dim in ave
  !> @dcomp indices of density dimensions or zero for integral dimension
  !> @comp  startint component for each field in rsp
  !> @tc    total number of components. Dimensions of ave
  !> @nc    number of components. Dimensions of rsp
  subroutine perm_comp_add(neg, imag, ndim, perm, dcomp, comp, tc, nc, ave, rsp)
    logical,    intent(in)    :: neg, imag
    integer,    intent(in)    :: ndim, perm(ndim), dcomp(ndim)
    integer,    intent(in)    :: comp(ndim), nc(ndim), tc(ndim)
    real(8),    intent(in)    :: ave(product(tc)) !tc = total num comp
    complex(8), intent(inout) :: rsp(product(nc)) !nc = num comp
    
  end subroutine

end module
