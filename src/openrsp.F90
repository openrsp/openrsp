!!  openrsp: interface for solving response equations using atomic basis
!!  Copyright 2009 Bin Gao
!!
!!  This file is part of openrsp.
!!
!!  openrsp is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  openrsp is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with openrsp. If not, see <http://www.gnu.org/licenses/>.
!!
!!  openrsp is an interface for solving response equations using
!!  Andreas J. Thorvaldsen's codes, and the solver in DALTON
!!
!!  NOTICE: DALTON is distributed under its own License Agreement
!!          You need to first obtain such agreement and DALTON
!!          from http://www.kjemi.uio.no/software/dalton/dalton.html
!!
!!  NOTICE: Andreas J. Thorvaldsen's codes are also distributed under
!!          GNU Lesser General Public License, you may also first
!!          obtain his code
!!
!!  This file is the main module of openrsp, provides external user interfaces,
!!  also contains the interface of calling Andreas J. Thorvaldsen's codes
!!
!!  2009-12-08, Bin Gao:
!!  * first version

!> \brief main module of openrsp
!> \details contains the external user interfaces,
!>          and the interface of calling Andreas J. Thorvaldsen's codes
!> \author Bin Gao
!> \date 2009-12-08
module openrsp

  use interface_molecule
  use interface_io
  use interface_xc
  use interface_pcm
  use interface_scf
  use interface_f77_memory
  use interface_rsp_solver
  use interface_1el
  use interface_basis
  use interface_dirac_gen1int
  use interface_interest
  use rsp_functions
  use rsp_field_tuple, only: p_tuple, p_tuple_standardorder
  use rsp_general, only: rsp_prop
  use dalton_ifc
  use openrsp_cfg
  use matrix_defop
  use matrix_lowlevel,  only: mat_init
  use eri_contractions, only: ctr_arg
  use eri_basis_loops,  only: unopt_geodiff_loop
  use vib_prop_old, only: load_vib_modes
  use vib_pv_contribs

#ifndef PRG_DIRAC
! xcint
  use interface_ao_specific
  use xcint_main
#endif

  implicit none

  public openrsp_setup
  public openrsp_calc
  public openrsp_finalize

  private

  ! ------------- SCF state and settings ------------
  ! overlap matrix
  type(matrix) S
  ! AO density matrix
  type(matrix), target :: D
  ! Fock matrix
  type(matrix) F

  ! --------------- geometry and basis settings -----------
  ! number of atoms
  integer num_atoms
  ! number of cgto shell blocks
  integer num_cgto_blocks
  ! number of exponents and contraction coefficients
  integer num_exp_and_ctr
  ! array to hold exponents and contraction coefs
  real(8), allocatable :: exp_and_ctr(:)

  ! print level to be used here and in solver
  integer :: print_level = 0

contains

  subroutine openrsp_setup(WAVPCM, LWORK, WORK)
    logical, intent(in)    :: WAVPCM
    integer, intent(in)    :: LWORK
    integer                :: nbast, lupri
    real(8), intent(inout) :: WORK(LWORK)
    type(matrix)           :: H1 !one electron Hamiltonian
    type(matrix)           :: G  !two electron Hamiltonian
    real(8), allocatable   :: xc_dmat(:)
    real(8), allocatable   :: xc_fmat(:)
    integer                :: mat_dim
    real(8)                :: xc_energy
    real(8), target        :: temp(1)
    real(8)                :: ave(100)

#ifdef PRG_DIRAC
    type(matrix)           :: TX, TY, TZ, Dp(1)
#endif

    type(ctr_arg) :: arg(1)
#ifdef VAR_LSDALTON
    STOP 'openrsp_setup'
#else
    call interface_molecule_init()
    call interface_io_init()
    call interface_xc_init()
    call interface_pcm_init(wavpcm)
    call interface_scf_init()
    call interface_basis_init()
#endif

    nbast = get_nr_ao()
    lupri = get_print_unit()

    ! prints the header and license information
    call TITLER('OpenRSP: Response functions computed using AO basis', '*', -1)
    write (LUPRI,*) '>> -- solving equations with the MO response solver in DALTON.'
    write (LUPRI,*)
    write (LUPRI,*) '>> Original reference:'
    write (LUPRI,*) '>>  Andreas J. Thorvaldsen, Kenneth Ruud, Kasper Kristensen,'// &
                        ' Poul Jorgensen, and Sonia Coriani, '
    write (LUPRI,*) '>>  J. Chem. Phys. 129, 214108 (2008).'
    write (LUPRI,*)

    call interface_f77_memory_init(work_len=lwork, work=work)

    ! initialize the MO response solver
    call rsp_mosolver_init(openrsp_cfg_solver_maxitr, &
                           openrsp_cfg_solver_maxphp, &
                           openrsp_cfg_solver_maxred, &
                           openrsp_cfg_solver_thresh, &
                           openrsp_cfg_solver_optorb)

    ! initialize and allocate overlap matrix
    call mat_init(S, NBAST, NBAST, .false., .false., .false., .false., .false.)

    ! get the overlap
    call interface_scf_get_s(S)

    ! get the AO density matrix and divide by 2
    D = 0*S
    call mat_ensure_alloc(D, only_alloc=.true.)
    call interface_scf_get_d(D)
    D = 0.5d0*D

    ! get the one electron Hamiltonian
    H1 = 0*S
    call mat_ensure_alloc(H1, only_alloc=.true.)
    call interface_scf_get_h1(H1)

    ! get the two electron contribution (G) to Fock matrix
    G = 0*S
    call mat_ensure_alloc(G, only_alloc=.true.)
    call interface_scf_get_g(D, G)

    ! Fock matrix F = H1 + G
    F = H1 + G

#ifdef PRG_DIRAC
!   radovan: the contributions below are extremely useful
!            for debugging
!            after i get dirac fully interfaced i will remove/clean it up
    print *, 'nr of electrons   =', dot(D, S)
    print *, '1-el energy       =', dot(H1, D)
    print *, '2-el energy       =', 0.5d0*dot(G, D)
    print *, 'electronic energy =', dot(H1, D) + 0.5d0*dot(G, D)

    TX = 0*D
    call mat_ensure_alloc(TX, only_alloc=.true.)
    TY = 0*D
    call mat_ensure_alloc(TY, only_alloc=.true.)
    TZ = 0*D
    call mat_ensure_alloc(TZ, only_alloc=.true.)
    call get_1el_integrals(                                &
                           M=(/TX, TY, TZ/),               &
                           prop_name="INT_CART_MULTIPOLE", &
                           num_ints=3,                     &
                           order_mom=1,                    &
                           order_elec=0,                   &
                           order_geo_total=0,              &
                           max_num_cent=0,                 &
                           blocks=(/1, 1, 2, 2/),          &
                           print_unit=get_print_unit()     &
                          )

    print *, 'dipole z          =', dot(TZ, D)
!   call rsp_mosolver_exec((/TZ/), (/0.0d0/), Dp)
!   print *, 'polarizability zz =', -dot(TZ, Dp(1))

    TX = 0
    TY = 0
    TZ = 0
#endif /* ifdef PRG_DIRAC */

    H1 = 0
    G  = 0

#ifndef PRG_DIRAC
    if (get_is_ks_calculation()) then
       ! write xcint interface files
       call interface_ao_write()

       ! add xc contribution to the fock matrix
       mat_dim = D%nrow
       allocate(xc_dmat(mat_dim*mat_dim))
       xc_dmat = 0.0d0
       call daxpy(mat_dim*mat_dim, 1.0d0, D%elms, 1, xc_dmat, 1)
       allocate(xc_fmat(mat_dim*mat_dim))
       call xc_integrate(                  &
                         mat_dim=mat_dim,  &
                         nr_dmat=1,        &
                         dmat=xc_dmat,     &
                         energy=xc_energy, &
                         get_ave=.false.,  &
                         fmat=xc_fmat,     &
                         geo_coor=(/0/)    &
                        )
       call daxpy(mat_dim*mat_dim, 1.0d0, xc_fmat, 1, F%elms, 1)
       deallocate(xc_dmat)
       deallocate(xc_fmat)
    end if
#endif /* ifdef PRG_DIRAC */

    num_atoms = get_nr_atoms()

end subroutine

  subroutine openrsp_calc()

    integer       :: kn(2)
    integer       :: i
    type(p_tuple) :: perturbation_tuple
real(8), dimension(3) :: fld_dum

    integer       :: n_nm, j, k, ierr
    real(8), allocatable, dimension(:) :: nm_freq, nm_freq_b
    real(8), allocatable, dimension(:,:) :: T
    complex(8), allocatable, dimension(:,:) :: egf_cart, egf_nm, ff_pv
    complex(8), allocatable, dimension(:,:,:) :: egff_cart, egff_nm, fff_pv
    complex(8), allocatable, dimension(:,:,:,:) :: egfff_cart, egfff_nm, ffff_pv


    if (openrsp_cfg_gradient) then
       call prop_test_gradient(3*num_atoms, S, D, F)
    end if

    if (openrsp_cfg_hessian) then
       call prop_test_hessian(3*num_atoms, S, D, F)
    end if

    if (openrsp_cfg_dipole_gradient) then
       call prop_test_diphes(3*num_atoms, S, D, F, .true.)
    end if

    if (openrsp_cfg_dipole_hessian) then
       call prop_test_diphes(3*num_atoms, S, D, F, .false.)
    end if

    if (openrsp_cfg_cubic_ff) then
       call prop_test_cubicff(3*num_atoms, S, D, F)
    end if

    if (openrsp_cfg_quartic_ff) then
       call prop_test_quarticff(3*num_atoms, S, D, F)
    end if

    if (openrsp_cfg_pnc_gradient) then
       call prop_test_pnc_gradient(3*num_atoms, S, D, F)
    end if

    if (openrsp_cfg_general_f) then
       kn = (/0, 0/)

       perturbation_tuple%n_perturbations = 1
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3/)
       perturbation_tuple%plab = (/'EL  '/)
       perturbation_tuple%pid = (/1/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_g) then
       kn = (/0, 0/)

       perturbation_tuple%n_perturbations = 1
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms/)
       perturbation_tuple%plab = (/'GEO '/)
       perturbation_tuple%pid = (/1/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_ff) then
       kn = (/0, 1/)

       perturbation_tuple%n_perturbations = 2
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3, 3/)
       perturbation_tuple%plab = (/'EL  ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/ &
          (-1.0) * sum(openrsp_cfg_real_freqs), openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_gf) then
       kn = (/0, 1/)

       perturbation_tuple%n_perturbations = 2
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3/)
       perturbation_tuple%plab = (/'GEO ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_gg) then
       kn = (/0, 1/)

       perturbation_tuple%n_perturbations = 2
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms/)
       perturbation_tuple%plab = (/'GEO ', 'GEO '/)
       perturbation_tuple%pid = (/1, 2/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_fff) then
       kn = (/1, 1/)

       perturbation_tuple%n_perturbations = 3
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3, 3, 3/)
       perturbation_tuple%plab = (/'EL  ', 'EL  ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/ &
          (-1.0) * sum(openrsp_cfg_real_freqs), openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_gff) then
       kn = (/0, 2/)

       perturbation_tuple%n_perturbations = 3
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3, 3/)
       perturbation_tuple%plab = (/'GEO ', 'EL  ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, &
          (-1.0) * sum(openrsp_cfg_real_freqs), openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_ggf) then
       kn = (/1, 1/)

       perturbation_tuple%n_perturbations = 3
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3/)
       perturbation_tuple%plab = (/'GEO ', 'GEO ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_ggg) then
       kn = (/1, 1/)

       perturbation_tuple%n_perturbations = 3
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3*num_atoms/)
       perturbation_tuple%plab = (/'GEO ', 'GEO ', 'GEO '/)
       perturbation_tuple%pid = (/1, 2, 3/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0/)

       end if           

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_ffff) then
       kn = (/2, 1/)

       perturbation_tuple%n_perturbations = 4
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3, 3, 3, 3/)
       perturbation_tuple%plab = (/'EL  ', 'EL  ', 'EL  ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3, 4/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/ &
          (-1.0) * sum(openrsp_cfg_real_freqs), openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_gfff) then
       kn = (/0, 3/)

       perturbation_tuple%n_perturbations = 4
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3, 3, 3/)
       perturbation_tuple%plab = (/'GEO ', 'EL  ', 'EL  ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3, 4/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, &
          (-1.0) * sum(openrsp_cfg_real_freqs), openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_ggff) then
       kn = (/1, 2/)

       perturbation_tuple%n_perturbations = 4
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3, 3/)
       perturbation_tuple%plab = (/'GEO ', 'GEO ', 'EL  ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3, 4/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0, &
          (-1.0) * sum(openrsp_cfg_real_freqs), openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_gggf) then
       kn = (/2, 1/)

       perturbation_tuple%n_perturbations = 4
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3*num_atoms, 3/)
       perturbation_tuple%plab = (/'GEO ', 'GEO ', 'GEO ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3, 4/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_gggg) then
       kn = (/2, 1/)

       perturbation_tuple%n_perturbations = 4
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3*num_atoms, &
                                   3*num_atoms/)
       perturbation_tuple%plab = (/'GEO ', 'GEO ', 'GEO ', 'GEO '/)
       perturbation_tuple%pid = (/1, 2, 3, 4/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_fffff) then
       kn = (/2, 2/)

       perturbation_tuple%n_perturbations = 5
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3, 3, 3, 3, 3/)
       perturbation_tuple%plab = (/'EL  ', 'EL  ', 'EL  ', 'EL  ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3, 4, 5/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/ &
          (-1.0) * sum(openrsp_cfg_real_freqs), openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_gffff) then
       kn = (/0, 4/)

       perturbation_tuple%n_perturbations = 5
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3, 3, 3, 3/)
       perturbation_tuple%plab = (/'GEO ', 'EL  ', 'EL  ', 'EL  ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3, 4, 5/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, &
          (-1.0) * sum(openrsp_cfg_real_freqs), openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_ggfff) then
       kn = (/1, 3/)

       perturbation_tuple%n_perturbations = 5
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3, 3, 3/)
       perturbation_tuple%plab = (/'GEO ', 'GEO ', 'EL  ', 'EL  ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3, 4, 5/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0, &
          (-1.0) * sum(openrsp_cfg_real_freqs), openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_gggff) then
       kn = (/2, 2/)

       perturbation_tuple%n_perturbations = 5
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3*num_atoms, 3, 3/)
       perturbation_tuple%plab = (/'GEO ', 'GEO ', 'GEO ', 'EL  ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3, 4, 5/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, &
          (-1.0) * sum(openrsp_cfg_real_freqs), openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_ggggf) then
       kn = (/2, 2/)

       perturbation_tuple%n_perturbations = 5
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3*num_atoms, &
                                   3*num_atoms, 3/)
       perturbation_tuple%plab = (/'GEO ', 'GEO ', 'GEO ', 'GEO ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3, 4, 5/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_ggggg) then
       kn = (/2, 2/)

       perturbation_tuple%n_perturbations = 5
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3*num_atoms, &
                                   3*num_atoms, 3*num_atoms/)
       perturbation_tuple%plab = (/'GEO ', 'GEO ', 'GEO ', 'GEO ', 'GEO '/)
       perturbation_tuple%pid = (/1, 2, 3, 4, 5/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_ffffff) then
       kn = (/3, 2/)

       perturbation_tuple%n_perturbations = 6
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3, 3, 3, 3, 3, 3/)
       perturbation_tuple%plab = (/'EL  ', 'EL  ', 'EL  ', 'EL  ', 'EL  ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3, 4, 5, 6/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/ &
          (-1.0) * sum(openrsp_cfg_real_freqs), openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_gfffff) then
       kn = (/0, 5/)

       perturbation_tuple%n_perturbations = 6
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3, 3, 3, 3, 3/)
       perturbation_tuple%plab = (/'GEO ', 'EL  ', 'EL  ', 'EL  ', 'EL  ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3, 4, 5, 6/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, &
          (-1.0) * sum(openrsp_cfg_real_freqs), openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_ggffff) then
       kn = (/1, 4/)

       perturbation_tuple%n_perturbations = 6
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3, 3, 3, 3/)
       perturbation_tuple%plab = (/'GEO ', 'GEO ', 'EL  ', 'EL  ', 'EL  ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3, 4, 5, 6/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0, &
          (-1.0) * sum(openrsp_cfg_real_freqs), openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_gggfff) then
       kn = (/2, 3/)

       perturbation_tuple%n_perturbations = 6
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3*num_atoms, 3, 3, 3/)
       perturbation_tuple%plab = (/'GEO ', 'GEO ', 'GEO ', 'EL  ', 'EL  ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3, 4, 5, 6/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, &
          (-1.0) * sum(openrsp_cfg_real_freqs), openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_ggggff) then
       kn = (/3, 2/)

       perturbation_tuple%n_perturbations = 6
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3*num_atoms, &
                                   3*num_atoms, 3, 3/)
       perturbation_tuple%plab = (/'GEO ', 'GEO ', 'GEO ', 'GEO ', 'EL  ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3, 4, 5, 6/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, &
          (-1.0) * sum(openrsp_cfg_real_freqs), openrsp_cfg_real_freqs(:)/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_gggggf) then
       kn = (/3, 2/)

       perturbation_tuple%n_perturbations = 6
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3*num_atoms, &
                                   3*num_atoms, 3*num_atoms, 3/)
       perturbation_tuple%plab = (/'GEO ', 'GEO ', 'GEO ', 'GEO ', 'GEO ', 'EL  '/)
       perturbation_tuple%pid = (/1, 2, 3, 4, 5, 6/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if

    if (openrsp_cfg_general_gggggg) then

       kn = (/3, 2/)

       perturbation_tuple%n_perturbations = 6
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%pdim = (/3*num_atoms, 3*num_atoms, 3*num_atoms, &
                                   3*num_atoms, 3*num_atoms, 3*num_atoms/)
       perturbation_tuple%plab = (/'GEO ', 'GEO ', 'GEO ', 'GEO ', 'GEO ', 'GEO '/)
       perturbation_tuple%pid = (/1, 2, 3, 4, 5, 6/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       else

          perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       end if

       call rsp_prop(perturbation_tuple, kn, F, D, S)
    end if


    if (openrsp_cfg_general_specify) then

       kn = openrsp_cfg_specify_kn

       perturbation_tuple%n_perturbations = openrsp_cfg_specify_order
       allocate(perturbation_tuple%pdim(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%plab(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%pid(perturbation_tuple%n_perturbations))
       allocate(perturbation_tuple%freq(perturbation_tuple%n_perturbations))

       perturbation_tuple%plab = openrsp_cfg_specify_plab

       do i = 1, openrsp_cfg_specify_order

          if (openrsp_cfg_specify_plab(i) == 'GEO ') then

             perturbation_tuple%pdim(i) = 3*num_atoms

          elseif (openrsp_cfg_specify_plab(i) == 'EL  ') then

             perturbation_tuple%pdim(i) = 3

          elseif (openrsp_cfg_specify_plab(i) == 'MAG ') then

             perturbation_tuple%pdim(i) = 3

          else

             write(*,*) 'ERROR: Unrecognized field', openrsp_cfg_specify_plab(i)
             write(*,*) 'encountered when setting up customized response property calculation'

          end if

       end do
       
       perturbation_tuple%pid = (/(i, i = 1, openrsp_cfg_specify_order)/)

       if(allocated(openrsp_cfg_real_freqs)) then

          perturbation_tuple%freq = (/(i * 0.0d0 , i = 1, openrsp_cfg_specify_order - &
          openrsp_cfg_nr_real_freqs - 1), (-1.0) * sum(openrsp_cfg_real_freqs), &
          openrsp_cfg_real_freqs(:) /)

       else

          perturbation_tuple%freq = (/(i * 0.0d0 , i = 1, openrsp_cfg_specify_order)/)

       end if

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/(i, i = 1, openrsp_cfg_specify_order)/)


       call rsp_prop(perturbation_tuple, kn, F, D, S)

    end if

    if (openrsp_cfg_general_pv2f) then

       ! Get normal mode transformation matrix
       ! Get normal mode frequencies

fld_dum = 0.0

       allocate(T(3*num_atoms, 3*num_atoms))
       allocate(nm_freq_b(3*num_atoms))

nm_freq_b = 0.0

       call load_vib_modes(3*num_atoms, n_nm, nm_freq_b, T)

       allocate(nm_freq(n_nm))

nm_freq = 0.0

       nm_freq = nm_freq_b(1:n_nm)
       deallocate(nm_freq_b)

       allocate(ff_pv(3, 3))
       allocate(egf_cart(3*num_atoms, 3))
       allocate(egf_nm(3*n_nm, 3))

ff_pv = 0.0
egf_cart = 0.0
egf_nm = 0.0


       ! Calculate gradient of dipole moment

       kn = (/0,1/)

       perturbation_tuple%n_perturbations = 2
       allocate(perturbation_tuple%pdim(2))
       allocate(perturbation_tuple%plab(2))
       allocate(perturbation_tuple%pid(2))
       allocate(perturbation_tuple%freq(2))

       perturbation_tuple%plab = (/'GEO ', 'EL  '/)
       perturbation_tuple%pdim = (/3*num_atoms, 3/)
       perturbation_tuple%pid = (/1, 2/)
       perturbation_tuple%freq = (/0.0d0, 0.0d0/)

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/1, 2/)

       call rsp_prop(perturbation_tuple, kn, F, D, S)


       ! Read dipole moment gradient from file and transform to normal mode basis

       open(unit = 258, file='rsp_tensor', status='old', action='read', iostat=ierr)

       do i = 1, 3*num_atoms
          read(258,*) fld_dum
          egf_cart(i, :) = fld_dum
       end do

       close(258)

       egf_nm = trans_cartnc_1w1d(3*num_atoms, n_nm, egf_cart, T(:,1:n_nm))

       ! Calculate PV contribution to polarizability

write(*,*) 'egf nm', egf_nm
       ff_pv = alpha_pv(n_nm, nm_freq, (/ (-1.0)* openrsp_cfg_real_freqs(1), &
                        openrsp_cfg_real_freqs(1) /), dm_1d = egf_nm)


       deallocate(T)
       deallocate(ff_pv)
       deallocate(egf_cart)
       deallocate(egf_nm)

       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)

    end if

    if (openrsp_cfg_general_pv3f) then

fld_dum = 0.0

       ! Get normal mode transformation matrix
       ! Get normal mode frequencies

       allocate(T(3*num_atoms, 3*num_atoms))

       allocate(nm_freq_b(3*num_atoms))

nm_freq_b = 0.0

       call load_vib_modes(3*num_atoms, n_nm, nm_freq_b, T)

       allocate(nm_freq(n_nm))

nm_freq = 0.0

       nm_freq = nm_freq_b(1:n_nm)
       deallocate(nm_freq_b)

       allocate(ff_pv(3, 3))
       allocate(fff_pv(3, 3, 3))
       allocate(egf_cart(3*num_atoms, 3))
       allocate(egf_nm(3*n_nm, 3))
       allocate(egff_cart(3*num_atoms, 3, 3))
       allocate(egff_nm(3*n_nm, 3, 3))

ff_pv = 0.0
fff_pv = 0.0
egf_cart = 0.0
egf_nm = 0.0
egff_cart = 0.0
egff_nm = 0.0


       ! Calculate gradient of dipole moment

       kn = (/0,1/)

       perturbation_tuple%n_perturbations = 2
       allocate(perturbation_tuple%pdim(2))
       allocate(perturbation_tuple%plab(2))
       allocate(perturbation_tuple%pid(2))
       allocate(perturbation_tuple%freq(2))

       perturbation_tuple%plab = (/'GEO ', 'EL  '/)
       perturbation_tuple%pdim = (/3*num_atoms, 3/)
       perturbation_tuple%pid = (/1, 2/)
       perturbation_tuple%freq = (/0.0d0, 0.0d0/)

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/1, 2/)

       call rsp_prop(perturbation_tuple, kn, F, D, S)


       ! Read dipole moment gradient from file and transform to normal mode basis

       open(unit = 258, file='rsp_tensor', status='old', action='read', iostat=ierr)

       do i = 1, 3*num_atoms
          read(258,*) fld_dum
          egf_cart(i, :) = fld_dum
       end do

       close(258)

       egf_nm = trans_cartnc_1w1d(3*num_atoms, n_nm, egf_cart, T(:,1:n_nm))

       ! Calculate PV contribution to polarizability

       ff_pv = alpha_pv(n_nm, nm_freq, (/ (-1.0)* openrsp_cfg_real_freqs(1), &
                        openrsp_cfg_real_freqs(1) /), dm_1d = egf_nm) 

       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)

      ! Calculate gradient of polarizability

       kn = (/0,2/)

       perturbation_tuple%n_perturbations = 3
       allocate(perturbation_tuple%pdim(3))
       allocate(perturbation_tuple%plab(3))
       allocate(perturbation_tuple%pid(3))
       allocate(perturbation_tuple%freq(3))

       perturbation_tuple%plab = (/'GEO ', 'EL  ', 'EL  '/)
       perturbation_tuple%pdim = (/3*num_atoms, 3, 3/)
       perturbation_tuple%pid = (/1, 2, 3/)
       perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0/)

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/1, 2, 3/)

       call rsp_prop(perturbation_tuple, kn, F, D, S)

       ! Read polarizability gradient from file and transform to normal mode basis

       open(unit = 258, file='rsp_tensor', status='old', action='read', iostat=ierr)

       do i = 1, 3*num_atoms
          do j = 1, 3
             read(258,*) fld_dum
             egff_cart(i, j, :) = fld_dum
          end do
       end do

       close(258)

       egff_nm = trans_cartnc_2w1d(3*num_atoms, n_nm, egff_cart, T(:,1:n_nm))

       ! Calculate PV contribution to 1st hyperpolarizability

       fff_pv = beta_pv(n_nm, nm_freq, (/ (-1.0)* sum(openrsp_cfg_real_freqs), &
                openrsp_cfg_real_freqs(1), openrsp_cfg_real_freqs(2) /), dm_1d = egf_nm, &
                po_1d = egff_nm)

       deallocate(T)
       deallocate(nm_freq)
       deallocate(ff_pv)
       deallocate(fff_pv)
       deallocate(egf_cart)
       deallocate(egf_nm)
       deallocate(egff_cart)
       deallocate(egff_nm)

       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)

    end if

    if (openrsp_cfg_general_pv4f) then

fld_dum = 0.0

       ! Get normal mode transformation matrix
       ! Get normal mode frequencies

       allocate(T(3*num_atoms, 3*num_atoms))

       allocate(nm_freq_b(3*num_atoms))

nm_freq_b = 0.0

       call load_vib_modes(3*num_atoms, n_nm, nm_freq_b, T)

       allocate(nm_freq(n_nm))

nm_freq = 0.0

       nm_freq = nm_freq_b(1:n_nm)
       deallocate(nm_freq_b)

       allocate(ff_pv(3, 3))
       allocate(fff_pv(3, 3, 3))
       allocate(ffff_pv(3, 3, 3, 3))
       allocate(egf_cart(3*num_atoms, 3))
       allocate(egf_nm(3*n_nm, 3))
       allocate(egff_cart(3*num_atoms, 3, 3))
       allocate(egff_nm(3*n_nm, 3, 3))
       allocate(egfff_cart(3*num_atoms, 3, 3, 3))
       allocate(egfff_nm(3*n_nm, 3, 3, 3))

ff_pv = 0.0
fff_pv = 0.0
ffff_pv = 0.0
egf_cart = 0.0
egf_nm = 0.0
egff_cart = 0.0
egff_nm = 0.0
egfff_cart = 0.0
egfff_nm = 0.0

       ! Calculate gradient of dipole moment

       kn = (/0,1/)

       perturbation_tuple%n_perturbations = 2
       allocate(perturbation_tuple%pdim(2))
       allocate(perturbation_tuple%plab(2))
       allocate(perturbation_tuple%pid(2))
       allocate(perturbation_tuple%freq(2))

       perturbation_tuple%plab = (/'GEO ', 'EL  '/)
       perturbation_tuple%pdim = (/3*num_atoms, 3/)
       perturbation_tuple%pid = (/1, 2/)
       perturbation_tuple%freq = (/0.0d0, 0.0d0/)

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/1, 2/)

       call rsp_prop(perturbation_tuple, kn, F, D, S)


       ! Read dipole moment gradient from file and transform to normal mode basis

       open(unit = 258, file='rsp_tensor', status='old', action='read', iostat=ierr)

       do i = 1, 3*num_atoms
          read(258,*) fld_dum
          egf_cart(i, :) = fld_dum
       end do

       close(258)

       egf_nm = trans_cartnc_1w1d(3*num_atoms, n_nm, egf_cart, T(:,1:n_nm))

       ! Calculate PV contribution to polarizability

       ff_pv = alpha_pv(n_nm, nm_freq, (/ (-1.0)* openrsp_cfg_real_freqs(1), &
                        openrsp_cfg_real_freqs(1) /), dm_1d = egf_nm) 


       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)


       ! Calculate gradient of polarizability

       kn = (/0,2/)

       perturbation_tuple%n_perturbations = 3
       allocate(perturbation_tuple%pdim(3))
       allocate(perturbation_tuple%plab(3))
       allocate(perturbation_tuple%pid(3))
       allocate(perturbation_tuple%freq(3))

       perturbation_tuple%plab = (/'GEO ', 'EL  ', 'EL  '/)
       perturbation_tuple%pdim = (/3*num_atoms, 3, 3/)
       perturbation_tuple%pid = (/1, 2, 3/)
       perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0/)

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/1, 2, 3/)

       call rsp_prop(perturbation_tuple, kn, F, D, S)

       ! Read polarizability gradient from file and transform to normal mode basis

       open(unit = 258, file='rsp_tensor', status='old', action='read', iostat=ierr)

       do i = 1, 3*num_atoms
          do j = 1, 3
             read(258,*) fld_dum
             egff_cart(i, j, :) = fld_dum
          end do
       end do

       close(258)

       egff_nm = trans_cartnc_2w1d(3*num_atoms, n_nm, egff_cart, T(:,1:n_nm))

       ! Calculate PV contribution to 1st hyperpolarizability

       fff_pv = beta_pv(n_nm, nm_freq, (/ (-1.0)* sum(openrsp_cfg_real_freqs), &
                openrsp_cfg_real_freqs(1), openrsp_cfg_real_freqs(2) /), dm_1d = egf_nm, &
                po_1d = egff_nm)

       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)

       ! Calculate gradient of first hyperpolarizability

       kn = (/0,3/)

       perturbation_tuple%n_perturbations = 4
       allocate(perturbation_tuple%pdim(4))
       allocate(perturbation_tuple%plab(4))
       allocate(perturbation_tuple%pid(4))
       allocate(perturbation_tuple%freq(4))

       perturbation_tuple%plab = (/'GEO ', 'EL  ', 'EL  ', 'EL  '/)
       perturbation_tuple%pdim = (/3*num_atoms, 3, 3, 3/)
       perturbation_tuple%pid = (/1, 2, 3, 4/)
       perturbation_tuple%freq = (/0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       perturbation_tuple = p_tuple_standardorder(perturbation_tuple)
       perturbation_tuple%pid = (/1, 2, 3, 4/)

       call rsp_prop(perturbation_tuple, kn, F, D, S)

       ! Read 1st hyperpolarizability gradient from file and transform to normal mode basis

       open(unit = 258, file='rsp_tensor', status='old', action='read', iostat=ierr)

       do i = 1, 3*num_atoms
          do j = 1, 3
             do k = 1, 3
                read(258,*) fld_dum
                egfff_cart(i, j, k, :) = fld_dum
             end do
          end do
       end do

       close(258)

       egfff_nm = trans_cartnc_3w1d(3*num_atoms, n_nm, egfff_cart, T(:,1:n_nm))

       ! Calculate PV contribution to 2nd hyperpolarizability

       ffff_pv = gamma_pv(n_nm, nm_freq, (/ (-1.0)* sum(openrsp_cfg_real_freqs), &
                          openrsp_cfg_real_freqs(1), openrsp_cfg_real_freqs(2), &
                          openrsp_cfg_real_freqs(3)/), dm_1d = egf_nm, po_1d = egff_nm, &
                          hp_1d = egfff_nm) 

       deallocate(T)
       deallocate(nm_freq)
       deallocate(ff_pv)
       deallocate(fff_pv)
       deallocate(ffff_pv)
       deallocate(egf_cart)
       deallocate(egf_nm)
       deallocate(egff_cart)
       deallocate(egff_nm)
       deallocate(egfff_cart)
       deallocate(egfff_nm)

       deallocate(perturbation_tuple%pdim)
       deallocate(perturbation_tuple%plab)
       deallocate(perturbation_tuple%pid)
       deallocate(perturbation_tuple%freq)

    end if




  end subroutine

  subroutine openrsp_finalize()
    integer lupri
    lupri = get_print_unit()
    ! free

    S = 0
    D = 0
    F = 0

    call rsp_mosolver_finalize()
#ifdef USE_WAVPCM
    if (get_is_pcm_calculation()) then
       call pcm_finalize()
    end if
#endif
    call interface_f77_memory_finalize()
    ! stamp with date, time and hostname
    call TSTAMP(' ', lupri)
    write (lupri,*)
    write (lupri,*) '>> ** End of OpenRSP Section'
    write (lupri,*)
  end subroutine

end module

!> \brief driver used for DALTON
!> \author Bin Gao
!> \date 2009-12-08
!> \param WORK contains the work memory
!> \param LWORK is the size of the work memory
subroutine openrsp_driver(work, lwork, wavpcm)

   use openrsp

   implicit none

   integer, intent(in)    :: lwork
   real(8), intent(inout) :: work(lwork)
   logical, intent(in)    :: wavpcm

   call openrsp_setup(wavpcm, lwork, work)
   call openrsp_calc()
   call openrsp_finalize()

end subroutine
