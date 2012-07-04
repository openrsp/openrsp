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
  use rsp_functions
  use rsp_general, only: p_tuple, rsp_prop
  use dalton_ifc
  use openrsp_cfg
  use matrix_defop

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
  type(matrix) D
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
    type(matrix)           :: T  !temp
    real(8), allocatable   :: xc_dmat(:)
    real(8), allocatable   :: xc_fmat(:)
    integer                :: mat_dim
    integer                :: algebra
    real(8)                :: xc_energy

    call interface_molecule_init()
    call interface_io_init()
    call interface_xc_init()
    call interface_pcm_init(wavpcm)
    call interface_scf_init()
    call interface_basis_init()

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

    algebra = 1
#ifdef PRG_DIRAC
    algebra = 4 !fixme hardcoded
#endif

    ! initialize and allocate overlap matrix
    call mat_init(S, nrow=NBAST, ncol=NBAST, closed_shell=.true., algebra=algebra)

    ! get the overlap
    call interface_scf_get_s(S)

    ! get the AO density matrix and divide by 2
    D  = mat_alloc_like(S)
    call interface_scf_get_d(D)
    D = 0.5d0*D

    ! get the one electron Hamiltonian
    H1 = mat_alloc_like(S)
    call interface_scf_get_h1(H1)

    ! get the two electron contribution (G) to Fock matrix
    G  = mat_alloc_like(S)
    call interface_scf_get_g(D, G)

    ! Fock matrix F = H1 + G
    F = H1 + G

    print *, 'nr of electrons from dot(D, S) =', dot(D, S)
    print *, '1-el electronic energy from dot(H1, D) =', dot(H1, D)
    T = H1*D
    print *, '1-el electronic energy from tr(H1*D) =', tr(T)
    T = 0
    print *, '2-el electronic energy from 0.5d0*dot(G, D)=', 0.5d0*dot(G, D)
    print *, 'electronic energy from dot(H1, D) + 0.5d0*dot(G, D) =', dot(H1, D) + 0.5d0*dot(G, D)

    H1 = 0
    G  = 0

#ifdef PRG_DIRAC
!   stop here, nothing below can work on the dirac side
    stop 1
#endif

#ifndef PRG_DIRAC
    if (get_is_ks_calculation()) then
       ! write xcint interface files
       call interface_ao_write()

       ! add xc contribution to the fock matrix
       mat_dim = D%nrow
       allocate(xc_dmat(mat_dim*mat_dim))
       xc_dmat = 0.0d0
       call daxpy(mat_dim*mat_dim, 1.0d0, D%elms_alpha, 1, xc_dmat, 1)
       allocate(xc_fmat(mat_dim*mat_dim))
       call xc_integrate(                     &
                         xc_mat_dim=mat_dim,  &
                         xc_nr_dmat=1,        &
                         xc_dmat=xc_dmat,     &
                         xc_energy=xc_energy, &
                         xc_fmat=xc_fmat      &
                        )
       call daxpy(mat_dim*mat_dim, 1.0d0, xc_fmat, 1, F%elms_alpha, 1)
       deallocate(xc_dmat)
       deallocate(xc_fmat)
    end if
#endif

    num_atoms = get_nr_atoms()

end subroutine

  subroutine openrsp_calc()

!   logical, intent(in) :: dryrun
!   integer             :: LUCMD
!   character(80) word
!   integer       num_lines_read, err, l, i
    integer       :: kn(2)
    type(p_tuple) :: perturbation_tuple


    if (openrsp_cfg_gradient) then
       call prop_test_gradient(3*num_atoms, S, D, F)
    end if

    if (openrsp_cfg_hessian) then
       call prop_test_hessian(3*num_atoms, S, D, F)
    end if

    if (openrsp_cfg_dipole_hessian) then
       call prop_test_diphes(3*num_atoms, S, D, F)
    end if

    if (openrsp_cfg_cubic_ff) then
       call prop_test_cubicff(3*num_atoms, S, D, F)
    end if

    if (openrsp_cfg_quartic_ff) then
       call prop_test_quarticff(3*num_atoms, S, D, F)
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
       perturbation_tuple%freq = (/0.0/)

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
       perturbation_tuple%freq = (/0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0/)             

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

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
       perturbation_tuple%freq = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

       call rsp_prop(perturbation_tuple, kn, F, D, S)
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
