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

  use matrix_backend
  use dalton_ifc
  use rsp_functions
  use rsp_backend
  use rsp_contribs, only: rsp_cfg

#ifndef OPENRSP_STANDALONE
! xcint
  use num_grid
  use interface_ao_specific
  use xcint_main
#endif /* OPENRSP_STANDALONE */

  implicit none

  public openrsp_setup
  public openrsp_calc
  public openrsp_finalize
  private

  ! ------------ MO response solver settings -----------
  ! maximum number of micro iterations in the iterative solution of
  ! the frequency independent linear response functions
  integer :: solver_maxit = 100
  ! maximum dimension of the sub-block of the configuration Hessian
  integer :: solver_maxphp = 0
  ! maximum dimension of the reduced space to which new basis vectors are added
  integer :: solver_mxrm = 400
  ! convergence threshold for the solution of the frequency-independent response equations
  real(8) :: solver_thresh = 1.0D-10
  ! true for optimal orbital trial vectors in the iterative solution of
  ! the frequency-dependent linear response equations
  logical :: solver_optorb = .false.
  
  ! ------------- SCF state and settings ------------
  ! config
  type(rsp_cfg) cfg
  ! overlap matrix
  type(matrix) S
  ! AO density matrix
  type(matrix) D
  ! Fock matrix
  type(matrix) F

  ! -------------- response function settings -----------
  ! number of frequencies
  integer num_freq
  ! (real parts of) frequencies
  real(8), allocatable :: real_freqs(:)
  ! imaginary frequencies
  real(8), allocatable :: imag_freqs(:)

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


  subroutine openrsp_setup(NBAST, WAVPCM, LUCMD, LUPRI, LWORK, WORK)
    use rsp_backend,  only: rsp_backend_setup, get_natom
    logical, intent(in)    :: WAVPCM
    integer, intent(in)    :: NBAST, LUCMD, LUPRI, LWORK
    real(8), intent(inout) :: WORK(LWORK)
    type(matrix)           :: H1 !one electron Hamiltonian
    integer                :: kfree
    type(matrix)           :: X(1) !XC contribution to fock matrix
    
    ! prints the header and license information
    call TITLER('OpenRSP: Response functions computed using AO basis', '*', -1)
    write (LUPRI,*) '>> -- solving equations with the MO response solver in DALTON.'
    write (LUPRI,*)
    write (LUPRI,*) '>> Original reference:'
    write (LUPRI,*) '>>  Andreas J. Thorvaldsen, Kenneth Ruud, Kasper Kristensen,'// &
                        ' Poul Jorgensen, and Sonia Coriani, '
    write (LUPRI,*) '>>  J. Chem. Phys. 129, 214108 (2008).'
    write (LUPRI,*)

    ! first parse input by running _calc in dryrun mode
    call openrsp_calc(LUCMD, .true.) !2nd arg: dryrun = .true.

    ! initialize the interface to DALTON
    call dal_ifc_init(WORK, LWORK, LUPRI, print_level, WAVPCM)

    ! initialize the MO response solver
    call rsp_mosolver_init(solver_maxit, solver_maxphp, solver_mxrm, &
                           solver_thresh, solver_optorb)

    ! initialize and allocate matrices
    call mat_nullify(S)
    S%nrow  = NBAST
    S%ncol  = S%nrow
    S%closed_shell = .true.
    S%magic_tag = mat_magic_setup
    call mat_alloc(S)
    call mat_nullify(D)
    call mat_setup(D, S)
    call mat_alloc(D)
    call mat_nullify(H1)
    call mat_setup(H1, S)
    call mat_alloc(H1)
    call mat_nullify(F)
    call mat_setup(F, S)
    call mat_alloc(F)

    ! get the overlap and one electron Hamiltonian matrices
    call di_get_overlap_and_H1(S, H1)
    ! get the AO density matrix, halving it
    call di_get_dens(D)
    call mat_axpy((1d0,0d0)/2, D, .false., .true., D)
    ! get the two electron contribution (G) to Fock matrix
    call di_get_gmat(D, F)
    ! Fock matrix F = H1 + G
    call mat_axpy((1d0,0d0), H1, .false., .false., F)
    call mat_free(H1)

    if (is_ks_calculation()) then
       ! write xcint interface files
       kfree = 1
       call write_num_grid_to_file(D%elms, work, kfree, lwork)
       call interface_ao_write()

#ifdef OPENRSP_STANDALONE
    print *, 'error: not part of standalone'
    stop 1
#else /* OPENRSP_STANDALONE */
       ! add xc contribution to the fock matrix
       X(1) = tiny(0.0d0)*F
       call xc_integrate(D=(/D/), F=X)
       F = F + X(1)
       X(1) = 0
#endif /* OPENRSP_STANDALONE */
    end if

    ! setup response config structure
    call mat_nullify(cfg%zeromat)
    call mat_setup(cfg%zeromat, S) !after setup, mat is zero to matrix_defop
    ! moved to rsp_backend: cfg%natom = num_atoms
    ! moved to rsp_backend: cfg%lupri = LUPRI
    ! moved to rsp_backend: cfg%hasxc = .false.
    ! moved to rsp_backend: call CHARGE_ifc(cfg%charge)
    ! moved to rsp_backend: call CORD_ifc(cfg%coord)
    ! setup basis descriptor cfg%basis
    call SHELLS_find_sizes(num_cgto_blocks, num_exp_and_ctr)
    allocate(cfg%basis(num_cgto_blocks))
    allocate(exp_and_ctr(num_exp_and_ctr))
    call SHELLS_to_type_cgto(num_cgto_blocks, num_exp_and_ctr, &
                             exp_and_ctr, cfg%basis)

    ! retrieve NATOMS from DALTON's common blocks
    num_atoms = get_natom()

end subroutine



  !> \brief If \param dryrun = .false., read commands from \param LUCMD
  !> and perform the requested calculations. If \param dryrun = .true.,
  !> only parse and validate the commands, but update the module's setup.
  !> \author Bin Gao
  !> \date 2009-12-08
  !> openrsp_setup must be called prior to this, and 
  !> \todo add some comments when calling Andreas' codes
  subroutine openrsp_calc(LUCMD, dryrun)

    integer, intent(in) :: LUCMD
    logical, intent(in) :: dryrun
    character(80) word
    integer       num_lines_read, err, l, i

    ! in case dryrun=T, we need to backspace LUCMD before returning,
    ! so keep track of number of lines read
    num_lines_read = 0

    do
       read (LUCMD, '(a)', end=888, err=999) word
       l = len_trim(word) !word/line's length
       num_lines_read = num_lines_read + 1

       ! Initial '!' or '#' are comments, skip those
       if (word(:1) == '!' .or. word(:1) == '#') cycle
       ! '*END OF' marks end of **AORSP section
       if (word(:7) == '*END OF') exit
       ! keywords start with '.', sections with '*'
       if (word(:1) /= '.' .and. word(:1) /= '*') &
          call quit("Failed to parse input. Expected word starting with" &
                    // " '.' or '*', got '" // word(:l) // "'")

       ! parse the keyword (or section)
       select case(word(:l))
       ! print level
       case ('.PRINT')
          read (LUCMD, *, end=888, err=999) print_level
          num_lines_read = num_lines_read + 1

       ! settings for MO response solver
       case ('.MAX IT')
          read (LUCMD, *, end=888, err=999) solver_maxit
          num_lines_read = num_lines_read + 1
       case ('.MAXPHP')
          read (LUCMD, *, end=888, err=999) solver_maxphp
          num_lines_read = num_lines_read + 1
       case ('.MAXRED')
          read (LUCMD, *, end=888, err=999) solver_mxrm
          num_lines_read = num_lines_read + 1
       case ('.THRESH')
          read (LUCMD, *, end=888, err=999) solver_thresh
          num_lines_read = num_lines_read + 1
       case ('.OPTORB')
          solver_optorb = .true.

       ! input (real part of) frequencies
       case ('.FREQ')
          ! read the number of frequencies first
          read (LUCMD, *, end=888, err=999) num_freq
          backspace(LUCMD)
          allocate(real_freqs(num_freq))
          read (LUCMD, *, end=888, err=999) num_freq, real_freqs
          num_lines_read = num_lines_read + 1

       ! input imaginary parts of frequencies
       case ('.IMFREQ')
          ! read the number of imaginary frequencies first
          read (LUCMD, *, end=888, err=999) num_freq
          backspace(LUCMD)
          allocate(imag_freqs(num_freq))
          read (LUCMD, *, end=888, err=999) num_freq, imag_freqs
          num_lines_read = num_lines_read + 1

       ! calculate properties (unless dryrun)
       !case ('.EFGB')
       !case ('.CME')
       !case ('.ROA')
       !case ('.CARS')
       !case ('.JONES')
       !case ('.POLARIZ')
       !case ('.HYPOLAR2')
       !case ('.HYPOLAR')
       !case ('.SECHYP3')
       !case ('.SECHYP')
       !case ('.SECHYP1')
       !case ('.VIBBETA')
       !case ('.VIBSHYP')
       case ('.GRADIENT')
          if (.not.dryrun) &
             call prop_test_gradient(cfg, 3*num_atoms, S, D, F)
       case ('.HESSIAN')
          if (.not.dryrun) &
             call prop_test_hessian(cfg, 3*num_atoms, S, D, F)
       case ('.DIPHES')
          if (.not.dryrun) &
             call prop_test_diphes(cfg, 3*num_atoms, S, D, F)
       case ('.CUBICFF')
          if (.not.dryrun) &
             call prop_test_cubicff(cfg, 3*num_atoms, S, D, F)
       case ('.QUARTICFF')
          if (.not.dryrun) &
             call prop_test_quarticff(cfg, 3*num_atoms, S, D, F)

       ! illegal keyword
       case default
          call quit("Error: Keyword '" // word(:l) // "' is not recognized in OpenRSP!")

       end select
    end do

    ! if dryrun, rewind to where we started
888 if (dryrun) then
       do i = 1, num_lines_read
          backspace(LUCMD)
       end do
    end if
    return
    
999 call quit('Failed to process input "' // word(:l) // '"!')

  end subroutine

  
  
  subroutine openrsp_finalize()
    use rsp_backend, only: rsp_backend_finalize, &
                           get_lupri
    integer lupri
    lupri = get_lupri()
    ! free
    call mat_free(S)
    call mat_free(D)
    call mat_free(F)
    call rsp_mosolver_free
    call dal_ifc_finalize
    call rsp_backend_finalize !deallocate internals
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
subroutine OPENRSP_DRIVER(WORK, LWORK, WAVPCM)
  use openrsp
  implicit none
  integer, intent(in)    :: LWORK
  real(8), intent(inout) :: WORK(LWORK)
  logical, intent(in)    :: WAVPCM
  ! uses LUPRI, LUCMD which are the pre-defined unit numbers
#include <priunit.h>
  ! uses NBAST, NNBAST, NNBASX, N2BASX
#include <inforb.h>
  ! first setup, then calc, then finalize
  call QENTER('OpenRSP ')
  call openrsp_setup(NBAST, WAVPCM, LUCMD, LUPRI, LWORK, WORK)
  call openrsp_calc(LUCMD, .false.) !2nd arg: dryrun = .false.
  call openrsp_finalize
  call QEXIT('OpenRSP ')
end subroutine
