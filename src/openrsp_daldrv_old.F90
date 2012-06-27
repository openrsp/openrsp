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
!!  This file contains the driver used for DALTON.
!!
!!  2010-02-27, Bin Gao:
!!  * adds the functionality of dumping molecule information
!!
!!  2009-12-08, Bin Gao:
!!  * first version

  !> \brief driver used for DALTON
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param WORK contains the work memory
  !> \param LWORK is the size of the work memory
  subroutine openrsp_daldrv_old( WORK, LWORK, WAVPCM )
    ! matrix
    use matrix_defop
    ! interface of DALTON
    use interface_f77_memory
    use interface_io
    use interface_xc
    use interface_pcm
    use interface_scf
    use interface_rsp_solver
    use interface_1el
    ! main module of openrsp
    use openrsp_old

    use interface_molecule

    implicit none

    integer LWORK
    real(8) WORK( LWORK )
    logical WAVPCM

    integer lupri, lucmd, nbast

    ! control information of openrsp
    type(rspinfo_t) openrsp_info
    ! read in characters from input file
    character*(80) word
    ! print level
    integer :: level_print = 10
    ! number of frequencies
    integer num_freq
    ! real frequencies
    real(8), allocatable :: real_freqs(:)
    ! imaginary frequencies
    real(8), allocatable :: imag_freqs(:)
    ! if calculates electric-field-gradient-induced (Buckingham) birefringence (EFGB)
    logical :: openrsp_efgb = .false.
    ! if calculates London Cotton-mouton constant
    logical :: openrsp_cme = .false.
    ! if calculates Raman optical activity (ROA) properties
    logical :: openrsp_roa = .false.
    ! if calculates coherent anti-Stokes Raman Scattering (CARS)
    logical :: openrsp_cars = .false.
    ! if calculates magnetoelectric Jones spectroscopy
    logical :: openrsp_jones = .false.
    ! if calculates linear polarizability
    logical :: openrsp_polariz = .false.
    ! if calculates 1st hyperpolarizability using n+1 rule
    logical :: openrsp_hypolar2 = .false.
    ! if calculates 1st hyperpolarizability using 2n+1 rule
    logical :: openrsp_hypolar = .false.
    ! if calculates 2nd hyperpolarizability using n+1 rule
    logical :: openrsp_sechyp3 = .false.
    ! if calculates 2nd hyperpolarizability using 2n+1 (1+2+1) rule
    logical :: openrsp_sechyp = .false.
    ! if calculates 2nd hyperpolarizability using 2n+1 (2+1+1) rule
    logical :: openrsp_sechyp1 = .false.
    ! if calculates vibrational hyperpolarizability
    logical :: openrsp_vibbeta = .false.
    ! if calculates vibrational 2nd hyperpolarizability
    logical :: openrsp_vibshyp = .false.

    ! the followings are for MO response solver
    !
    ! maximum number of micro iterations in the iterative solution of
    ! the frequency independent linear response functions
    integer :: solver_maxit = 100
    ! maximum dimension of the sub-block of the configuration Hessian
    integer :: solver_maxphp = 0
    ! maximum dimension of the reduced space to which new basis vectors are added
    integer :: solver_mxrm = 400
    ! convergence threshold for the solution of the frequency-independent response equations
    real(8) :: solver_thresh = 1.0D-07
    ! true for optimal orbital trial vectors in the iterative solution of
    ! the frequency-dependent linear response equations
    logical :: solver_optorb = .false.

    ! overlap matrix
    type(matrix) S
    ! AO density matrix
    type(matrix) D
    ! Fock matrix
    type(matrix) F
    ! one electron Hamiltonian
    type(matrix) H1
    ! two electron Hamiltonian
    type(matrix) G

    ! number of atoms
    integer num_atoms
    ! name of atoms
    character*4, allocatable :: aname(:)
    ! charge of atoms
    real(8), allocatable :: acharge(:)
    ! coordinates of atoms
    real(8), allocatable :: acoord(:,:)


    ! some constants
    real(8), parameter :: one = 1.0D+00
    real(8), parameter :: half = 5.0D-01
    ! temporary stuff
    integer i, j, k
    ! error information
    integer ierr
    integer mat_dim
    real(8), allocatable :: xc_fmat(:)
    real(8), allocatable :: xc_dmat(:)
    real(8)              :: xc_energy
    real(8), allocatable :: eigval(:)

    call interface_molecule_init()
    call interface_io_init()
    call interface_xc_init()
    call interface_scf_init()
    call interface_pcm_init(wavpcm)

    nbast = get_nr_ao()
    lupri = get_print_unit()
    lucmd = get_input_unit()

    call QENTER( 'OpenRSP ' )

    ! prints the header and license information
    call TITLER( 'OpenRSP: solve response equation using AO basis', '*', -1 )
    write( LUPRI, 100 ) '-- with Andreas J. Thorvaldsen''s code '// &
                        'and the MO response solver in DALTON.'
    write( LUPRI, '()' )
    write( LUPRI, 100 ) 'Reference:'
    write( LUPRI, 100 ) ' Andreas J. Thorvaldsen, Kenneth Ruud, Kasper Kristensen,'// &
                        ' Poul Jorgensen, and Sonia Coriani, '
    write( LUPRI, 100 ) ' J. Chem. Phys. 129, 214108 (2008).'
    write( LUPRI, '()' )
    if( WAVPCM ) then
       write( LUPRI, 100 ) '* Wavelet-PCM non-nonequilibrium calculation'
       write( LUPRI, '()' )
    end if

    ! dumps the molecule information
    num_atoms = get_nr_atoms()
    allocate( aname(num_atoms), stat=ierr )
    if ( ierr /= 0 ) call QUIT( 'Failed to allocate aname!' )
    allocate( acharge(num_atoms), stat=ierr )
    if ( ierr /= 0 ) call QUIT( 'Failed to allocate acharge!' )
    allocate( acoord(3,num_atoms), stat=ierr )
    if ( ierr /= 0 ) call QUIT( 'Failed to allocate acoord!' )

    do i = 1, num_atoms
       aname(i)     = get_nuc_name(i)
       acharge(i)   = get_nuc_charge(i)
       acoord(1, i) = get_nuc_xyz(1, i)
       acoord(2, i) = get_nuc_xyz(2, i)
       acoord(3, i) = get_nuc_xyz(3, i)
    end do

    ! dumps geometry for PyVib2
    write( LUPRI, '(2X,A)' ) 'Cartesian Coordinates (a.u.)'
    write( LUPRI, '(2X,A)' ) '----------------------------'
    write( LUPRI, '()' )
    write( LUPRI, '(2X,A,I5)' ) 'Total number of coordinates:', 3*num_atoms
    do i = 1, num_atoms
      write( LUPRI, '(2X,A,3X," : ",3(I4,2X,A,F15.10))' ) &
        aname(i), 3*i-2, 'x', acoord(1,i), 3*i-1, 'y', acoord(2,i), 3*i, 'z', acoord(3,i)
    end do
    write( LUPRI, '()' )
    deallocate( aname )
    deallocate( acharge )
    deallocate( acoord )

    ! processes the input information
    read( LUCMD, 110, err=999, end=999 ) word
    do while( word(1:12) /= '$END RESPONS' )
      ! reads comments
      if ( word(1:1) == '!' .or. word(1:1) == '#' ) then
        read( LUCMD, 110, err=999, end=999 ) word
        cycle
      end if
      select case ( trim( word ) )
      ! reads print level
      case ( '.PRINT' )
        read( LUCMD, *, err=999, end=999 ) level_print
      ! reads control information of MO response solver
      case ( '.MAX IT' )
        read( LUCMD, *, err=999, end=999 ) solver_maxit
      case ( '.MAXPHP' )
        read( LUCMD, *, err=999, end=999 ) solver_maxphp
      case ('.MAXRED')
        read( LUCMD, *, err=999, end=999 ) solver_mxrm
      case ( '.THRESH')
        read( LUCMD, *, err=999, end=999 ) solver_thresh
      case ( '.OPTORB')
        solver_optorb = .true.
      ! reads real frequencies
      case ( '.FREQ' )
        ! reads the number of real frequencies
        read( LUCMD, * ) num_freq
        allocate( real_freqs( num_freq ), stat=ierr )
        if ( ierr /= 0 ) call QUIT( 'Failed to allocate real_freqs!' )
        backspace( LUCMD )
        read( LUCMD, *, err=999, end=999 ) num_freq, ( real_freqs(i), i = 1, num_freq )
      ! reads imaginary frequencies
      case ( '.IMFREQ' )
        ! reads the number of imaginary frequencies
        read( LUCMD, * ) num_freq
        allocate( imag_freqs( num_freq ), stat=ierr )
        if ( ierr /= 0 ) call QUIT( 'Failed to allocate imag_freqs!' )
        backspace( LUCMD )
        read( LUCMD, *, err=999, end=999 ) num_freq, ( imag_freqs(i), i = 1, num_freq )
      ! reads calcualted properties
      case ( '.EFGB' )
        openrsp_efgb = .true.
      case ( '.CME' )
        openrsp_cme = .true.
      case ( '.ROA' )
        openrsp_roa = .true.
      case ( '.CARS' )
        openrsp_cars = .true.
      case ( '.JONES' )
        openrsp_jones = .true.
      case ( '.POLARIZ' )
        openrsp_polariz = .true.
      case ( '.HYPOLAR2' )
        openrsp_hypolar2 = .true.
      case ( '.HYPOLAR' )
        openrsp_hypolar = .true.
      case ( '.SECHYP3' )
        openrsp_sechyp3 = .true.
      case ( '.SECHYP' )
        openrsp_sechyp = .true.
      case ( '.SECHYP1' )
        openrsp_sechyp1 = .true.
      case ( '.VIBBETA' )
        openrsp_vibbeta = .true.
      case ( '.VIBSHYP' )
        openrsp_vibshyp = .true.
      ! illegal keyword
      case default
        call QUIT( ' Keyword "'//trim(word)//'" is not recognized in OpenRSP!' )
      end select
      ! reads next line
      read( LUCMD, 110, err=999, end=999 ) word
    end do

    ! default is static properties
    if ( .not. allocated(real_freqs) ) then
      allocate( real_freqs( 1 ), stat=ierr )
      if ( ierr /= 0 ) call QUIT( 'Failed to allocate real_freqs!' )
      real_freqs = 0.0D+00
    end if
    if ( .not. allocated(imag_freqs) ) then
      allocate( imag_freqs( size(real_freqs) ), stat=ierr )
      if ( ierr /= 0 ) call QUIT( 'Failed to allocate imag_freqs!' )
      imag_freqs = 0.0D+00
    end if

    call interface_f77_memory_init(work_len=lwork, work=work)

    ! initializes the MO response solver
    call rsp_mosolver_init( solver_maxit, solver_maxphp, solver_mxrm, &
                            solver_thresh, solver_optorb )

    ! sets the control information of openrsp
    call openrsp_info_set( openrsp_info,       &
                           LUPRI, level_print, &
                           real_freqs,         &
                           imag_freqs,         &
                           openrsp_efgb,       &
                           openrsp_cme,        &
                           openrsp_roa,        &
                           openrsp_cars,       &
                           openrsp_jones,      &
                           openrsp_polariz,    &
                           openrsp_hypolar2,   &
                           openrsp_hypolar,    &
                           openrsp_sechyp3,    &
                           openrsp_sechyp,     &
                           openrsp_sechyp1,    &
                           openrsp_vibbeta,    &
                           openrsp_vibshyp )
    deallocate( real_freqs )
    deallocate( imag_freqs )

    ! dumps the control information of openrsp
    call openrsp_info_dump( openrsp_info, LUPRI )

    ! dumps the control information of MO response solver
    call rsp_mosolver_dump( LUPRI )

    ! initialize and allocate matrices
    call mat_init(S, nrow=NBAST, ncol=NBAST, closed_shell=.true.)

    D  = mat_alloc_like(S)
    H1 = mat_alloc_like(S)
    G  = mat_alloc_like(S)

    ! get the overlap and one electron Hamiltonian matrices
    call di_get_overlap_and_H1(S, H1)

    ! get the AO density matrix, halving it
    call di_get_dens(D)

    D = 0.5d0*D

    ! get the two electron contribution (G) to Fock matrix
    call di_get_gmat(D, G)

    ! Fock matrix F = H1 + G
    F = H1 + G

    H1 = 0
    G  = 0

    if (get_is_ks_calculation()) then
       ! write xcint interface files
       call interface_ao_write()

       ! add xc contribution to the fock matrix
       mat_dim = D%nrow
       allocate(xc_dmat(mat_dim*mat_dim))
       xc_dmat = 0.0d0
       call daxpy(mat_dim*mat_dim, 1.0d0, D%elms, 1, xc_dmat, 1)
       allocate(xc_fmat(mat_dim*mat_dim))
       call xc_integrate(                     &
                         xc_mat_dim=mat_dim,  &
                         xc_nr_dmat=1,        &
                         xc_dmat=xc_dmat,     &
                         xc_energy=xc_energy, &
                         xc_fmat=xc_fmat      &
                        )
       call daxpy(mat_dim*mat_dim, 1.0d0, xc_fmat, 1, F%elms, 1)
       deallocate(xc_dmat)
       deallocate(xc_fmat)
    end if

    ! performs the calculations
    call openrsp_prop_calc(S, D, F, openrsp_info )

    ! cleans
    S = 0
    D = 0
    F = 0
    call openrsp_info_clean(openrsp_info)
    call rsp_mosolver_free
#ifdef USE_WAVPCM
    if (get_is_pcm_calculation()) then
       call pcm_finalize()
    end if
#endif
    call interface_f77_memory_finalize()

    ! stamps the date, time and hostname
    call TSTAMP( ' ', LUPRI )
    write( LUPRI, '()' )
    write( LUPRI, 100 ) '** End of OpenRSP Section'
    write( LUPRI, '()' )

    call QEXIT( 'OpenRSP ' )
    return

999 call QUIT( 'Failed to process input after reading '//trim( word )//'!' )

100  format(1X,'>> ',A,I6)
110  format(80A)
1010  format(1X,'>> ',A,4E12.6)
1020  format(1X,'>> ',A,4F10.6)
  end subroutine
