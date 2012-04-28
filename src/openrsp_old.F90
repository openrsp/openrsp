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
module openrsp_old
  ! matrix
  use matrix_backend
  ! response-related testing routines and some calculations
  use prop_test_old
  ! calculation and outputting of vibrational (optical) properties
  use vib_prop_old
  ! calculation and outputting of optical birefringences
  use birefring_old
  use interface_host

! xcint
#ifndef OPENRSP_STANDALONE
  use interface_ao_specific
  use xcint_main
#endif /* OPENRSP_STANDALONE */

  implicit none

  !> control information of openrsp
  type, public :: rspinfo_t
    private
    !> IO unit of log file
    integer :: log_io = 6
    !> print level
    integer :: level_print = 10
    !> real frequencies
    real(8), allocatable :: real_freqs(:)
    !> imaginary frequencies
    real(8), allocatable :: imag_freqs(:)
    !> if calculates electric-field-gradient-induced (Buckingham) birefringence (EFGB)
    logical :: openrsp_efgb = .false.
    !> if calculates London Cotton-mouton constant
    logical :: openrsp_cme = .false.
    !> if calculates Raman optical activity (ROA) properties
    logical :: openrsp_roa = .false.
    !> if calculates coherent anti-Stokes Raman Scattering (CARS)
    logical :: openrsp_cars = .false.
    !> if calculates magnetoelectric Jones spectroscopy
    logical :: openrsp_jones = .false.
    !> if calculates linear polarizability
    logical :: openrsp_polariz = .false.
    !> if calculates 1st hyperpolarizability using n+1 rule
    logical :: openrsp_hypolar2 = .false.
    !> if calculates 1st hyperpolarizability using 2n+1 rule
    logical :: openrsp_hypolar  = .false.
    !> if calculates 2nd hyperpolarizability using n+1 rule
    logical :: openrsp_sechyp3  = .false.
    !> if calculates 2nd hyperpolarizability using 2n+1 (1+2+1) rule
    logical :: openrsp_sechyp = .false.
    !> if calculates 2nd hyperpolarizability using 2n+1 (2+1+1) rule
    logical :: openrsp_sechyp1 = .false.
    !> if calculates vibrational hyperpolarizability
    logical :: openrsp_vibbeta = .false.
    !> if calculates vibrational 2nd hyperpolarizability
    logical :: openrsp_vibshyp = .false.
  end type rspinfo_t

  public :: openrsp_info_set
  public :: openrsp_info_dump
  public :: openrsp_info_clean
  public :: openrsp_prop_calc

  contains

  !> \brief sets the control information of openrsp
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param log_io is the IO unit of log file
  !> \param level_print is the print level
  !> \param real_freqs contains the real frequencies
  !> \param imag_freqs contains the imaginary frequencies
  !> \param openrsp_efgb indicates if calculating electric-field-gradient-induced (Buckingham)
  !>        birefringence (EFGB)
  !> \param openrsp_cme indicates if calculating London Cotton-mouton constant
  !> \param openrsp_roa indicates if calculating Raman optical activity (ROA) properties
  !> \param openrsp_cars indicates if calculating coherent anti-Stokes Raman Scattering (CARS)
  !> \param openrsp_jones indicates if calculating magnetoelectric Jones spectroscopy
  !> \param openrsp_polariz indicates if calculating linear polarizability
  !> \param openrsp_hypolar2 indicates if calculating 1st hyperpolarizability using n+1 rule
  !> \param openrsp_hypolar indicates if calculating 1st hyperpolarizability using 2n+1 rule
  !> \param openrsp_sechyp3 indicates if calculating 2nd hyperpolarizability using n+1 rule
  !> \param openrsp_sechyp indicates if calculating 2nd hyperpolarizability using 2n+1 (1+2+1) rule
  !> \param openrsp_sechyp1 indicates if calculating 2nd hyperpolarizability using 2n+1 (2+1+1) rule
  !> \param openrsp_vibbeta indicates if calculating vibrational hyperpolarizability
  !> \param openrsp_vibshyp indicates if calculating vibrational 2nd hyperpolarizability
  !> \return this_info contains the control information
  subroutine openrsp_info_set( this_info,           &
                               log_io, level_print, &
                               real_freqs,          &
                               imag_freqs,          &
                               openrsp_efgb,        &
                               openrsp_cme,         &
                               openrsp_roa,         &
                               openrsp_cars,        &
                               openrsp_jones,       &
                               openrsp_polariz,     &
                               openrsp_hypolar2,    &
                               openrsp_hypolar,     &
                               openrsp_sechyp3,     &
                               openrsp_sechyp,      &
                               openrsp_sechyp1,     &
                               openrsp_vibbeta,     &
                               openrsp_vibshyp )
    type(rspinfo_t), intent(inout) :: this_info
    integer, optional, intent(in) :: log_io
    integer, optional, intent(in) :: level_print
    real(8), optional, intent(in) :: real_freqs(:)
    real(8), optional, intent(in) :: imag_freqs(:)
    logical, optional, intent(in) :: openrsp_efgb
    logical, optional, intent(in) :: openrsp_cme
    logical, optional, intent(in) :: openrsp_roa
    logical, optional, intent(in) :: openrsp_cars
    logical, optional, intent(in) :: openrsp_jones
    logical, optional, intent(in) :: openrsp_polariz
    logical, optional, intent(in) :: openrsp_hypolar2
    logical, optional, intent(in) :: openrsp_hypolar
    logical, optional, intent(in) :: openrsp_sechyp3
    logical, optional, intent(in) :: openrsp_sechyp
    logical, optional, intent(in) :: openrsp_sechyp1
    logical, optional, intent(in) :: openrsp_vibbeta
    logical, optional, intent(in) :: openrsp_vibshyp
    ! error information
    integer ierr
    ! sets the IO unit of log file
    if ( present( log_io ) ) this_info%log_io = log_io
    ! sets the print level
    if ( present( level_print ) ) this_info%level_print = level_print
    ! sets the real frequencies
    if ( present( real_freqs) ) then
      allocate( this_info%real_freqs( size( real_freqs ) ), stat=ierr )
      if ( ierr /= 0 ) then
         print *, 'Failed to allocate real_freqs!'
         stop 1
      end if
      this_info%real_freqs = real_freqs
    end if
    ! sets the imaginary frequencies
    if ( present( imag_freqs ) ) then
      allocate( this_info%imag_freqs( size( imag_freqs ) ), stat=ierr )
      if ( ierr /= 0 ) then
         print *, 'Failed to allocate imag_freqs!'
         stop 1
      end if
      this_info%imag_freqs = imag_freqs
    end if
    ! if calculates electric-field-gradient-induced (Buckingham) birefringence (EFGB)
    if ( present( openrsp_efgb ) ) this_info%openrsp_efgb = openrsp_efgb
    ! if calculates London Cotton-mouton constant
    if ( present( openrsp_cme ) ) this_info%openrsp_cme = openrsp_cme
    ! if calculates Raman optical activity (ROA) properties
    if ( present( openrsp_roa ) ) this_info%openrsp_roa = openrsp_roa
    ! if calculates coherent anti-Stokes Raman Scattering (CARS)
    if ( present( openrsp_cars ) ) this_info%openrsp_cars = openrsp_cars
    ! if calculates magnetoelectric Jones spectroscopy
    if ( present( openrsp_jones ) ) this_info%openrsp_jones = openrsp_jones
    ! if calculates linear polarizability
    if ( present( openrsp_polariz ) ) this_info%openrsp_polariz = openrsp_polariz
    ! if calculates 1st hyperpolarizability using n+1 rule
    if ( present( openrsp_hypolar2 ) ) this_info%openrsp_hypolar2 = openrsp_hypolar2
    ! if calculates 1st hyperpolarizability using 2n+1 rule
    if ( present( openrsp_hypolar ) ) this_info%openrsp_hypolar = openrsp_hypolar
    ! if calculates 2nd hyperpolarizability using n+1 rule
    if ( present( openrsp_sechyp3 ) ) this_info%openrsp_sechyp3 = openrsp_sechyp3
    ! if calculates 2nd hyperpolarizability using 2n+1 (1+2+1) rule
    if ( present( openrsp_sechyp ) ) this_info%openrsp_sechyp = openrsp_sechyp
    ! if calculates 2nd hyperpolarizability using 2n+1 (2+1+1) rule
    if ( present( openrsp_sechyp1 ) ) this_info%openrsp_sechyp1 = openrsp_sechyp1
    ! if calculates vibrational hyperpolarizability
    if ( present( openrsp_vibbeta ) ) this_info%openrsp_vibbeta = openrsp_vibbeta
    ! if calculates vibrational 2nd hyperpolarizability
    if ( present( openrsp_vibshyp ) ) this_info%openrsp_vibshyp = openrsp_vibshyp
  end subroutine openrsp_info_set

  !> \brief dumps the control information of openrsp
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param this_info contains the control information
  !> \param io_dump is the IO unit to dump
  subroutine openrsp_info_dump( this_info, io_dump )
    type(rspinfo_t), intent(in) :: this_info
    integer, optional, intent(in) :: io_dump
    ! local IO unit
    integer l_io_dump
    ! number of frequencies
    integer num_freq
    ! temporary stuff
    integer i, j, k
    if ( present( io_dump ) ) then
      l_io_dump = io_dump
    ! uses IO unit of log file to dump the control information
    else
      l_io_dump = this_info%log_io
    end if
    write( l_io_dump, 100 ) 'IO unit of log file: ', this_info%log_io
    write( l_io_dump, 100 ) 'Print level:         ', this_info%level_print
    ! real frequencies
    if ( allocated( this_info%real_freqs ) ) then
      num_freq = size( this_info%real_freqs )
      write( l_io_dump, 100 ) 'Number of specified real frequencies:      ', num_freq
      i = 1
      do while ( i <= num_freq )
        j = min( num_freq, i + 4 )
        write( l_io_dump, 110 ) 'Real frequencies:      ', ( this_info%real_freqs(k), k = i, j )
        i = i + 5
      end do
    end if
    ! imaginary frequencies
    if ( allocated( this_info%imag_freqs ) ) then
      num_freq = size( this_info%imag_freqs )
      write( l_io_dump, 100 ) 'Number of specified imaginary frequencies: ', num_freq
      i = 1
      do while ( i <= num_freq )
        j = min( num_freq, i + 4 )
        write( l_io_dump, 110 ) 'Imaginary frequencies: ', ( this_info%imag_freqs(k), k = i, j )
        i = i + 5
      end do
    end if
    if ( this_info%openrsp_efgb ) write( l_io_dump, 100 )     &
      'Calculate electric-field-gradient-induced (Buckingham) birefringence'
    if ( this_info%openrsp_cme ) write( l_io_dump, 100 )      &
      'Calculate London Cotton-mouton constant'
    if ( this_info%openrsp_roa ) write( l_io_dump, 100 )      &
      'Calculate Raman optical activity (ROA) properties'
    if ( this_info%openrsp_cars ) write( l_io_dump, 100 )     &
      'Calculate coherent anti-Stokes Raman Scattering (CARS)'
    if ( this_info%openrsp_jones ) write( l_io_dump, 100 )    &
      'Calculate magnetoelectric Jones spectroscopy'
    if ( this_info%openrsp_polariz ) write( l_io_dump, 100 )  &
      'Calculate linear polarizability'
    if ( this_info%openrsp_hypolar2 ) write( l_io_dump, 100 ) &
      'Calculate 1st hyperpolarizability using n+1 rule'
    if ( this_info%openrsp_hypolar ) write( l_io_dump, 100 )  &
      'Calculate 1st hyperpolarizability using 2n+1 rule'
    if ( this_info%openrsp_sechyp3 ) write( l_io_dump, 100 )  &
      'Calculate 2nd hyperpolarizability using n+1 rule'
    if ( this_info%openrsp_sechyp ) write( l_io_dump, 100 )   &
      'Calculate 2nd hyperpolarizability using 2n+1 (1+2+1) rule'
    if ( this_info%openrsp_sechyp1 ) write( l_io_dump, 100 )  &
      'Calculate 2nd hyperpolarizability using 2n+1 (2+1+1) rule'
    if ( this_info%openrsp_vibbeta ) write( l_io_dump, 100 )  &
      'Calculate vibrational hyperpolarizability'
    if ( this_info%openrsp_vibshyp ) write( l_io_dump, 100 )  &
      'Calculate vibrational second hyperpolarizability'
100 format('INFO ',A,I6)
110 format('INFO ',A,10F12.6)
  end subroutine openrsp_info_dump

  !> \brief cleans the control information of openrsp
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param this_info contains the control information
  subroutine openrsp_info_clean( this_info )
    type(rspinfo_t), intent(inout) :: this_info
    ! cleans the real frequencies
    if ( allocated( this_info%real_freqs ) ) deallocate( this_info%real_freqs )
    ! cleans the imaginary frequencies
    if ( allocated( this_info%imag_freqs ) ) deallocate( this_info%imag_freqs )
  end subroutine openrsp_info_clean

  !> \brief performs the calculations asked
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param molcfg contains molecule, energy, integral and solver config
  !> \param S contains the overlap matrix
  !> \param D contains the unperturbed density matrix
  !> \param F contains the unperturbed Fock matrix
  !> \param this_info contains the control information of openrsp
  !> \todo adds some comments when calling Andreas' codes
  subroutine openrsp_prop_calc( molcfg, S, D, F, this_info )
    !> \todo
    !> \todo type(molecule_t), intent(in) :: this_mol
    use prop_contribs_old, only: prop_molcfg
    use dalton_ifc
    type(prop_molcfg), intent(in) :: molcfg
    type(matrix), intent(in) :: S
    type(matrix), intent(in) :: D
    type(matrix), intent(in) :: F
    type(rspinfo_t), intent(in) :: this_info
    ! number of atoms
    integer num_atoms
    ! number of coordinates
    integer num_coord
    ! scratch for property tensors
    complex(8), allocatable :: tsr(:)
    ! error information
    integer ierr
    ! imaginary unit
    complex(8), parameter :: imag_one = (0.0D+00,1.0D+00)
    ! zero
    real(8), parameter :: zero = 0.0D+00
    ! complex one
    complex(8), parameter :: cplx_one = (1.0D+00,0.0D+00)
    ! incremental recorder over frequencies
    integer iw
    ! temporary stuff
    integer i, j, k, l

#ifdef OPENRSP_STANDALONE
    print *, 'error: not part of standalone'
    stop 1
#else /* OPENRSP_STANDALONE */
    if (get_is_ks_calculation()) then
       call interface_ao_write()
    end if
#endif /* OPENRSP_STANDALONE */

    ! calculates electric-field-gradient-induced (Buckingham) birefringence
    if ( this_info%openrsp_efgb ) then
      do iw = 1, size( this_info%real_freqs )
        if ( this_info%level_print >= 10 ) &
          write( this_info%log_io, 100 ) 'Calling efgb_Jpri_Bten_Bcal w=', this_info%real_freqs(iw)
        allocate( tsr(3+6+9+9+9+9+18+27+27+27+27+54+54), stat=ierr ) !ajt: This tsr() thing is a hack.
        if ( ierr /= 0 ) then
           print *, 'Failed to allocate tsr!'
           stop 1
        end if
        ! response functions should perhaps eventually be stored in a dedicated module
        call efgb_Jpri_Bten_Bcal( molcfg, S, D, F,                                        &
                                  this_info%real_freqs(iw)*cplx_one*(/1,-1,0/),           &
                                  tsr(1:3), tsr(4:9), tsr(10:18), tsr(19:27), tsr(28:36), &
                                  tsr(37:45), tsr(46:63), tsr(64:90), tsr(91:117),        &
                                  tsr(118:144), tsr(145:171), tsr(172:225), tsr(226:279) )
        if ( this_info%level_print >= 10 ) then
          write( this_info%log_io, 100 ) 'Backing from efgb_Jpri_Bten_Bcal w=', this_info%real_freqs(iw)
          write( this_info%log_io, 100 ) 'Calling efgb_output w=', this_info%real_freqs(iw)
        end if
        call efgb_output( this_info%real_freqs(iw),                               &
                          dreal(tsr(1:3)), dreal(tsr(4:9)), dreal(tsr(10:18)),    &
                          (/((dreal(tsr(19+i+3*j)),j=0,2),i=0,2)/),               &
                          (/((dreal(tsr(28+i+3*j)),j=0,2),i=0,2)/),               &
                          (/((dreal(tsr(37+i+3*j)),j=0,2),i=0,2)/),               &
                          (/((dreal(tsr(46+i+6*j)),j=0,2),i=0,5)/),               &
                          (/(((dreal(tsr(64 +i+3*j+9*k)),k=0,2),j=0,2),i=0,2)/),  &
                          (/(((dreal(tsr(91 +i+3*j+9*k)),k=0,2),j=0,2),i=0,2)/),  &
                          (/(((dreal(tsr(118+i+3*j+9*k)),k=0,2),j=0,2),i=0,2)/),  &
                          (/(((dreal(tsr(145+i+3*j+9*k)),k=0,2),j=0,2),i=0,2)/),  &
                          (/(((dreal(tsr(172+i+6*j+18*k)),k=0,2),j=0,2),i=0,5)/), &
                          dreal(tsr(226:279)), this_info%log_io )
        deallocate( tsr )
        if ( this_info%level_print >= 10 ) &
          write( this_info%log_io, 100 ) 'Backing from efgb_output w=', this_info%real_freqs(iw)
      end do
    end if
    ! calculates London Cotton-mouton constant
    if ( this_info%openrsp_cme ) then
    end if
    ! calculates Raman optical activity (ROA) properties
    if ( this_info%openrsp_roa ) then
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Calling roa_pol_Gpri_Aten_grads, freq =', this_info%real_freqs(1)
      call get_natoms( num_atoms )
      num_coord = 3*num_atoms
      allocate( tsr(3*3+3*3+6*3+num_coord*3*3+num_coord*3*3+num_coord*6*3), stat=ierr ) !ajt: hack
      if ( ierr /= 0 ) then
         print *, 'Failed to allocate tsr!'
         stop 1
      end if
      call roa_pol_Gpri_Aten_grads( molcfg, S, D, F, num_coord,          &
                                    this_info%real_freqs(1)*cplx_one,    &
                                    tsr(1:9), tsr(10:18), tsr(19:36),    &
                                    tsr(37:36+9*num_coord),              &
                                    tsr(37+9*num_coord:36+18*num_coord), &
                                    tsr(37+18*num_coord:36+36*num_coord) )
      if ( this_info%level_print >= 10 ) then
        write( this_info%log_io, 100 ) 'Backing from roa_pol_Gpri_Aten_grads ...'
        write( this_info%log_io, 100 ) 'Calling VIBCTL_ifc ...'
      end if
      call VIBCTL_ifc( num_coord, this_info%real_freqs(1),                           &
                       (/((dreal(tsr(1 +i+3*j)),j=0,2),i=0,2)/),                     &
                       (/((dreal(tsr(10+i+3*j)),j=0,2),i=0,2)/),                     &
                       (/((dreal(tsr(19+i+6*j)),j=0,5),i=0,2)/),                     &
                       (/(zero,i=1,3*num_coord)/),                                   &
                       (/(((dreal(tsr(37+i+num_coord*j+num_coord*3*k)),              &
                            k=0,2),j=0,2),i=0,num_coord-1)/),                        &
                       (/(((dreal(tsr(37+9*num_coord+i+num_coord*j+num_coord*3*k)),  &
                            k=0,2),j=0,2),i=0,num_coord-1)/),                        &
                       (/(((dreal(tsr(37+18*num_coord+i+num_coord*j+num_coord*6*k)), &
                            k=0,2),j=0,5),i=0,num_coord-1)/) )
      deallocate( tsr )
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Backing from VIBCTL_ifc ...'
    end if
    ! calculates coherent anti-Stokes Raman Scattering (CARS)
    if ( this_info%openrsp_cars ) then
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Calling cars_pol_shyp_polgra, freq =', this_info%real_freqs(1)
      ! gets the number of atoms
      !> \todo call xmolecule_get_natom( this_mol, num_atoms )
      call get_natoms( num_atoms )
      num_coord = 3*num_atoms
      allocate( tsr(3+3*3+3*3*3*3+num_coord*3*3), stat=ierr ) !ajt: hack
      if ( ierr /= 0 ) then
         print *, 'Failed to allocate tsr!'
         stop 1
      end if
      call cars_pol_shyp_polgra( molcfg, S, D, F, num_coord, this_info%real_freqs(1)*cplx_one, &
                                 tsr(1:3), tsr(4:12), tsr(13:93), tsr(94:93+9*num_coord) )
      if ( this_info%level_print >= 10 ) then
        write( this_info%log_io, 100 ) 'Backing from cars_pol_shyp_polgra ...'
        write( this_info%log_io, 100 ) 'Calling print_shypol ...'
      end if
      call print_shypol( this_info%real_freqs(1)*(/-1,1,-1,1/), dreal(tsr(1:3)),  &
                         (/dreal(tsr(4:12)),(zero,i=1,9*3)/),(/(zero,i=1,27*6)/), &
                         dreal(tsr(13:93)), this_info%log_io )
      if ( this_info%level_print >= 10 ) then
        write( this_info%log_io, 100 ) 'Backing from print_shypol ...'
        write( this_info%log_io, 100 ) 'Calling vib_ana_polari ...'
      end if
      call vib_ana_polari( molcfg, this_info%real_freqs(1)*(/-1,1,-1,1/),             &
                           dreal(tsr(1:3)), num_coord, (/(zero,i=1,3*num_coord*4)/),  &
                          (/((dreal(tsr(94+i+num_coord*j)),j=0,8),i=0,num_coord-1),   &
                            (zero,i=1,9*num_coord),                                   &
                            ((dreal(tsr(94+i+num_coord*j)),j=0,8),i=0,num_coord-1),   &
                            ((dreal(tsr(94+i+num_coord*j)),j=0,8),i=0,num_coord-1),   &
                            (zero,i=1,9*num_coord),                                   &
                            ((dreal(tsr(94+i+num_coord*j)),j=0,8),i=0,num_coord-1)/), &
                           (/(zero,i=1,27*num_coord*4)/), this_info%log_io )
      deallocate( tsr )
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Backing from vib_ana_polari ...'
    end if
    ! calculates magnetoelectric Jones spectroscopy
    if ( this_info%openrsp_jones ) then
      allocate( tsr(3+6+3*3+3*3+3*3+3*3+6*3+3*3+3*3*3+3*3*3+6*3*3+6*3*3 &
                    +3*3*3*3+3*3*3*3+6*3*3*3+6*3*3*3), stat=ierr )
      if ( ierr /= 0 ) then
         print *, 'Failed to allocate tsr!'
         stop 1
      end if
      do iw = 1, size( this_info%real_freqs )
        if ( this_info%level_print >= 10 ) &
          write( this_info%log_io, 100 ) 'Calling cme_jones_eta_apri, freq =', this_info%real_freqs(1)
         call cme_jones_eta_apri( molcfg, S, D, F,                                        &
                                  this_info%real_freqs(iw)*cplx_one*(/1,0,-1,0/),         &
                                  tsr(1:3), tsr(4:9), tsr(10:18), tsr(19:27), tsr(28:36), &
                                  tsr(37:45), tsr(46:63), tsr(64:72), tsr(73:99),         &
                                  tsr(100:126), tsr(127:180), tsr(181:234),               &
                                  tsr(235:315), tsr(316:396), tsr(397:558), tsr(559:720) )
          if ( this_info%level_print >= 10 ) then
            write( this_info%log_io, 100 ) 'Backing from cme_jones_eta_apri ...'
            write( this_info%log_io, 100 ) 'Calling jones_output no-London output ...'
          end if
         !no-London output (.false.)
         call jones_output( this_info%real_freqs(iw)*cplx_one, .false., tsr(1:3),         &
                            (/(((tsr(73 +i+3*j +9*k),k=0,2),j=0,2),i=0,2)/),              &
                            (/(((tsr(127+i+6*j+18*k),k=0,2),j=0,2),i=0,5)/),              &
                            (/((((tsr(235+i+3*j +9*k+27*l),l=0,2),k=0,2),j=0,2),i=0,2)/), &
                            (/((((tsr(397+i+6*j+18*k+54*l),l=0,2),k=0,2),j=0,2),i=0,5)/), &
                            this_info%log_io )
          if ( this_info%level_print >= 10 ) then
            write( this_info%log_io, 100 ) 'Calling jones_output London output ...'
          end if
         !London output (.true.)
         call jones_output( this_info%real_freqs(iw)*cplx_one, .true., tsr(1:3),          &
                            (/(((tsr(100+i+3*j +9*k),k=0,2),j=0,2),i=0,2)/),              &
                            (/(((tsr(181+i+6*j+18*k),k=0,2),j=0,2),i=0,5)/),              &
                            (/((((tsr(316+i+3*j +9*k+27*l),l=0,2),k=0,2),j=0,2),i=0,2)/), &
                            (/((((tsr(559+i+6*j+18*k+54*l),l=0,2),k=0,2),j=0,2),i=0,5)/), &
                            this_info%log_io )
      end do
      deallocate( tsr )
    end if
    ! calculates linear polarizability
    if ( this_info%openrsp_polariz ) then
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Calling elec_polariz ...'
      allocate( tsr(3*3), stat=ierr )
      if ( ierr /= 0 ) then
         print *, 'Failed to allocate tsr!'
         stop 1
      end if
      call elec_polariz( molcfg, S, D, F, this_info%real_freqs(1) &
                         + imag_one*this_info%imag_freqs(1), tsr )
      deallocate( tsr )
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Backing from elec_polariz ...'
    end if
    ! calculates 1st hyperpolarizability using n+1 rule
    if ( this_info%openrsp_hypolar2 ) then
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Calling elec_hypolar ...'
      allocate( tsr(3*3*3), stat=ierr )
      if ( ierr /= 0 ) then
         print *, 'Failed to allocate tsr!'
         stop 1
      end if
      call elec_hypolar( molcfg, S, D, F,                                             &
                         (/-this_info%real_freqs(1)-imag_one*this_info%imag_freqs(1)  &
                           -this_info%real_freqs(2)-imag_one*this_info%imag_freqs(2), &
                            this_info%real_freqs(1)+imag_one*this_info%imag_freqs(1), &
                            this_info%real_freqs(2)+imag_one*this_info%imag_freqs(2) /), tsr )
      deallocate( tsr )
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Backing from elec_hypolar ...'
    end if
    ! calculates 1st hyperpolarizability using 2n+1 rule
    if ( this_info%openrsp_hypolar ) then
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Calling alt_elec_hypol ...'
      allocate( tsr(3*3*3), stat=ierr )
      if ( ierr /= 0 ) then
         print *, 'Failed to allocate tsr!'
         stop 1
      end if
      call alt_elec_hypol( molcfg, S, D, F,                                             &
                           (/-this_info%real_freqs(1)-imag_one*this_info%imag_freqs(1)  &
                             -this_info%real_freqs(2)-imag_one*this_info%imag_freqs(2), &
                              this_info%real_freqs(1)+imag_one*this_info%imag_freqs(1), &
                              this_info%real_freqs(2)+imag_one*this_info%imag_freqs(2) /), tsr )
      deallocate( tsr )
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Backing from alt_elec_hypol ...'
    end if
    ! calculates 2nd hyperpolarizability using n+1 rule
    if ( this_info%openrsp_sechyp3 ) then
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Calling elec_sechyp ...'
           allocate( tsr(3*3*3*3), stat=ierr )
      if ( ierr /= 0 ) then
         print *, 'Failed to allocate tsr!'
         stop 1
      end if
      call elec_sechyp( molcfg, S, D, F,                                             &
                        (/-this_info%real_freqs(1)-imag_one*this_info%imag_freqs(1)  &
                          -this_info%real_freqs(2)-imag_one*this_info%imag_freqs(2)  &
                          -this_info%real_freqs(3)-imag_one*this_info%imag_freqs(3), &
                           this_info%real_freqs(1)+imag_one*this_info%imag_freqs(1), &
                           this_info%real_freqs(2)+imag_one*this_info%imag_freqs(2), &
                           this_info%real_freqs(3)+imag_one*this_info%imag_freqs(3) /), tsr )
      deallocate( tsr )
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Backing from elec_sechyp ...'
    end if
    ! calculates 2nd hyperpolarizability using 2n+1 (1+2+1) rule
    if ( this_info%openrsp_sechyp ) then
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Calling alt_elec_sechyp ...'
           allocate( tsr(3*3*3*3), stat=ierr )
      if ( ierr /= 0 ) then
         print *, 'Failed to allocate tsr!'
         stop 1
      end if
      call alt_elec_sechyp( molcfg, S, D, F,                                             &
                            (/-this_info%real_freqs(1)-imag_one*this_info%imag_freqs(1)  &
                              -this_info%real_freqs(2)-imag_one*this_info%imag_freqs(2)  &
                              -this_info%real_freqs(3)-imag_one*this_info%imag_freqs(3), &
                               this_info%real_freqs(1)+imag_one*this_info%imag_freqs(1), &
                               this_info%real_freqs(2)+imag_one*this_info%imag_freqs(2), &
                               this_info%real_freqs(3)+imag_one*this_info%imag_freqs(3) /), tsr )
      deallocate( tsr )
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Backing from alt_elec_sechyp ...'
    end if
    ! calculates 2nd hyperpolarizability using 2n+1 (2+1+1) rule
    if ( this_info%openrsp_sechyp1 ) then
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Calling alt2_elec_sechyp ...'
           allocate( tsr(3*3*3*3), stat=ierr )
      if ( ierr /= 0 ) then
         print *, 'Failed to allocate tsr!'
         stop 1
      end if
      call alt2_elec_sechyp( molcfg, S, D, F,                                             &
                             (/-this_info%real_freqs(1)-imag_one*this_info%imag_freqs(1)  &
                               -this_info%real_freqs(2)-imag_one*this_info%imag_freqs(2)  &
                               -this_info%real_freqs(3)-imag_one*this_info%imag_freqs(3), &
                                this_info%real_freqs(1)+imag_one*this_info%imag_freqs(1), &
                                this_info%real_freqs(2)+imag_one*this_info%imag_freqs(2), &
                                this_info%real_freqs(3)+imag_one*this_info%imag_freqs(3) /), tsr )
      deallocate( tsr )
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Backing from alt2_elec_sechyp ...'
    end if
    ! calculates vibrational hyperpolarizability
    if ( this_info%openrsp_vibbeta ) then
      call get_natoms( num_atoms )
      num_coord = 3*num_atoms
      allocate( tsr(3+3*3*3+3*3*3+num_coord*3*3+num_coord*3*3*3), stat=ierr )
      if ( ierr /= 0 ) then
         print *, 'Failed to allocate tsr!'
         stop 1
      end if
      ! static case
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Calling vibhyp_hyp_dipgra_polgra ...'
      call vibhyp_hyp_dipgra_polgra( molcfg, S, D, F, num_coord, cplx_one* &
                                     !(/-this_info%real_freqs(iw+1)    & !-imag_one*this_info%imag_freqs(1)    &
                                     !  -this_info%real_freqs(iw+2),   & !-imag_one*this_info%imag_freqs(2),   &
                                     !   this_info%real_freqs(iw+1),   & !+imag_one*this_info%imag_freqs(1),   &
                                     !   this_info%real_freqs(iw+2)/), & !+imag_one*this_info%imag_freqs(2)/), &
                                     (/zero,zero,zero/),              &
                                     tsr(1:3), tsr(4:30), tsr(31:57), &
                                     tsr(58:57+9*num_coord),          &
                                     tsr(58+9*num_coord:57+36*num_coord) )
      if ( this_info%level_print >= 10 ) then
        write( this_info%log_io, 100 ) 'Backing from vibhyp_hyp_dipgra_polgra ...'
        write( this_info%log_io, 100 ) 'Calling print_shypol ...'
      end if
      !call print_shypol( (/-this_info%real_freqs(iw+1)          &
      !                     -this_info%real_freqs(iw+2),         &
      !                      this_info%real_freqs(iw+1),         &
      !                      this_info%real_freqs(iw+2),zero/),  &
      call print_shypol( (/zero,zero,zero,zero/),               &
                         dreal(tsr(1:3)),                       &
                         (/dreal(tsr(4:30)),(zero,i=1,9)/),     &
                         (/dreal(tsr(31:57)),(zero,i=1,27*5)/), &
                         (/(zero,i=1,81)/), this_info%log_io )
      if ( this_info%level_print >= 10 ) &
        write( this_info%log_io, 100 ) 'Backing from print_shypol ...'
      do iw = 0, size( this_info%real_freqs )-2, 2
        if ( this_info%level_print >= 10 ) &
          write( this_info%log_io, 100 ) 'Calling vib_ana_polari freq=',          &
                                         -sum( this_info%real_freqs(iw+1:iw+2) ), &
                                         this_info%real_freqs(iw+1:iw+2)
        call vib_ana_polari( molcfg, &
                             (/-this_info%real_freqs(iw+1)  & !-imag_one*this_info%imag_freqs(1)  &
                               -this_info%real_freqs(iw+2), & !-imag_one*this_info%imag_freqs(2), &
                                this_info%real_freqs(iw+1), & !+imag_one*this_info%imag_freqs(1), &
                                this_info%real_freqs(iw+2), & !+imag_one*this_info%imag_freqs(2), &
                                zero/),                     & !(zero,zero)/),                     &
                             dreal(tsr(1:3)), num_coord,                                  &
! results from \fn(vibhyp_hyp_dipgra_polgra) is (ng,3,3)
! to \fn(vib_ana_polari) should be (3,ng,3)
                             (/(((dreal(tsr(58+i+num_coord*j+num_coord*3*k)),             &
                                  j=0,2),i=0,num_coord-1),k=0,2),                         &
                               (zero,i=1,num_coord*3*1)/),                                &
! results from \fn(vibhyp_hyp_dipgra_polgra) is (ng,3,3,3)
! to \fn(vib_ana_polari) should be (3,3,ng,3)
                             (/(((dreal(tsr(58+num_coord*9+i+num_coord*j+num_coord*9*k)), &
                                  j=0,8),i=0,num_coord-1),k=0,2),                         &
                               (zero,i=1,num_coord*9*3)/),                                &
                             (/(zero,i=1,num_coord*27*4)/), this_info%log_io )
        if ( this_info%level_print >= 10 ) &
          write( this_info%log_io, 100 ) 'Backing from vib_ana_polari ...'
      end do
      deallocate( tsr )
    end if
    ! calculates vibrational 2nd hyperpolarizability
    if ( this_info%openrsp_vibshyp ) then
      call get_natoms( num_atoms )
      num_coord = 3*num_atoms
      allocate( tsr(3+3*3*3*3+(3*4+3*3*6+3*3*3*4)*num_coord), stat=ierr )
      if ( ierr /= 0 ) then
         print *, 'Failed to allocate tsr!'
         stop 1
      end if
      do iw = 0, size( this_info%real_freqs )-4, 4
        if ( this_info%level_print >= 10 ) &
          write( this_info%log_io, 100 ) 'Calling vibshyp_shyp_dipg_polg_hypg ...'
        call vibshyp_shyp_dipg_polg_hypg( molcfg, S, D, F, num_coord,                &
                                          (1d0,0d0)*this_info%real_freqs(iw+1:iw+4)  &
                                        + (0d0,1d0)*this_info%imag_freqs(iw+1:iw+4), &
                                          tsr(1:3), tsr(4:84),                       &
                                          tsr(85:84+12*num_coord),                   &
                                          tsr(85+12*num_coord:84+66*num_coord),      &
                                          tsr(85+66*num_coord:84+174*num_coord) )
        if ( this_info%level_print >= 10 ) then
          write( this_info%log_io, 100 ) 'Returned from vibshyp_shyp_dipg_polg_hypg ...'
          write( this_info%log_io, 100 ) 'Calling print_shypol ...'
        end if
        !ajt FIXME if freq complex, the nonzero imaginary part is never printed
        call print_shypol( this_info%real_freqs(iw+1:iw+4), &
                           dreal(tsr(1:3)),                 &
                           (/(zero,i=1,3*3*4)/),            &
                           (/(zero,i=1,3*3*3*6)/),          &
                           (/dreal(tsr(4:84))/), this_info%log_io )
        if ( this_info%level_print >= 10 ) &
          write( this_info%log_io, 100 ) 'Returned from print_shypol ...'
        if ( this_info%level_print >= 10 ) &
          write( this_info%log_io, 100 ) 'Calling vib_ana_polari freq =', &
                                         this_info%real_freqs(iw+1:iw+4)
        !ajt FIXME reverse index ordering in vib_ana_polari
        call vib_ana_polari( molcfg, this_info%real_freqs(iw+1:iw+4),        &
                             dreal(tsr(1:3)), num_coord,                     &
                             (/(((dreal(tsr(85+i+num_coord*(j+3*k))),        &
                                  j=0,3-1),i=0,num_coord-1),k=0,4-1)/),      &
                             (/(((dreal(tsr(85+i+num_coord*(12+j+3*3*k))),   &
                                  j=0,3*3-1),i=0,num_coord-1),k=0,6-1)/),    &
                             (/(((dreal(tsr(85+i+num_coord*(66+j+3*3*3*k))), &
                                  j=0,3*3*3-1),i=0,num_coord-1),k=0,4-1)/),  &
                             this_info%log_io )
        if ( this_info%level_print >= 10 ) &
          write( this_info%log_io, 100 ) 'Returned from vib_ana_polari ...'
      end do
      deallocate( tsr )
    end if
100 format('RSPC ',A,10F11.8)
  end subroutine openrsp_prop_calc

end module openrsp_old

