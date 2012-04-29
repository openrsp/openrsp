!!  gen1int: compute derivatives of one-electron integrals using Hermite Gaussians
!!  Copyright 2009 Bin Gao, Andreas Thorvaldsen, Radovan Bast, and Kenneth Ruud
!!
!!  This file is part of gen1int.
!!
!!  gen1int is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  gen1int is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!!  
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with gen1int. If not, see <http://www.gnu.org/licenses/>.
!!
!!  This file contains the interface of DALTON.
!!
!!  2010-09-27, ajt:
!!  * added subroutine NUCAAT_ifc while updating against repo
!!  * split DIPNUC_ifc into DIPNUC_ifc and DPGNUC_ifc
!!
!!  2010-02-27, Bin Gao:
!!  * adds subroutine get_molecule to get the information of molecule
!!
!!  2010-01-16, Bin Gao:
!!  * adds active part of one-electron density matrix (MO) in MO response solver
!!
!!  2010-01-14, Bin Gao:
!!  * fix the dimension of MO coefficient matrix
!!
!!  2009-12-12, Bin Gao:
!!  * add subroutine NUCLEI_ifc
!!  * add subroutine VIBCTL_ifc
!!
!!  2009-12-10, Bin Gao:
!!  * add subroutine get_natoms to get the number of atoms
!!  * add subroutine DIPNUC_ifc to get the nuclear contributions to electric
!!    dipole moments and dipole gradients
!!  * add subroutine QDRNUC_ifc to get the nuclear contribution to quadrupole moments
!!  * add subroutine GRADNN_ifc to get the nuclear repulsion contribution
!!  * add subroutine HESSNN_ifc to get the nuclear contribution to Hessian
!!  * add subroutine di_select_wrk and di_deselect_wrk to allocate and release the memory
!!
!!  2009-12-08, Bin Gao:
!!  * first version

!> \brief interface of DALTON
!> \author Bin Gao
!> \date 2009-12-08
module dalton_ifc

  use matrix_defop   !type matrix with operators

  implicit none

  public dal_ifc_init
  public dal_ifc_finalize
  public di_select_wrk
  public di_deselect_wrk

  public di_get_overlap_and_H1
  public di_read_operator_int
  public di_get_dens
  public di_get_gmat
  public di_twofck
  public di_get_cmo

  public rsp_mosolver_init
  public rsp_mosolver_free
  public rsp_mosolver_dump
  public rsp_mosolver_exec

  public NUCLEI_ifc
  public DIPNUC_ifc
  public QDRNUC_ifc
  public GRADNN_ifc
  public HESSNN_ifc
  public DPGNUC_ifc
  public AATNUC_ifc
  public VIBCTL_ifc

  !> length of the work array
  integer, private, save :: total_f77_work = 0
  !> position of the non-used work array
  integer, private, save :: next_f77_work = 0
  !> amount of the left work array
  integer, private, save :: left_f77_work = 0
  !> DALTON work array pointer
  real(8), public, pointer, save :: dal_work(:)
  !> true for RHF - closed shell or one electron in one active orbital
  logical, save, private :: restrict_scf = .true.
  !> print level
  integer, save, private :: lprt_dal = 5
  !> IO unit of log file
  integer, save, private :: log_dal = 6
  !> PCM
  logical, save, private :: dal_pcm = .false.

  !> control information of MO response solver
  !>
  !> -# maximum number of micro iterations in the iterative solution of
  !>    the frequency independent linear response functions
  integer, save, private :: solver_maxit = 100
  !> -# maximum dimension of the sub-block of the configuration Hessian
  integer, save, private :: solver_maxphp = 0
  !> -# maximum dimension of the reduced space to which new basis vectors are added
  integer, save, private :: solver_mxrm = 400
  !> -# convergence threshold for the solution of the frequency-independent response equations
  real(8), save, private :: solver_thresh = 1.0D-07
  !> -# true for optimal orbital trial vectors in the iterative solution of
  !>    the frequency-dependent linear response equations
  logical, save, private :: solver_optorb = .false.
  !> -# true for CI calculations
  logical, save, private :: solver_ci = .false.
  !> -# true for triplet perturbation operators
  logical, save, private :: solver_triplet = .false.
  !> -# true for excitation energy calculations, false for linear response equations
  logical, save, private :: solver_excit = .false.

  !> coefficients of molecular orbitals
  type(matrix), private, save :: solver_CMO
  !> coefficients of occupied molecular orbitals
  type(matrix), private, save :: solver_CMO_OCC
  !> coefficients of virtual molecular orbitals
  type(matrix), private, save :: solver_CMO_VIR
  !> active part of one-electron density matrix (MO)
  real(8), allocatable, save :: solver_DV(:)

  contains

  !> \brief initializes the interface of DALTON
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param WORK contains the work memory
  !> \param LWORK is the size of the work memory
  !> \param log_io is the IO unit of log file
  !> \param level_print is the print level
  subroutine dal_ifc_init( WORK, LWORK, log_io, level_print, WAVPCM )
    implicit integer (i,m-n)
#include <implicit.h>
    real(8), target :: WORK(:)
    integer LWORK
    integer, optional, intent(in) :: log_io
    integer, optional, intent(in) :: level_print
    logical, optional, intent(in) :: WAVPCM
    ! uses NASHT
#include <inforb.h>
    ! uses LUPRI: pre-defined unit numbers
#include <priunit.h>
    if ( present( WAVPCM ) ) dal_pcm = WAVPCM
    total_f77_work = LWORK
    next_f77_work = 1
    left_f77_work = LWORK
    dal_work => WORK
    restrict_scf = ( NASHT == 0 )
    if ( present( log_io ) ) then
      log_dal = log_io
    else
      log_dal = LUPRI
    end if
    if ( present( level_print ) ) then
      lprt_dal = level_print
    else
      lprt_dal = 10
    end if
  end subroutine dal_ifc_init


  !> \brief finalizes the interface of DALTON
  !> \author Bin Gao
  !> \date 2009-12-08
  subroutine dal_ifc_finalize
    if ( dal_pcm ) then
#ifdef USE_WAVPCM
      call pcm_finalize
#endif
    end if
    total_f77_work = 0
    next_f77_work = 0
    left_f77_work = 0
    nullify( dal_work )
  end subroutine dal_ifc_finalize


  !> \brief allocates the memory asked
  !> \author Bin Gao
  !> \date 2009-12-10
  !> \param mem_req specifies the amount of memory asked
  !> \return wrk is the asked memory if sucessfully allocated
  subroutine di_select_wrk( wrk, mem_req )
    real(8), pointer, intent(inout) :: wrk(:)
    integer, intent(in) :: mem_req
    if ( mem_req > left_f77_work ) then
      !> \todo call di_create_wrk( wrk, mem_req )
      call STOPIT( 'DALTON_IFC', 'di_select_wrk', next_f77_work+mem_req-1, total_f77_work )
    else
      wrk => dal_work
      next_f77_work = next_f77_work + mem_req
      left_f77_work = total_f77_work - next_f77_work + 1
    end if
  end subroutine di_select_wrk


  !> \brief releases the memory asked
  !> \author Bin Gao
  !> \date 2009-12-10
  !> \param mem_req specifies the amount of memory
  !> \return wrk is the  memory to be released
  subroutine di_deselect_wrk( wrk, mem_req )
    real(8), pointer, intent(inout) :: wrk(:)
    integer, intent(in) :: mem_req
    !> \todo if ( mem_req > left_f77_work ) then
    !> \todo   call di_delete_wrk( wrk )
    !> \todo else
      nullify( wrk )
      left_f77_work = left_f77_work + mem_req
      next_f77_work = next_f77_work - mem_req
    !> \todo end if
  end subroutine di_deselect_wrk


  !> \brief gets the overlap and one-electron Hamiltonian matrices
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \return S contains the overlap matrix
  !> \return H1 contains the one-electron Hamiltonian matrix
  subroutine di_get_overlap_and_H1( S, H1 )
    implicit integer (i,m-n)
#include <implicit.h>
    type(matrix), intent(inout) :: S
    type(matrix), intent(inout) :: H1
    ! uses NBAST, NNBASX, N2BASX
#include <inforb.h>
    ! IO units, use LUPROP for file AOPROPER
#include<inftap.h>
    ! starts of the overlap and one electron Hamiltonian matrices in work memory
    integer work_ovlp, work_ham1
    ! PCM one-electron contributions
    integer work_pcm
    ! integer constants
    real(8), parameter :: one = 1.0D+00
    ! dummy stuff
    integer idummy
    ! external DALTON function finding the corresponding label
    logical FNDLAB
    ! reads overlap and one electron Hamiltonian matrices by calling RDONEL
    work_ovlp = next_f77_work
    work_ham1 = work_ovlp + NNBASX
    next_f77_work = work_ham1 + NNBASX
    left_f77_work = total_f77_work - next_f77_work + 1

    if ( left_f77_work < 0 ) call STOPIT( 'DALTON_IFC', 'di_get_SH1', next_f77_work-1, total_f77_work )
    call RDONEL( 'OVERLAP', .true., dal_work(work_ovlp), NNBASX )
    call RDONEL( 'ONEHAMIL', .true., dal_work(work_ham1), NNBASX )
    ! PCM one-electron contributions    
    if ( dal_pcm ) then
      work_pcm = next_f77_work
      next_f77_work = work_pcm + NNBASX
      left_f77_work = total_f77_work - next_f77_work + 1
      if ( left_f77_work < 0 ) call STOPIT( 'DALTON_IFC', 'di_get_SH1', next_f77_work-1, total_f77_work )
#ifdef USE_WAVPCM      
      call pcm_ao_rsp_1elfock( dal_work(work_pcm) )
      ! adds to one-electron Hamiltonian
      dal_work( work_ham1 : work_ham1 + NNBASX - 1 ) &
                    = dal_work( work_ham1 : work_ham1 + NNBASX - 1 ) &
                    + dal_work( work_pcm  : work_pcm  + NNBASX - 1 )
#endif
      ! cleans
      next_f77_work = work_pcm
      left_f77_work = left_f77_work + NNBASX
    end if
    ! fills the data into matrices S and H1
    !N N2BASX = NBAST * NBAST
    left_f77_work = left_f77_work - N2BASX
    if ( left_f77_work < 0 ) call STOPIT( 'DALTON_IFC', 'DSPTSI', next_f77_work+N2BASX-1, total_f77_work )
    ! gets S
    call DSPTSI( NBAST, dal_work(work_ovlp), S%elms )
    ! gets H1
    call DSPTSI( NBAST, dal_work(work_ham1), H1%elms )
    ! clean
    next_f77_work = work_ovlp
    left_f77_work = total_f77_work - next_f77_work + 1
  end subroutine di_get_overlap_and_H1



  !> \brief gets property integrals from file AOPROPER
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param prop_lab is the label of integral
  !> \param init_prop indicates if initialize the integral matrix
  !> \return prop_int contains the integral matrix
  subroutine di_read_operator_int( prop_lab, prop_int )
    implicit integer (i,m-n)
#include <implicit.h>
    character*(8), intent(in) :: prop_lab
    type(matrix), intent(inout) :: prop_int
    ! uses NBAST, NNBAST, NNBASX, N2BASX
#include <inforb.h>
    ! IO units, use LUPROP for file AOPROPER
#include<inftap.h>
    ! external DALTON function finding the corresponding label
    logical FNDLB2
    ! information when calling subroutine FNDLB2
    character*8 RTNLBL(2)
    ! dummy stuff
    integer IDUMMY
    ! one-electron Hamiltonian
    if ( prop_lab == 'ONEHAMIL' ) then
      call QUIT( 'Not implemented!' )
      call RDONEL( 'ONEHAMIL', ANTI, dal_work(next_f77_work), NNBASX )
      call DSPTSI( NBAST, dal_work(next_f77_work), prop_int%elms )
    else
      ! closes file AOPROPER first
      if ( LUPROP > 0 ) call GPCLOSE( LUPROP, 'KEEP' )
      call GPOPEN( LUPROP, 'AOPROPER', 'OLD', ' ', 'UNFORMATTED', IDUMMY, .false. )
      ! finds the label
      if ( FNDLB2( prop_lab, RTNLBL, LUPROP ) ) then
        ! square matrix
        if ( RTNLBL(2) == 'SQUARE' ) then
          call READT( LUPROP, N2BASX, prop_int%elms )
        ! symmetric matrix
        else if ( RTNLBL(2) == 'SYMMETRI' ) then
          call READT( LUPROP, NNBASX, dal_work(next_f77_work) )
          call DSPTSI( NBAST, dal_work(next_f77_work), prop_int%elms )
        ! anti-symmetric matrix
        else if ( RTNLBL(2) == 'ANTISYMM' ) then
          call READT( LUPROP, NNBASX, dal_work(next_f77_work) )
          call DAPTGE( NBAST, dal_work(next_f77_work), prop_int%elms )
        else
          call QUIT( 'Error: No symmetry label on AOPROPER!' )
        end if
        ! closes file AOPROPER
        call GPCLOSE( LUPROP, 'KEEP' )
      else
        call GPCLOSE( LUPROP, 'KEEP' )
        call QUIT( 'Integrals with label '''//prop_lab// &
                   ''' not found on file ''AOPROPER''!' )
      end if
    end if
  end subroutine di_read_operator_int



  !> \brief gets the AO density matrix
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \return D contains the AO density matrix
  subroutine di_get_dens( D )
    implicit integer (i,m-n)
#include <implicit.h>
    type(matrix), intent(inout) :: D
    ! uses NCMOT, NASHT, NNASHX, N2BASX
#include <inforb.h>
    ! uses LUSIFC
#include <inftap.h>
    ! start of coefficients of molecular orbitals in the work array
    integer strt_cmo
    ! start of active part of one-electron density matrix (MO and AO) in the work array
    integer strt_dv, strt_dvao
    ! if calculating DCAO and DVAO
    logical GETDC, GETDV
    ! dummy stuff
    integer idummy
    ! indicates if found required data from SIRIFC
    logical found
    ! start of coefficients of molecular orbitals
    strt_cmo = next_f77_work
    ! start of active part of one-electron density matrix (MO)
    strt_dv = strt_cmo + NCMOT
    if ( restrict_scf ) then
      ! start of active part of one-electron density matrix (AO)
      strt_dvao = strt_dv + 1
      ! start of left workspace
      next_f77_work = strt_dvao + 1
      ! only calculates DCAO
      GETDC = .true.
      GETDV = .false.
    else
      ! start of active part of one-electron density matrix (AO)
      strt_dvao = strt_dv + NNASHX
      ! start of left workspace
      next_f77_work = strt_dvao + N2BASX
      ! we calculate DCAO and DVAO
      GETDC = .true.
      GETDV = .true.
    end if
    ! size of the left work
    left_f77_work = total_f77_work - next_f77_work + 1
    ! checks if the memory is enough
    if ( left_f77_work < 0 ) call STOPIT( 'DALTON_IFC', 'di_get_dens', next_f77_work-1, total_f77_work )
    ! opens SIRIFC
    if ( LUSIFC <= 0 ) &
      call GPOPEN( LUSIFC, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED', idummy, .false. )
    rewind( LUSIFC )
    ! reads the molecular orbital coefficients
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
    call SZERO( dal_work(strt_cmo), NCMOT )
#else
    call DZERO( dal_work(strt_cmo), NCMOT )
#endif
    call rd_sirifc( 'CMO', found, dal_work(strt_cmo), dal_work(next_f77_work), left_f77_work )
    if ( .not. found ) call QUIT( 'CMO not found on SIRIFC!' )
    ! reads active part of one-electron density matrix (MO)
    if ( GETDV ) then
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      call SZERO( dal_work(strt_dv), NNASHX )
#else
      call DZERO( dal_work(strt_dv), NNASHX )
#endif
      call rd_sirifc( 'DV', found, dal_work(strt_dv), dal_work(next_f77_work), left_f77_work )
      if ( .not. found ) call QUIT( 'DV not found on SIRIFC!' )
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      call SZERO( dal_work(strt_dvao), N2BASX )
#else
      call DZERO( dal_work(strt_dvao), N2BASX )
#endif
    end if
    ! gets the AO density matrix, using
    !
    ! FCKDEN(GETDC,GETDV,DCAO,DVAO,CMO,DV,WRK,LFRSAV)
    ! Input:
    !   GETDC   if true calculate DCAO
    !   GETDV   if true calculate DVAO
    !   CMO(*)  molecular orbital coefficients
    !   DV(*)   active part of one-electron density matrix (over MO's)
    ! Scratch:
    !   WRK(LFRSAV)
    call FCKDEN( GETDC, GETDV, D%elms, dal_work(strt_dvao), dal_work(strt_cmo), &
                 dal_work(strt_dv), dal_work(next_f77_work), left_f77_work )
    ! sums DCAO and DVAO
    if ( GETDV ) &
      D%elms = D%elms + reshape( dal_work(strt_dvao : strt_dvao+N2BASX-1), &
                                 (/D%nrow, D%ncol/) )
    ! clean
    next_f77_work = strt_cmo
    left_f77_work = total_f77_work - next_f77_work + 1
  end subroutine di_get_dens



  !> \brief gets the two electron contribution (G) to Fock matrix
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param D contains the AO density matrix
  !> \return G contains the two electron contribution
  subroutine di_get_gmat( D, G )
    implicit integer (i,m-n)
#include <implicit.h>
    type(matrix), intent(in) :: D
    type(matrix), intent(inout) :: G
    ! uses NBAST, NNBAST, NNBASX, N2BASX
#include <inforb.h>
    ! uses NODC, NODV
    ! since we here prefer "implicit none", we have to
    ! define some variables even we do not need them
    integer IPRINT, IDCOOR, MAXDIF, IDATOM
#include <cbione.h>
    ! use PCM
    ! start of total density matrix (AO) in work memory
    integer work_ao_dens
    ! start of PCM two-electron contributions
    integer work_pcm, work_pcm2
    ! parameters for SIRFCK
    integer NDMAT, ISYMDM, IFCTYP
    ! integer constants
    real(8), parameter :: two = 2.0D+00
    ! dummy stuff
    integer idummy
    real(8) xdummy
    ! assigns the work memory
    work_ao_dens = next_f77_work
    next_f77_work = work_ao_dens + N2BASX
    left_f77_work = total_f77_work - next_f77_work + 1
    if ( left_f77_work < 0 ) call STOPIT( 'DALTON_IFC', 'di_get_gmat', next_f77_work-1, total_f77_work )
    ! sets the total density matrix
    !> \todo this may fail for unrestricted calculations
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
    call SCOPY( N2BASX, D%elms, 1, dal_work(work_ao_dens), 1 )
    call SSCAL( N2BASX, two, dal_work(work_ao_dens), 1 )
#else
    call DCOPY( N2BASX, D%elms, 1, dal_work(work_ao_dens), 1 )
    call DSCAL( N2BASX, two, dal_work(work_ao_dens), 1 )
#endif
    ! outputs the total density matrix to check
    if ( lprt_dal >= 20 ) then
      !> \todo call xdump_array2( dims = (/NBAST,NBAST/),          &
      !> \todo                    darray = dal_work(work_ao_dens), &
      !> \todo                    iout = log_dal,                  &
      !> \todo                    label = 'Total density matrix (AO) in DALTON_IFC' )
    end if
    ! only one density matrix
    NDMAT = 1
    ISYMDM = 1
    !> \todo determines IFCTYP run-time
    IFCTYP = 3
    ! calculates two electron contribution by calling SIRFCK
    call SIRFCK( G%elms, dal_work(work_ao_dens), NDMAT, &
                 ISYMDM, IFCTYP, .true., dal_work(next_f77_work), left_f77_work )
    ! PCM two-electron contributions
    if ( dal_pcm ) then
      ! assigns the work memory
      work_pcm = next_f77_work
      work_pcm2 = work_pcm + NNBASX
      next_f77_work = work_pcm2 + N2BASX
      left_f77_work = total_f77_work - next_f77_work + 1
      if ( left_f77_work < 0 ) call STOPIT( 'DALTON_IFC', 'di_pcmfck2', next_f77_work-1, total_f77_work )
      call WAVPCM_2EL( dal_work(work_pcm), dal_work(work_ao_dens), dal_work(next_f77_work), left_f77_work )
!      call pcm_ao_rsp_2elfock( dal_work(work_pcm), dal_work(work_ao_dens) )

      ! transforms to square matrix
      call DZERO( dal_work(work_pcm2), N2BASX )
      call DSPTSI( NBAST, dal_work(work_pcm), dal_work(work_pcm2) )

      ! adds to G
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      call SAXPY( N2BASX, 1.0, dal_work(work_pcm2), 1, G%elms, 1 )
#else
      call DAXPY( N2BASX, 1D0, dal_work(work_pcm2), 1, G%elms, 1 )
#endif
    end if
    !N if ( .not. restrict_scf ) G%elmsb = G%elms
    ! cleans
    next_f77_work = work_ao_dens
    left_f77_work = total_f77_work - next_f77_work + 1
  end subroutine di_get_gmat



  !> \brief gets the two electron contribution (G) to Fock matrix
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param D contains the AO density matrix
  !> \return G contains the two electron contribution
  subroutine di_twofck( D, G )
    implicit integer (i,m-n)
#include <implicit.h>
    type(matrix), intent(in) :: D
    type(matrix), intent(inout) :: G
    ! uses NBAST, NNBAST, NNBASX, N2BASX
#include <inforb.h>
    ! uses NODC, NODV
    ! since we here prefer "implicit none", we have to
    ! define some variables even we do not need them
    integer IPRINT, IDCOOR, MAXDIF, IDATOM
#include <cbione.h>
    ! start of total density matrix (AO) in work memory
    integer work_ao_dens
    ! parameters for SIRFCK
    integer NDMAT, ISYMDM(2), IFCTYP(2)
    ! integer constants
    real(8), parameter :: two = 2.0D+00
    ! number of elements
    integer nelms
    call di_get_gmat( D, G )
  end subroutine di_twofck



  !> \brief gets the coefficients of molecular orbitals
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \return CMO contains the coefficients of molecular orbitals
  !> \return CMO_OCC contains the coefficients of occupied molecular orbitals
  !> \return CMO_VIR contains the coefficients of virtual molecular orbitals
  subroutine di_get_cmo( CMO, CMO_OCC, CMO_VIR )
    implicit integer (i,m-n)
#include <implicit.h>
    type(matrix), intent(inout) :: CMO
    type(matrix), intent(inout) :: CMO_OCC
    type(matrix), intent(inout) :: CMO_VIR
    ! uses LUSIFC, unit number for SIRFC
#include <inftap.h>
    ! uses NBAST, NNBAST, NNBASX, N2BASX
#include <inforb.h>
    ! uses NODC, NODV
    ! since we here prefer "implicit none", we have to
    ! define some variables even we do not need them
    integer IPRINT, IDCOOR, MAXDIF, IDATOM
#include <cbione.h>
    ! dummy stuff
    integer idummy
    ! temparary stuff
    logical found
    integer i
    ! opens SIRIFC
    if ( LUSIFC <= 0 ) &
      call GPOPEN( LUSIFC, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED', idummy, .false. )
    rewind( LUSIFC )
    ! reads the molecular orbital coefficients
    call rd_sirifc( 'CMO', found, CMO%elms, dal_work(next_f77_work), left_f77_work )
    if ( .not. found ) call QUIT( 'CMO not found on SIRIFC!' )
    !N if ( .not. restrict_scf ) CMO%elmsb = CMO%elms
    ! generates the occupied and virtual molecular orbitals
    CMO_OCC%elms(:,:NOCCT)   = CMO%elms(:,:NOCCT)
    CMO_OCC%elms(:,NOCCT+1:) = 0
    CMO_VIR%elms(:,:NOCCT)   = 0
    CMO_VIR%elms(:,NOCCT+1:) = CMO%elms(:,NOCCT+1:)
    ! closes SIRIFC
    if ( LUSIFC > 0 ) call GPCLOSE( LUSIFC, 'KEEP' )
  end subroutine di_get_cmo



  !> \brief initializes the MO response solver
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param max_num_iterat is the maximum number of micro iterations of solving linear response functions
  !> \param max_dim_hess is the maximum dimension of the sub-block of the configuration Hessian
  !> \param max_dim_reduc is the maximum dimension of the reduced space to which new basis vectors are added
  !> \param threshold is the convergence threshold of solving response equations
  !> \param optimal_orb indicates if using optimal orbital trial vectors of solving response equations
  subroutine rsp_mosolver_init( max_num_iterat, max_dim_hess, max_dim_reduc, &
                                threshold, optimal_orb )
    ! need certain low-level matrix routines
    ! radovan: accessing low lever routines here is bad
    use matrix_backend, only: mat_nullify, &
                              mat_setup, &
                              mat_magic_setup, &
                              mat_alloc
    implicit integer (i,m-n)
#include <implicit.h>
    integer, optional, intent(in) :: max_num_iterat
    integer, optional, intent(in) :: max_dim_hess
    integer, optional, intent(in) :: max_dim_reduc
    real(8), optional, intent(in) :: threshold
    logical, optional, intent(in) :: optimal_orb
    ! uses NBAST
#include <inforb.h>
    ! uses LUSIFC
#include <inftap.h>
    ! error information
    integer ierr
    ! if found information on SIRIFC
    logical found
#ifdef G1INT_DEBUG
    ! pushes current subroutine into the stack
    call xsub_enter('rsp_mosolver_init')
#endif
    ! sets the maximum number of micro iterations of solving the response equations
    if ( present( max_num_iterat ) ) solver_maxit = max_num_iterat
    ! sets the maximum dimension of the sub-block of the configuration Hessian
    if ( present( max_dim_hess ) ) solver_maxphp = max_dim_hess
    ! sets the maximum dimension of the reduced space to which new basis vectors are added
    if ( present( max_dim_reduc ) ) solver_mxrm = max_dim_reduc
    ! sets the convergence threshold of solving the response equations
    if ( present( threshold ) ) solver_thresh = threshold
    ! if using the optimal orbital trial vectors in solving response equations
    if ( present( optimal_orb ) ) solver_optorb = optimal_orb
    ! initializes the coefficients of molecular orbitals matrices
    call mat_nullify( solver_CMO )
    solver_CMO%nrow  = NBAST
    solver_CMO%ncol  = solver_CMO%nrow !also virtual
    solver_CMO%closed_shell = .true.
    solver_CMO%magic_tag = mat_magic_setup
    call mat_alloc( solver_CMO )
    call mat_nullify( solver_CMO_OCC )
    call mat_setup( solver_CMO_OCC, solver_CMO )
    call mat_alloc( solver_CMO_OCC )
    call mat_nullify( solver_CMO_VIR )
    call mat_setup( solver_CMO_VIR, solver_CMO )
    call mat_alloc( solver_CMO_VIR )
    ! gets the coefficients of molecular orbitals
    call di_get_cmo( solver_CMO, solver_CMO_OCC, solver_CMO_VIR )
    ! reads active part of one-electron density matrix (MO)
    if ( restrict_scf ) then
      allocate( solver_DV( NNASHX ), stat=ierr )
      if ( ierr /= 0 ) call QUIT( 'Failed to allcoate solver_DV!' )
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      call SZERO( solver_DV, NNASHX )
#else
      call DZERO( solver_DV, NNASHX )
#endif
      ! opens SIRIFC
      if ( LUSIFC <= 0 ) &
        call GPOPEN( LUSIFC, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED', ierr, .false. )
      rewind( LUSIFC )
      call rd_sirifc( 'DV', found, solver_DV, dal_work(next_f77_work), left_f77_work )
      if ( .not. found ) call QUIT( 'DV not found on SIRIFC!' )
    else
      allocate( solver_DV(1), stat=ierr )
      if ( ierr /= 0 ) call QUIT( 'Failed to allcoate solver_DV!' )
      solver_DV(1) = 0.0D+00
    end if
#ifdef G1INT_DEBUG
    ! pops the stack
    call xsub_leave
#endif
  end subroutine rsp_mosolver_init



  !> \brief cleans the MO response solver
  !> \author Bin Gao
  !> \date 2009-12-08
  subroutine rsp_mosolver_free
    solver_CMO = 0
    solver_CMO_OCC = 0
    solver_CMO_VIR = 0
    deallocate( solver_DV )
  end subroutine rsp_mosolver_free


  !> \brief dumps the control information of MO response solver
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param io_dump is the IO unit to dump
  subroutine rsp_mosolver_dump( io_dump )
    integer, intent(in) :: io_dump
#ifdef G1INT_DEBUG
    ! pushes current subroutine into the stack
    call xsub_enter('rsp_mosolver_dump')
#endif
    write( io_dump, 100 ) &
      'Maximum number of micro iterations of solving the response equations: ', solver_maxit
    write( io_dump, 100 ) &
      'Maximum dimension of the sub-block of the configuration Hessian:      ', solver_maxphp
    write( io_dump, 100 ) &
      'Maximum dimension of the reduced space:                               ', solver_mxrm
    write( io_dump, 110 ) &
      'Convergence threshold of solving the response equations:              ', solver_thresh
    if ( solver_optorb ) write( io_dump, 100 ) &
      'Using the optimal orbital trial vectors in solving the response equations'
    if ( restrict_scf )  write( io_dump, 100 ) &
      'Restricted Hartree-Fock (RHF) calculations'
    write( io_dump, 100 ) 'IO unit of log file: ', log_dal
    write( io_dump, 100 ) 'Print level:         ', lprt_dal
    ! outputs matrices to check
#ifdef G1INT_DEBUG
    ! pops the stack
    call xsub_leave
#endif
100 format('INFO ',A,I6)
110 format('INFO ',A,E16.8)
  end subroutine rsp_mosolver_dump



  !> \brief calls MO response solver in DALTON
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param GD is the right hand side vectors (property gradients?)
  !> \param eigval contains the frequencies
  !> \return eigvec contains the solution vectors (AO)
  subroutine rsp_mosolver_exec( GD, eigval, eigvec )
    ! need certain low-level matrix routines
    use matrix_backend, only: mat_nullify, &
                              mat_setup, &
                              mat_magic_setup, &
                              mat_alloc
    implicit integer (i,m-n)
#include <implicit.h>
    type(matrix), intent(in) :: GD(*)
    real(8), intent(in)  :: eigval(*)
    type(matrix), intent(inout) :: eigvec(*)
    ! right hand side vector (MO)
    type(matrix) :: GD_MO
    ! for MO response solver of DALTON, 1 for real operator, -1 for imaginary operator
    integer, parameter :: solver_nabaty = 2
    ! solution vectors of MO solver
    type(matrix), allocatable :: mo_eigvec(:)
    ! LUGDVE - unit number for right-hand sides
    ! LUSOVE - unit number for solutions
    ! LUREVE - unit number for residuals
    integer LUSOVE, LUGDVE, LUREVE
    ! dummy stuff
    integer idummy
    ! uses NSYM
#include <inforb.h>
    ! uses NVARPT (number of solution vectors?)
    integer IPRLIN, LSYMRF, LSYMPT, LSYMST, &
            NCONRF, NCONST, NWOPPT, NVARPT
#include <inflin.h>
    ! uses LURSP
#include <inftap.h>
    ! uses JWOP(2,MAXWOP) for the orbital rotation
#include <infvar.h>
    !
    integer KMJWOP
    ! uses IRAT
#include <iratdef.h>
    !
    character*8 :: LAB1 = 'MOSOLVER', LAB2 = '        '
    ! uses KZVAR
#include <wrkrsp.h>
    ! constants
    real(8), parameter :: half = 5.0D-01
    real(8), parameter :: zero = 0.0D+00
    ! tempary stuff
    integer :: ISYM = 1
    integer IRHS
    integer ISOL
    integer ierr
    !
    integer :: rsp2_number_of_omegas = 1
    integer :: rsp2_number_of_rhs = 1

    KZWOPT = NOCCT*NVIRT
    KZVAR = KZWOPT
    KZYVAR = KZWOPT + KZWOPT
    KMJWOP = next_f77_work + KZYVAR
    if ( KMJWOP+(16*MAXWOP+1)/IRAT > total_f77_work ) &
      call STOPIT( 'DALTON_IFC', 'SOL(MO)', KMJWOP+(16*MAXWOP+1)/IRAT, total_f77_work )

    ! open files
    LUSOVE = -1
    LUGDVE = -1
    LUREVE = -1
    call GPOPEN( LUSOVE, ' ', 'UNKNOWN', ' ', ' ', idummy, .false. )
    call GPOPEN( LUGDVE, ' ', 'UNKNOWN', ' ', ' ', idummy, .false. )
    call GPOPEN( LUREVE, ' ', 'UNKNOWN', ' ', ' ', idummy, .false. )

    !N ! loops over symmetry
    !N do ISYM = 1, NSYM

    ! transforms from AO to MO, and writes RHS (MO) into file

    call mat_nullify(GD_MO)
    GD_MO%nrow  = NORBT
    GD_MO%ncol  = GD_MO%nrow
    GD_MO%closed_shell = .true.
    GD_MO%magic_tag = mat_magic_setup
    call mat_alloc( GD_MO )
    do IRHS = 1, rsp2_number_of_rhs
      ! TRANSFORM (ISYM,JSYM) SYMMETRY BLOCK OF THE MATRIX PRPAO
      ! FROM AO SYMMETRY ORBITALS TO MO BASIS
      call UTHV( solver_CMO%elms, GD(IRHS)%elms, solver_CMO%elms, &
                 ISYM, ISYM, NBAST, NBAST, GD_MO%elms, dal_work(next_f77_work) )
      ! DISTRIBUTE PROPERTY MO INTEGRALS INTO GP VECTORS
      call PRPORB( GD_MO%elms, solver_DV, dal_work(next_f77_work) )
      !FIXME: why multiplied by -1
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      call SSCAL( KZVAR, -1.0, dal_work(next_f77_work), 1 )
#else
      call DSCAL( KZVAR, -1.0D+00, dal_work(next_f77_work), 1 )
#endif
      ! writes out right hand side vector
      call WRITT( LUGDVE, KZYVAR, dal_work( next_f77_work : next_f77_work+KZYVAR-1 ) )
      !> \todo ! outputs to check
      !> \todo if ( lprt_dal >= 20 ) then
      !> \todo   call xdump_array2( dims = (/KZVAR,KZVAR/),       &
      !> \todo                      darray = dal_work(next_f77_work), &
      !> \todo                      iout = log_dal,               &
      !> \todo                      label = 'GP Vector (MO) in DALTON_IFC' )
      !> \todo end if
      call DZERO( GD_MO%elms, NORBT*NORBT )
    end do

    ! calculates the linear response vector and writes to file
    !
    ! when doing response calculations, we need to close RSPVEC
    if ( LURSP > 0 ) call GPCLOSE( LURSP, 'KEEP' )
    ! calls ABACUS solver
    call ABARSP( solver_ci, restrict_scf, solver_triplet, solver_optorb, &
                 ISYM, solver_excit, eigval, rsp2_number_of_omegas,      &
                 solver_nabaty, rsp2_number_of_rhs, LAB1,                &
                 LUGDVE, LUSOVE, LUREVE, solver_thresh, solver_maxit,    &
                 lprt_dal, solver_mxrm, solver_maxphp,                   &
                 dal_work(next_f77_work), left_f77_work )

    ! reads the MO solutions and residuals
    rewind( LUSOVE )
    !N rewind( LUREVE )
    allocate( mo_eigvec(rsp2_number_of_omegas), stat=ierr )
    if ( ierr /= 0 ) call QUIT( 'Failed to allocate MO solutions!' )
    ! loops over solution vectors
    do ISOL = 1, rsp2_number_of_omegas
      call mat_nullify( mo_eigvec(ISOL) )
      call mat_setup( mo_eigvec(ISOL), GD_MO )
      call mat_alloc( mo_eigvec(ISOL) )
      ! reads the solution
      call READT( LUSOVE, KZYVAR, dal_work(next_f77_work) )

      ! JWOP(1,i): inactive (i)
      ! JWOP(2,i): secondary (a)
      !    i  a
      ! i [0  k*]
      ! a [k  0 ]
      call SETZY( dal_work(KMJWOP) )
      ! This subroutine unpacks the ZY matrix from the vector.
      ! It uses the Z and the Y part of the vector.
      call GTZYMT( 1, dal_work(next_f77_work), KZYVAR, ISYM, mo_eigvec(ISOL)%elms, dal_work(KMJWOP) )
      ! divides solution by 2 in accordance with ABACUS solver, or
      ! because Andreas' code does not use total density matrix
      mo_eigvec(ISOL)%elms = mo_eigvec(ISOL)%elms / 2
      ! transforms from MO to AO
      eigvec(ISOL) = - solver_CMO_OCC*( mo_eigvec(ISOL)*trps( solver_CMO_VIR ) ) &
                     - solver_CMO_VIR*( mo_eigvec(ISOL)*trps( solver_CMO_OCC ) )
    end do ! loops over solution vectors

    ! closes and deletes files
    if ( LUSOVE > 0 ) call GPCLOSE( LUSOVE, 'DELETE' )
    if ( LUGDVE > 0 ) call GPCLOSE( LUGDVE, 'DELETE' )
    if ( LUREVE > 0 ) call GPCLOSE( LUREVE, 'DELETE' )
    ! cleans
    GD_MO = 0
    do ISOL = 1, rsp2_number_of_omegas
      mo_eigvec(ISOL) = 0
    end do
    deallocate( mo_eigvec )
  end subroutine rsp_mosolver_exec



  !> \brief gets the information of atoms
  !> \author Bin Gao
  !> \date 2009-12-12
  !> \param n is the number of atoms
  !> \return Z contains the charges
  !> \return IS contains the ...
  !> \return G contains the coordinates
  !> \note modified on linsca/linears/VIBCTL_interface.F
  subroutine NUCLEI_ifc( n, Z, IS, G )
    implicit integer (i,m-n)
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
    integer, intent(in) :: n
    integer, intent(out) :: IS(n)
    real(8), intent(out) :: Z(n)
    real(8), intent(out) :: G(3,n)
    if ( n /= NATOMS ) call QUIT( 'NUCLEI_ifc: n/=NATOMS!' )
    Z = CHARGE(1:n)
    IS = ISOTOP(1:n)
    G = CORD(:,1:n)
  end subroutine NUCLEI_ifc


  !> \brief gets the nuclear contribution to electric dipole moment
  !> \author Bin Gao
  !> \date 2009-12-10
  !> \return DN contains the nuclear contribution to electric dipole moment
  !> \note modified on linsca/linears/VIBCTL_interface.F
  subroutine DIPNUC_ifc( DN )
    implicit integer (i,m-n)
#include <implicit.h>
    ! uses MXCOOR
#include <mxcent.h>
    ! uses DIPMN and DDIPN
#include <dipole.h>
    real(8), intent(out) :: DN(3)
    real(8), parameter :: zero = 0.0D+00
    call DIPNUC( (/zero/), (/zero/), 0, .true. )
    DN = DIPMN
  end subroutine DIPNUC_ifc


  !> \brief gets the nuclear contribution to quadrupole moments
  !> \author Bin Gao
  !> \date 2009-12-10
  !> \return Q contains the nuclear contribution to quadrupole moments
  !> \note modified on linsca/linears/VIBCTL_interface.F
  subroutine QDRNUC_ifc( Q )
    implicit integer (i,m-n)
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <quadru.h>
#include <priunit.h>
    real(8), intent(out) :: Q(6)
    real(8), parameter :: zero = 0.0D+00
    call NUCQDR( CORD(:,1:NUCDEP), (/zero,zero,zero/), LUPRI, 0 )
    Q = (/QDRNUC(1,1),QDRNUC(1:2,2),QDRNUC(1:3,3)/)
  end subroutine QDRNUC_ifc


  !> \brief gets the nuclear repulsion contribution
  !> \author Bin Gao
  !> \date 2009-12-10
  !> \param n is the number of atoms
  !> \return G contains the nuclear repulsion contribution
  !> \note modified on linsca/linears/VIBCTL_interface.F
  subroutine GRADNN_ifc( n, G )
    implicit integer (i,m-n)
#include <implicit.h>
    ! uses MXCOOR
#include <mxcent.h>
    ! uses IPRINT and MAXDIF
#include <cbinuc.h>
    ! uses GRADNN
#include <energy.h>
    integer, intent(in) :: n
    real(8), intent(out) :: G( 3*n )
    real(8), parameter :: zero = 0.0D+00
    IPRINT = 0
    MAXDIF = 1
    call NUCREP( (/zero/), (/zero/), (/zero/) )
    G(1:3*n) = GRADNN(1:3*n)
  end subroutine GRADNN_ifc


  !> \brief gets the nuclear contribution to Hessian
  !> \author Bin Gao
  !> \date 2009-12-10
  !> \param na is the number of atoms
  !> \return H contains the nuclear contribution to Hessian
  !> \note modified on linsca/linears/VIBCTL_interface.F
  subroutine HESSNN_ifc( na, H )
    implicit integer (i,m-n)
#include <implicit.h>
    integer, intent(in) :: na
    real(8), intent(inout) :: H( 3*na, 3*na )
    ! uses MXCOOR
#include <mxcent.h>
    ! uses IPRINT and MAXDIF
#include <cbinuc.h>
    integer i, j
    real(8) HESSNN( MXCOOR, MXCOOR )
    real(8), parameter :: zero = 0.0D+00
    IPRINT = 0
    MAXDIF = 2
    ! second and third arg only used when IPRINT > 1
    call NUCREP( HESSNN, (/zero/), (/zero/) )
    ! ajt This might come out with halved diagonal
    do j = 1, 3*na
      do i = 1, 3*na
        H(i,j) = HESSNN( max(i,j), min(i,j) )
      end do
    end do
  end subroutine HESSNN_ifc


  !> \brief gets the nuclear contribution to electric dipole gradient
  !> \author Bin Gao
  !> \date 2009-12-10
  !> \param na is the number of atoms
  !> \return DGN contains the nuclear contributions to dipole gradient
  !> \note modified on linsca/linears/VIBCTL_interface.F
  subroutine DPGNUC_ifc( na, DGN )
    implicit integer (i,m-n)
#include <implicit.h>
    ! uses MXCOOR
#include <mxcent.h>
    ! uses DDIPN
#include <dipole.h>
    integer, intent(in) :: na
    real(8), intent(out) :: DGN( 3*na, 3 )
    real(8), parameter :: zero = 0.0D+00
    call DIPNUC( (/zero/), (/zero/), 0, .true. )
    DGN = transpose( DDIPN( :, 1:3*na ) )
  end subroutine DPGNUC_ifc


  !> Nuclear contribution to the atomic axial tenaor (AAT),
  !> needed for vibrational circular dichroism (VCD)
  !> In the quasienergy formalism, the AAT is:
  !> d^3E/dR(-w)dB(w)dw |w=0
  subroutine AATNUC_ifc( na, AATN )
    implicit integer (i,m-n)
#include <implicit.h>
#include <mxcent.h>
#include <aatens.h>
    integer, intent(in) :: na
    real(8), intent(out) :: AATN( 3, 3*na )
    real(8) :: CSTRA( 3*na, 3*na ), SCTRA( 3*na, 3*na )
    call NUCAAT( CSTRA, SCTRA, 0 )
         AATN(:,:) = AATNUC( :, :3*na )
  end subroutine AATNUC_ifc


  !> \brief gets the 
  !> \author Bin Gao
  !> \date 2009-12-12
  !> \param 
  !> \return 
  !> \note modified on linsca/linears/VIBCTL_interface.F
  subroutine VIBCTL_ifc( nc, w, ALPHA, GPRIME, THETA, DIPGRAD, &
                         dALPHAdR, dGPRIMEdR, dTHETAdR )
    implicit integer (i,m-n)
#include <implicit.h>
#include <mxcent.h>
#include <cbilnr.h>
#include <cbivib.h>
#include <abainf.h>
#include <moldip.h>
    dimension ALPHA(3,3), GPRIME(3,3), THETA(3,6)
    dimension DIPGRAD(3,nc), dALPHAdR(3,3,nc)
    dimension dGPRIMEdR(3,3,nc), dTHETAdR(3,6,nc)
    integer lunit, i, j, k
    NFRVAL = 1 !NFRVAL in cbilnr
    FRVAL(1) = w  !FRVAL in cbilnr
    FRVAL(2:MXFR) = 0
    lunit = -1
    call GPOPEN( lunit, 'DALTON.WLK', ' ', 'NEW', 'UNFORMATTED', (0), .false. )
    write (lunit)
    write (lunit) !three blank records
    write (lunit)
    write (lunit) NFRVAL,FRVAL,                                       &
                     ALPHA(:,:)  ,(0d0,i=1,3*3*(MXFR-1)),             & !alpha(x,y,w)
                 (dALPHAdR(:,:,i),(0d0,k=1,3*3*(MXFR-1)),i=1,nc),     & !dalpha/dg(x,y,w,c(1:MXCOOR))
                                  (0d0,i=1,3*3*MXFR*(MXCOOR-nc)),     &
                    GPRIME(:,:)  ,(0d0,i=1,3*3*(MXFR-1)),             & !gprime(x,z,w)
                (dGPRIMEdR(:,:,i),(0d0,k=1,3*3*(MXFR-1)),i=1,nc),     & !dgprime/dg(x,z,w,c(1:MXCOOR))
                                  (0d0,i=1,3*3*MXFR*(MXCOOR-nc)),     &
                    (0d0,i=1,3*3),(0d0,i=1,3*3*(MXFR-1)),             & !gprime_lon(x,z,w)
                   ((0d0,k=1,3*3),(0d0,k=1,3*3*(MXFR-1)),i=1,nc),     & !dgprime_lon/dg(x,z,w,c(1:MXCOOR))
                                  (0d0,i=1,3*3*MXFR*(MXCOOR-nc)),     &
                     THETA(:,1:2),THETA(:,4),THETA(:,2:3),THETA(:,5), &
                     THETA(:,4:6),(0d0,i=1,9*3*(MXFR-1)),             & !theta(x,xy,w)
               (dTHETAdR(:,1:2,i),dTHETAdR(:,4,i),dTHETAdR(:,2:3,i),  &
                dTHETAdR(:,5,i),dTHETAdR(:,4:6,i),                    & !dtheta/dg(x,xy,w,c(1:MXCOOR))
                                  (0d0,k=1,9*3*(MXFR-1)),i=1,nc),     &
                                  (0d0,i=1,9*3*MXFR*(MXCOOR-nc)),     &
                                  (0d0,i=1,9*3*MXFR)
    call GPCLOSE( lunit, 'KEEP' )
    call VIBINI()
    NUMHES = .false.
    HESFIL = .true. !DALTON.HES must be there!
    NINTCM = 0
    RAMAN  = any( dALPHAdR(:,:,:) /= 0 )
    VROA   = ( RAMAN .and. any( dGPRIMEdR(:,:,:) /= 0 ) .and. any( dTHETAdR(:,:,:) /= 0 ) )
    DIPDER = any( DIPGRAD(:,:) /= 0 )
    DIPFLT(:,1:nc) = DIPGRAD
    DIPFLT(:,nc+1:MXCOOR) = 0
    DOSYM(1) = .true.
    IPRINT = 6
    call VIBCTL( dal_work(next_f77_work), left_f77_work )
  end subroutine VIBCTL_ifc


  !> Count the number of contracted Gaussian-type orbital shells
  !> \param ncgto, and number of exponents and contraction coefficents
  !> \param nectr, so that memory for data structures can be allocated
  subroutine SHELLS_find_sizes(ncgto, nectr)
    implicit integer (i,m-n)
#include <implicit.h>
    integer, intent(out) :: ncgto, nectr
    ! need MXSHEL
#include <maxorb.h>
    ! need NLRGSH NBCH NUCO NRCO
#include <shells.h>
    logical haveit(MXSHEL)
    integer i, j
    ! count the number of cgto blocks and the number of cgto
    ! exponents (from +1 below) and contraction coefficients
    ncgto = 0
    nectr = 0
    haveit(:) = .false.
    do i = 1, NLRGSH
       ! if not the first contacted in this block, skip
       if (NUMCF(i) /= 1) cycle
       ! index of AO shell block
       j = NBCH(i)
       if (j <= 0 .or. j > MXSHEL) &
          call quit('SHELLS_sizes error: unexpected NBCH, <=0 or >MXSHEL')
       ! count as cgto
       ncgto = ncgto + 1
       ! if first occurance of this block, count exponents and
       ! contraction coefficients
       if (haveit(j)) cycle
       haveit(j) = .true.
       nectr = nectr + NUCO(i)*(NRCO(i)+1)
    end do
  end subroutine


  subroutine SHELLS_to_type_cgto(ncgto, nectr, ectr, bas)
    use basis_set, only: cgto
    implicit integer (i,m-n)
#include <implicit.h>
    integer,          intent(in)  :: ncgto, nectr
    type(cgto),       intent(out) :: bas(ncgto)
    real(8), target, intent(out) :: ectr(nectr)
    ! need MXSHEL
#include <maxorb.h>
    ! need NLRGSH NBCH NUCO NRCO NUMCF CENT NHKT NSTRT
#include <shells.h>
    ! need MXCONT
#include <aovec.h>
    ! need PRIEXP PRICCF
#include <primit.h>
    ! need MAXQNM
#include <maxmom.h>
    integer alreadyis(MXSHEL)
    real(8) rescal(MXQNM+1)
    integer i, j, k, l, m, ne, nc, maxm
    ! coefficient rescaling factors
    maxm = 0
    rescal(maxm+1) = 2*sqrt(acos(-1d0))
    k = 0 !index in bas
    l = 0 !offset in ectr
    alreadyis(:) = 0
    do i = 1, NLRGSH
       if (NBCH(i) <= 0 .or. NBCH(i) > MXSHEL) &
          call quit('SHELLS_sizes error: unexpected NBCH<=0 or >MXSHEL')
       ! if not the first contracted in this block, skip
       if (NUMCF(i) /= 1) cycle
       ! index of AO shell block
       j = NBCH(i)
       ! fill cgto entry
       k = k + 1
       ne = NUCO(i) !number of exponents
       nc = NRCO(i) !number of contracted
       bas(k)%cent   = CENT(i,:,1)
       bas(k)%mom    = NHKT(i) - 1
       bas(k)%charge = huge(1) !ajt FIXME wrong
       bas(k)%icent  = NCENT(i)
       bas(k)%nbas   = nc * (2*bas(k)%mom + 1)
       bas(k)%ibas   = NSTRT(i)
       if (alreadyis(j) /= 0) then
          bas(k)%exp => bas(alreadyis(j))%exp(:)
          bas(k)%ctr => bas(alreadyis(j))%ctr(:,:)
       else
          alreadyis(j) = k
          ! copy exponents from PRIEXP to ectr
          j = JSTRT(i)
          ectr(l+1:l+ne) = PRIEXP(j+1:j+ne)
          bas(k)%exp => ectr(l+1:l+ne)
          l = l + ne
          ! copy and rescale contraction coefs from PRICCF to ectr
          do while (maxm < bas(k)%mom)
             maxm = maxm + 1
             rescal(maxm+1) = rescal(maxm) / sqrt(2*maxm+1d0)
          end do
          ectr(l+1:l+ne*nc) = rescal(bas(k)%mom+1) * reshape( &
                              transpose(PRICCF(j+1:j+ne,:nc)),(/ne*nc/))
          call point(ectr(l+1:l+ne*nc), bas(k)%ctr)
          l = l + ne*nc
       end if
    end do
  contains
    ! point to re-ranked array
    subroutine point(ctr, ctr_pt)
      real(8), target,  intent(in)  :: ctr(nc,ne)
      real(8), pointer, intent(out) :: ctr_pt(:,:)
      ctr_pt => ctr
    end subroutine
  end subroutine


  !> Move the nuclei and basis functions, as seen from the integral
  !> programs. For doing finite difference differentiation.
  subroutine SHELLS_NUCLEI_displace(ic, dc)
    implicit integer (i,m-n)
#include <implicit.h>
    integer,  intent(in) :: ic
    real(8), intent(in) :: dc !step
    ! need MXSHEL
#include <maxorb.h>
    ! need CENT
#include <shells.h>
    ! need MXCOOR
#include <mxcent.h>
    ! need CORD
#include <nuclei.h>
    integer a, r, i
    a = 1 + (ic-1)/3
    r = 1 + mod(ic-1,3)
    if (a <= 0 .or. a > NATOMS) &
       call quit("SHELLS_NUCLEI_displace error: arg. 'ic' out of range")
    ! move atom in CORD, then in CENT
    CORD(r,a)  = CORD(r,a)  + dc
    do i = 1, NLRGSH
       if (NCENT(i) == a) &
          CENT(i,r,:) = CENT(i,r,:) + dc
    end do
  end subroutine

end module
