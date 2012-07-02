module interface_rsp_solver

   use matrix_defop
   use interface_scf
   use interface_f77_memory
   use interface_io

   implicit none

   public rsp_mosolver_init
   public rsp_mosolver_splash
   public rsp_mosolver_exec
   public rsp_mosolver_finalize

   private

   !> control information of MO response solver
   !>
   !> -# maximum number of micro iterations in the iterative solution of
   !>    the frequency independent linear response functions
   integer   :: solver_maxit = 100
   !> -# maximum dimension of the sub-block of the configuration Hessian
   integer   :: solver_maxphp = 0
   !> -# maximum dimension of the reduced space to which new basis vectors are added
   integer   :: solver_mxrm = 400
   !> -# convergence threshold for the solution of the frequency-independent response equations
   real(8)   :: solver_thresh = 1.0D-07
   !> -# true for optimal orbital trial vectors in the iterative solution of
   !>    the frequency-dependent linear response equations
   logical   :: solver_optorb = .false.
   !> -# true for CI calculations
   logical   :: solver_ci = .false.
   !> -# true for triplet perturbation operators
   logical   :: solver_triplet = .false.
   !> -# true for excitation energy calculations, false for linear response equations
   logical   :: solver_excit = .false.

   !> coefficients of molecular orbitals
   type(matrix) :: solver_CMO
   !> coefficients of occupied molecular orbitals
   type(matrix) :: solver_CMO_OCC
   !> coefficients of virtual molecular orbitals
   type(matrix) :: solver_CMO_VIR
   !> active part of one-electron density matrix (MO)
   real(8), allocatable :: solver_DV(:)

contains

  !> \brief initializes the MO response solver
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param max_num_iterat is the maximum number of micro iterations of solving linear response functions
  !> \param max_dim_hess is the maximum dimension of the sub-block of the configuration Hessian
  !> \param max_dim_reduc is the maximum dimension of the reduced space to which new basis vectors are added
  !> \param threshold is the convergence threshold of solving response equations
  !> \param optimal_orb indicates if using optimal orbital trial vectors of solving response equations
  subroutine rsp_mosolver_init(max_num_iterat, &
                               max_dim_hess,   &
                               max_dim_reduc,  &
                               threshold,      &
                               optimal_orb)

    integer, optional, intent(in) :: max_num_iterat
    integer, optional, intent(in) :: max_dim_hess
    integer, optional, intent(in) :: max_dim_reduc
    real(8), optional, intent(in) :: threshold
    logical, optional, intent(in) :: optimal_orb

#ifdef PRG_DALTON
    ! uses NBAST
#include "inforb.h"
    ! uses LUSIFC
#include "inftap.h"
    ! error information
    integer ierr
    ! if found information on SIRIFC
    logical found
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

    call mat_init(solver_CMO, nrow=NBAST, ncol=NBAST, closed_shell=.true.)

    solver_CMO_OCC = mat_alloc_like(solver_CMO)
    solver_CMO_VIR = mat_alloc_like(solver_CMO)
    ! gets the coefficients of molecular orbitals
    call di_get_cmo( solver_CMO, solver_CMO_OCC, solver_CMO_VIR )
    ! reads active part of one-electron density matrix (MO)
    if ( get_is_restricted_scf_calculation() ) then
      allocate( solver_DV( NNASHX ), stat=ierr )
      if ( ierr /= 0 ) call QUIT( 'Failed to allcoate solver_DV!' )
      call DZERO( solver_DV, NNASHX )
      ! opens SIRIFC
      if ( LUSIFC <= 0 ) &
        call GPOPEN( LUSIFC, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED', ierr, .false. )
      rewind( LUSIFC )
      call rd_sirifc( 'DV', found, solver_DV, f77_memory(get_f77_memory_next()), get_f77_memory_left() )
      if ( .not. found ) call QUIT( 'DV not found on SIRIFC!' )
    else
      allocate( solver_DV(1), stat=ierr )
      if ( ierr /= 0 ) call QUIT( 'Failed to allcoate solver_DV!' )
      solver_DV(1) = 0.0D+00
    end if
#endif /* ifdef PRG_DALTON */

  end subroutine



  !> \brief cleans the MO response solver
  !> \author Bin Gao
  !> \date 2009-12-08
  subroutine rsp_mosolver_finalize()

#ifdef PRG_DALTON
    solver_CMO = 0
    solver_CMO_OCC = 0
    solver_CMO_VIR = 0
    deallocate( solver_DV )
#endif /* ifdef PRG_DALTON */

  end subroutine


  !> \brief dumps the control information of MO response solver
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param io is the IO unit to dump
  subroutine rsp_mosolver_splash( io )
    integer, intent(in) :: io

#ifdef PRG_DALTON
    write( io, 100 ) &
      'Maximum number of micro iterations of solving the response equations: ', solver_maxit
    write( io, 100 ) &
      'Maximum dimension of the sub-block of the configuration Hessian:      ', solver_maxphp
    write( io, 100 ) &
      'Maximum dimension of the reduced space:                               ', solver_mxrm
    write( io, 110 ) &
      'Convergence threshold of solving the response equations:              ', solver_thresh
    if ( solver_optorb ) write( io, 100 ) &
      'Using the optimal orbital trial vectors in solving the response equations'
    if ( get_is_restricted_scf_calculation() )  write( io, 100 ) &
      'Restricted Hartree-Fock (RHF) calculations'
    write( io, 100 ) 'IO unit of log file: ', get_print_unit()
    ! outputs matrices to check
100 format('INFO ',A,I6)
110 format('INFO ',A,E16.8)
#endif /* ifdef PRG_DALTON */

  end subroutine



  !> \brief calls MO response solver in DALTON
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param RHS is the right hand side vectors (property gradients?)
  !> \param eigval contains the frequencies
  !> \return eigvec contains the solution vectors (AO)
  subroutine rsp_mosolver_exec(RHS,    &
                               eigval, &
                               eigvec)

    type(matrix), intent(in)    :: RHS(*)
    real(8),      intent(in)    :: eigval(*)
    type(matrix), intent(inout) :: eigvec(*)

#ifdef PRG_DALTON
    ! right hand side vector (MO)
    type(matrix) :: RHS_MO
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
#include "inforb.h"
    ! uses NVARPT (number of solution vectors?)
    integer IPRLIN, LSYMRF, LSYMPT, LSYMST, &
            NCONRF, NCONST, NWOPPT, NVARPT

#include "inflin.h"
    ! uses LURSP
#include "inftap.h"
    ! uses JWOP(2,MAXWOP) for the orbital rotation
#include "infvar.h"
    !
    integer KMJWOP
    ! uses IRAT
#include "iratdef.h"
    !
    character*8 :: LAB1 = 'MOSOLVER', LAB2 = '        '
    ! uses KZVAR
#include "wrkrsp.h"
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
    KMJWOP = get_f77_memory_next() + KZYVAR
    if ( KMJWOP+(16*MAXWOP+1)/IRAT > get_f77_memory_total() ) &
      call STOPIT( 'DALTON_IFC', 'SOL(MO)', KMJWOP+(16*MAXWOP+1)/IRAT, get_f77_memory_total() )

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

    call mat_init(RHS_MO, nrow=NORBT, ncol=NORBT, closed_shell=.true.)
    do IRHS = 1, rsp2_number_of_rhs
      ! TRANSFORM (ISYM,JSYM) SYMMETRY BLOCK OF THE MATRIX PRPAO
      ! FROM AO SYMMETRY ORBITALS TO MO BASIS
      call UTHV( solver_CMO%elms, RHS(IRHS)%elms, solver_CMO%elms, &
                 ISYM, ISYM, NBAST, NBAST, RHS_MO%elms, f77_memory(get_f77_memory_next()) )
      ! DISTRIBUTE PROPERTY MO INTEGRALS INTO GP VECTORS
      call PRPORB( RHS_MO%elms, solver_DV, f77_memory(get_f77_memory_next()) )
      !FIXME: why multiplied by -1
      call DSCAL( KZVAR, -1.0D+00, f77_memory(get_f77_memory_next()), 1 )
      ! writes out right hand side vector
      call WRITT( LUGDVE, KZYVAR, f77_memory( get_f77_memory_next() : get_f77_memory_next()+KZYVAR-1 ) )
      !> \todo ! outputs to check
      !> \todo if ( lprt_dal >= 20 ) then
      !> \todo   call xdump_array2( dims = (/KZVAR,KZVAR/),       &
      !> \todo                      darray = f77_memory(get_f77_memory_next()), &
      !> \todo                      iout = log_dal,               &
      !> \todo                      label = 'GP Vector (MO) in DALTON_IFC' )
      !> \todo end if
      call DZERO( RHS_MO%elms, NORBT*NORBT )
    end do

    ! calculates the linear response vector and writes to file
    !
    ! when doing response calculations, we need to close RSPVEC
    if ( LURSP > 0 ) call GPCLOSE( LURSP, 'KEEP' )
    ! calls ABACUS solver
    call ABARSP( solver_ci, get_is_restricted_scf_calculation(), solver_triplet, solver_optorb, &
                 ISYM, solver_excit, eigval, rsp2_number_of_omegas,      &
                 solver_nabaty, rsp2_number_of_rhs, LAB1,                &
                 LUGDVE, LUSOVE, LUREVE, solver_thresh, solver_maxit,    &
                 1, solver_mxrm, solver_maxphp,                   &
                 f77_memory(get_f77_memory_next()), get_f77_memory_left() )

    ! reads the MO solutions and residuals
    rewind( LUSOVE )
    !N rewind( LUREVE )
    allocate( mo_eigvec(rsp2_number_of_omegas), stat=ierr )
    if ( ierr /= 0 ) call QUIT( 'Failed to allocate MO solutions!' )
    ! loops over solution vectors
    do ISOL = 1, rsp2_number_of_omegas
      mo_eigvec(ISOL) = mat_alloc_like(RHS_MO)
      ! reads the solution
      call READT( LUSOVE, KZYVAR, f77_memory(get_f77_memory_next()) )

      ! JWOP(1,i): inactive (i)
      ! JWOP(2,i): secondary (a)
      !    i  a
      ! i [0  k*]
      ! a [k  0 ]
      call SETZY( f77_memory(KMJWOP) )
      ! This subroutine unpacks the ZY matrix from the vector.
      ! It uses the Z and the Y part of the vector.
      call GTZYMT( 1, f77_memory(get_f77_memory_next()), KZYVAR, ISYM, mo_eigvec(ISOL)%elms, f77_memory(KMJWOP) )
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
    RHS_MO = 0
    do ISOL = 1, rsp2_number_of_omegas
      mo_eigvec(ISOL) = 0
    end do
    deallocate( mo_eigvec )
#endif /* ifdef PRG_DALTON */

  end subroutine

#ifdef PRG_DALTON
  !> \brief gets the coefficients of molecular orbitals
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \return CMO contains the coefficients of molecular orbitals
  !> \return CMO_OCC contains the coefficients of occupied molecular orbitals
  !> \return CMO_VIR contains the coefficients of virtual molecular orbitals
  subroutine di_get_cmo( CMO, CMO_OCC, CMO_VIR )

    type(matrix), intent(inout) :: CMO
    type(matrix), intent(inout) :: CMO_OCC
    type(matrix), intent(inout) :: CMO_VIR
    ! uses LUSIFC, unit number for SIRFC
#include "inftap.h"
    ! uses NBAST, NNBAST, NNBASX, N2BASX
#include "inforb.h"
    ! uses NODC, NODV
#include "cbione.h"
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
    call rd_sirifc( 'CMO', found, CMO%elms, f77_memory(get_f77_memory_next()), get_f77_memory_left() )
    if ( .not. found ) call QUIT( 'CMO not found on SIRIFC!' )
    !N if ( .not. restrict_scf ) CMO%elmsb = CMO%elms
    ! generates the occupied and virtual molecular orbitals
    CMO_OCC%elms(:,:NOCCT)   = CMO%elms(:,:NOCCT)
    CMO_OCC%elms(:,NOCCT+1:) = 0
    CMO_VIR%elms(:,:NOCCT)   = 0
    CMO_VIR%elms(:,NOCCT+1:) = CMO%elms(:,NOCCT+1:)
    ! closes SIRIFC
    if ( LUSIFC > 0 ) call GPCLOSE( LUSIFC, 'KEEP' )
  end subroutine
#endif /* ifdef PRG_DALTON */

end module
