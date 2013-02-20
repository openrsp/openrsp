module interface_rsp_solver

   use matrix_defop
   use matrix_lowlevel, only: mat_init
   use interface_scf
   use interface_f77_memory
   use interface_io

   implicit none

   public rsp_mosolver_init
   public rsp_mosolver_splash
#ifndef VAR_LSDALTON
   public rsp_solver_exec
#endif
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

    call mat_init(solver_CMO, NBAST, NBAST)

    solver_CMO_OCC = 0*solver_CMO
    call mat_ensure_alloc(solver_CMO_OCC, only_alloc=.true.)
    solver_CMO_VIR = 0*solver_CMO
    call mat_ensure_alloc(solver_CMO_VIR, only_alloc=.true.)
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

#ifdef PRG_DALTON
  !> \brief calls MO response solver in DALTON
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param RHS is the right hand side vectors (property gradients?)
  !> \param eigval contains the frequencies
  !> \return eigvec contains the solution vectors (AO)
  subroutine rsp_solver_exec(RHS,    &
                               eigval, &
                               eigvec)

    type(matrix), intent(in)    :: RHS(*)
    real(8),      intent(in)    :: eigval(*)
    type(matrix), intent(inout) :: eigvec(*)

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

    call mat_init(RHS_MO, NORBT, NORBT)
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

      ! MaR: Next line replaced by mat_init call to avoid strange nonallocation error msg
      ! mo_eigvec(ISOL) = mat_alloc_like(RHS_MO)
      ! ASSUMES CLOSED SHELL
      call mat_init(mo_eigvec(ISOL), RHS_MO%nrow, RHS_MO%ncol)

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
      eigvec(ISOL) = - solver_CMO_OCC*( mo_eigvec(ISOL)*trans( solver_CMO_VIR ) ) &
                     - solver_CMO_VIR*( mo_eigvec(ISOL)*trans( solver_CMO_OCC ) )

      ! radovan: strange factor, was previously further up
      eigvec(ISOL) = -2.0d0*eigvec(ISOL)
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

  end subroutine
#endif /* ifdef PRG_DALTON */

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
    !N if ( .not. restrict_scf ) CMO%elms = CMO%elms
    ! generates the occupied and virtual molecular orbitals
    CMO_OCC%elms(:,:NOCCT,1)   = CMO%elms(:,:NOCCT,1)
    CMO_OCC%elms(:,NOCCT+1:,1) = 0
    CMO_VIR%elms(:,:NOCCT,1)   = 0
    CMO_VIR%elms(:,NOCCT+1:,1) = CMO%elms(:,NOCCT+1:,1)
    ! closes SIRIFC
    if ( LUSIFC > 0 ) call GPCLOSE( LUSIFC, 'KEEP' )
  end subroutine
#endif /* ifdef PRG_DALTON */


#ifdef PRG_DIRAC
  subroutine set_orbrot_indices(irep,       &
                                          length,     &
                                          from_block, &
                                          to_block,   &
                                          orbital_rotation_indices)

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: irep
    integer, intent(in)  :: length
    integer, intent(in)  :: from_block
    integer, intent(in)  :: to_block
    integer, intent(out) :: orbital_rotation_indices(2, length)
!   ----------------------------------------------------------------------------
    integer, parameter   :: ns = 1
    integer, parameter   :: pi = 2
    integer, parameter   :: ps = 3

    integer, parameter   :: g  = 1
    integer, parameter   :: u  = 2

    logical              :: gu_combine(2, 2)

    integer              :: nr_mo(2, 3)       = 0
    integer              :: range_mo(2, 3, 2) = 0

    integer              :: i, s, k
    integer              :: m, m1, m2
    integer              :: gu, gu1, gu2
    integer              :: block, block1, block2

    logical              :: debug_me = .false.
!   ----------------------------------------------------------------------------

#include "dcbbas.h"
#include "dcborb.h"
#include "dcbxpr.h"
#include "dcbxrs.h"
#include "dgroup.h"

    nr_mo(g, ns) = npsh(g)
    nr_mo(g, pi) = nish(g)
    nr_mo(g, ps) = nesh(g) - nish(g)

    nr_mo(u, ns) = npsh(u)
    nr_mo(u, pi) = nish(u)
    nr_mo(u, ps) = nesh(u) - nish(u)

    k = 0
    do gu = 1, 2
      do block = 1, 3
        do m = 1, nr_mo(gu, block)
          k = k + 1
          if (k > range_mo(gu, block, 1)) range_mo(gu, block, 1) = k
          range_mo(gu, block, 2) = k
        end do
      end do
    end do

    gu_combine = .false.
    if (jbtof(irep, 1) == 2) then
!     ungerade perturbation
      gu_combine(1, 2) = .true.
      gu_combine(2, 1) = .true.
    else
!     gerade perturbation
      gu_combine(1, 1) = .true.
      gu_combine(2, 2) = .true.
    end if

    k = 0
    i = 0
    do gu1 = 1, 2
      do block1 = 1, 3
        do m1 = 1, nr_mo(gu1, block1)
          i = i + 1
          s = 0
          do gu2 = 1, 2
            do block2 = 1, 3
              if (gu_combine(gu1, gu2)         &
                  .and. (block1 == from_block) &
                  .and. (block2 == to_block)) then
                do m2 = 1, nr_mo(gu2, block2)
                  k = k + 1
                  s = s + 1
                  orbital_rotation_indices(1, k) = i
                  orbital_rotation_indices(2, k) = s
                end do
              else
                s = s + nr_mo(gu2, block2)
              end if
            end do
          end do
        end do
      end do
    end do

    if (debug_me) then
      call header('debug orbital_rotation_indices', -1)
      do k = 1, length
        write(*, *) k, orbital_rotation_indices(1, k), &
                 '->', orbital_rotation_indices(2, k)
      end do
    end if

  end subroutine

  subroutine rsp_solver_exec(RHS,    &
                               eigval, &
                               Dp)

    type(matrix), intent(in)    :: RHS(*)
    real(8),      intent(in)    :: eigval(*)
    type(matrix), intent(inout) :: Dp(*)

!   1. RHS (AO) -> RHS (MO)
!   2. RHS (MO) -> property gradient
!   3. solve response equation, in:  property gradient
!                               out: response vector
!   4. response vector -> Wp (MO)
!   5. Wp (MO)         -> Dp (pert dens mat; AO)

!   ----------------------------------------------------------------------------
    type(matrix)                :: RHS_mo
    type(matrix)                :: Wp
    type(matrix)                :: C
    type(matrix)                :: Cig
    type(matrix)                :: Ciu
    type(matrix)                :: Csg
    type(matrix)                :: Csu

    real(8),      allocatable   :: mo_coef(:)

    integer,      allocatable   :: ibtyp(:)
    integer,      allocatable   :: ibtyp_pointer_pp(:)
    integer,      allocatable   :: ibtyp_pointer_pn(:)

    real(8),      allocatable   :: eigvec(:)
    real(8),      allocatable   :: convergence(:)

    real(8),      allocatable   :: response_vector_pph(:, :, :)
    real(8),      allocatable   :: response_vector_ppa(:, :, :)
    real(8),      allocatable   :: response_vector_pnh(:, :, :)
    real(8),      allocatable   :: response_vector_pna(:, :, :)

    real(8),      allocatable   :: prop_gradient_pp(:, :)
    real(8),      allocatable   :: prop_gradient_pn(:, :)

    integer,      allocatable   :: kappa_pp(:, :)
    integer,      allocatable   :: kappa_pn(:, :)

    integer                     :: length_pp
    integer                     :: length_pn

    integer                     :: iz
    integer                     :: k
    integer                     :: i
    integer                     :: s
    integer                     :: nr_freq
    integer                     :: kfree

    logical                     :: debug_me = .false.

    logical                     :: include_pp_rotations = .true.
    logical                     :: include_pn_rotations = .true.
!   ----------------------------------------------------------------------------

#include "dcbbas.h"
#include "dcborb.h"
#include "dcbxpr.h"
#include "dcbxrs.h"
#include "dgroup.h"

!   fixme: nssh not always properly set, change to nesh minus nish

    if (jbtof(RHS(1)%pg_sym-1, 1) == 2) then
!     ungerade perturbation
      length_pp = nish(1)*nssh(2) + nish(2)*nssh(1)
      length_pn = nish(1)*npsh(2) + nish(2)*npsh(1)
    else
!     gerade perturbation
      length_pp = nish(1)*nssh(1) + nish(2)*nssh(2)
      length_pn = nish(1)*npsh(1) + nish(2)*npsh(2)
    end if

    nzconf  = 0
    nzxope  = 0
    nzxopp  = 0

    if (include_pp_rotations) nzxope = length_pp
    if (include_pn_rotations) nzxopp = length_pn

    nzxopt  = nzxope + nzxopp

    nzconfq = nzconf*nz
    nzxopeq = nzxope*nz
    nzxoppq = nzxopp*nz
    nzxoptq = nzxopt*nz

    nzvar   = nzxopt
    nzvarq  = nzxoptq

    allocate(kappa_pp(2, length_pp))
    if (length_pp > 0) then
       call set_orbrot_indices(RHS(1)%pg_sym-1,  &
                                         length_pp, &
                                         2,         &
                                         3,         &
                                         kappa_pp)
    end if

    allocate(kappa_pn(2, length_pn))
    if (length_pn > 0) then
       call set_orbrot_indices(RHS(1)%pg_sym-1,  &
                                         length_pn, &
                                         2,         &
                                         1,         &
                                         kappa_pn)
    end if

!   fixme stop at spinfree
!   see old interface to see how it is done

       allocate(mo_coef(ncmotq))
       call read_mo_coef(mo_coef)
       call mat_init(C, ntbas(0), norbt)
       call get_C(C, mo_coef, i=1.0d0, s=1.0d0, g=1.0d0, u=1.0d0)
       deallocate(mo_coef)
       RHS_mo = trans(C)*(RHS(1)*C)
       C = 0

    RHS_mo%ih_sym = RHS(1)%ih_sym
    RHS_mo%pg_sym = RHS(1)%pg_sym


!   get positive -> positive property gradient
!   ==========================================

    if (include_pp_rotations) then

       allocate(prop_gradient_pp(length_pp, nz))

       do iz = 1, nz
          do k = 1, length_pp
             i = kappa_pp(1, k)
             s = kappa_pp(2, k)
             prop_gradient_pp(k, iz) = -2.0d0*RHS_mo%elms(s, i, iz)
          end do
       end do

       if (debug_me) then
          call header('debug prop_gradient_pp', -1)
          call print_q_vector(prop_gradient_pp, length_pp)
       end if

    end if


!   get positive -> negative property gradient
!   ==========================================

    if (include_pn_rotations) then

       allocate(prop_gradient_pn(length_pn, nz))

       do iz = 1, nz
          do k = 1, length_pn
             i = kappa_pn(1, k)
             s = kappa_pn(2, k)
             prop_gradient_pn(k, iz) = -2.0d0*RHS_mo%elms(s, i, iz)
          end do
       end do

       if (debug_me) then
          call header('debug prop_gradient_pn', -1)
          call print_q_vector(prop_gradient_pn, length_pn)
       end if

    end if


!   inherit shape
    Wp     = RHS_mo
    RHS_mo = 0

    call set_orbital_rotation_indices(                             &
                                      include_pp_rotations,        &
                                      include_pn_rotations,        &
                                      length_pp,                   &
                                      length_pn,                   &
                                      kappa_pp, &
                                      kappa_pn  &
                                     )


!   parameters
!   ==========

    call setrsp()

!   .URKBAL response with frequency=0.0 but static=.false. converges much worse
!   than frequency=0.0 with static=.true.
!   it is important to correctly set static
    static = .true.
    if (dabs(eigval(1)) > tiny(0.0d0)) then
       static = .false.
    end if

    lineq     = .true.
    lsvcfg(1) = .true.
    lsvcfg(2) = .true.
    tknorm    = .true.
    diaghe    = .true.
    iprxrs    = 0
!   thcxrs    = openrsp_cfg_threshold_response
    thcxrs    = 1.0d-8 !fixme hardcoded
    resfac    = 1.0d3
    maxitr    = 150
    nr_freq   = 1
    maxsim    = -1
    nredm     = 400
    n2redm    = nredm*nredm
    loffty    = 0
    uncoup    = .false.
    imfreq    = .false.
    triplet   = .false.
    cnvint(1) = 1.0d20
    cnvint(2) = 1.0d20
    itrint(1) = 1
    itrint(2) = 1
!   intdef    = integral_flag
    intdef    = 3 !fixme hardcoded
    sternh    = .false.
    sternc    = .false.
    e2chek    = .false.

!   nr_freq   = size(freq)
    jsymop    = RHS(1)%pg_sym
!   jtimop    = RHS%tr_sym
    jtimop    = RHS(1)%ih_sym
    jopsy     = jbtof(RHS(1)%pg_sym-1, 1)
    nfreq     = nr_freq
    nexsim    = nr_freq
    nexstv    = nr_freq
    nexcnv    = nr_freq
    ncred     = 0
    nered     = 0
    npred     = 0
    nzred     = ncred + nered + npred


!   solve response equation in the reduced space spanned by the trial vectors
!   =========================================================================

    allocate(ibtyp(            2*nredm))
    allocate(ibtyp_pointer_pp( nredm))
    allocate(ibtyp_pointer_pn( nredm))
    allocate(convergence(      nr_freq))
    allocate(eigvec(           nredm*nr_freq))

    kfree = 1
    call xrsctl(                                   &
                (/0.0d0/),                         &
                prop_gradient_pp,                  &
                prop_gradient_pn,                  &
                ibtyp,                             &
                (/0/),                             &
                ibtyp_pointer_pp,                  &
                ibtyp_pointer_pn,                  &
                convergence,                       &
                eigval,                            &
                eigvec,                            &
                f77_memory(get_f77_memory_next()), &
                kfree,                             &
                get_f77_memory_left()              &
               )

    deallocate(convergence)


!   construct response vector
!   =========================

    if (include_pp_rotations) then

      allocate(response_vector_pph(length_pp, nz, nr_freq))
      allocate(response_vector_ppa(length_pp, nz, nr_freq))

      response_vector_pph = 0.0d0
      response_vector_ppa = 0.0d0

      call construct_response_vector('pp',                &
                                     1,                   &
                                     response_vector_pph, &
                                     length_pp,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pp,    &
                                     eigvec)

      call construct_response_vector('pp',                &
                                     -1,                  &
                                     response_vector_ppa, &
                                     length_pp,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pp,    &
                                     eigvec)

      if (debug_me) then
        call header('debug response_vector_pph', -1)
        call prbvec(6, response_vector_pph, nr_freq, length_pp)
        call header('debug response_vector_ppa', -1)
        call prbvec(6, response_vector_ppa, nr_freq, length_pp)
      end if

    end if

    if (include_pn_rotations) then

      allocate(response_vector_pnh(length_pn, nz, nr_freq))
      allocate(response_vector_pna(length_pn, nz, nr_freq))

      response_vector_pnh = 0.0d0
      response_vector_pna = 0.0d0

      call construct_response_vector('pn',                &
                                     1,                   &
                                     response_vector_pnh, &
                                     length_pn,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pn,    &
                                     eigvec)

      call construct_response_vector('pn',                &
                                     -1,                  &
                                     response_vector_pna, &
                                     length_pn,           &
                                     nr_freq,             &
                                     ibtyp,               &
                                     ibtyp_pointer_pn,    &
                                     eigvec)

      if (debug_me) then
        call header('debug response_vector_pnh', -1)
        call prbvec(6, response_vector_pnh, nr_freq, length_pn)
        call header('debug response_vector_pna', -1)
        call prbvec(6, response_vector_pna, nr_freq, length_pn)
      end if

    end if

    deallocate(ibtyp)
    deallocate(ibtyp_pointer_pp)
    deallocate(ibtyp_pointer_pn)
    deallocate(eigvec)


!   scatter response vectors
!   ========================

    Wp%elms = 0.0d0

    if (include_pp_rotations) then
      call scatter_vector(length_pp,                   &
                          kappa_pp, &
                          1.0d0,                       &
                          response_vector_pph,         &
                          Wp%elms,                     &
                          Wp%pg_sym-1)
      call scatter_vector(length_pp,                   &
                          kappa_pp, &
                         -1.0d0,                       &
                          response_vector_ppa,         &
                          Wp%elms,                     &
                          Wp%pg_sym-1)

      deallocate(kappa_pp)
      deallocate(response_vector_pph)
      deallocate(response_vector_ppa)
    end if

    if (include_pn_rotations) then
      call scatter_vector(length_pn,                   &
                          kappa_pn, &
                          1.0d0,                       &
                          response_vector_pnh,         &
                          Wp%elms,                     &
                          Wp%pg_sym-1)
      call scatter_vector(length_pn,                   &
                          kappa_pn, &
                         -1.0d0,                       &
                          response_vector_pna,         &
                          Wp%elms,                     &
                          Wp%pg_sym-1)

      deallocate(kappa_pn)
      deallocate(response_vector_pnh)
      deallocate(response_vector_pna)
    end if


!   get coefficients
!   ================

    allocate(mo_coef(n2bbasxq))
    call read_mo_coef(mo_coef)
    call mat_init(Cig, ntbas(0), norbt)
    call mat_init(Csg, ntbas(0), norbt)
    call get_C(Cig, mo_coef, i=1.0d0, s=0.0d0, g=1.0d0, u=0.0d0)
    call get_C(Csg, mo_coef, i=0.0d0, s=1.0d0, g=1.0d0, u=0.0d0)
    if (nfsym == 2) then
      call get_C(Ciu, mo_coef, i=1.0d0, s=0.0d0, g=0.0d0, u=1.0d0)
      call get_C(Csu, mo_coef, i=0.0d0, s=1.0d0, g=0.0d0, u=1.0d0)
    end if
    deallocate(mo_coef)


!   construct perturbed AO density matrix
!   =====================================

    if (nfsym == 2) then
      if (jbtof(Wp%pg_sym-1, 1) == 2) then
!       ungerade perturbation
        Dp(1) = (Cig*(Wp*trans(Csu))) &
           + (Ciu*(Wp*trans(Csg))) &
           - (Csg*(Wp*trans(Ciu))) &
           - (Csu*(Wp*trans(Cig)))
      else
!       gerade perturbation
        Dp(1) = (Cig*(Wp*trans(Csg))) &
           + (Ciu*(Wp*trans(Csu))) &
           - (Csg*(Wp*trans(Cig))) &
           - (Csu*(Wp*trans(Ciu)))
      end if
      Ciu = 0
      Csu = 0
    else
      Dp(1) = (Cig*(Wp*trans(Csg))) &
         - (Csg*(Wp*trans(Cig)))
    end if
    Cig = 0
    Csg = 0

    Dp(1)%ih_sym = Wp%ih_sym
    Dp(1)%pg_sym = Wp%pg_sym

    Wp = 0

  end subroutine

  subroutine construct_response_vector(vector_type,     &
                                       ih,              &
                                       response_vector, &
                                       length,          &
                                       nr_freq,         &
                                       ibtyp,           &
                                       ibtyp_pointer,   &
                                       eigvec)

!   ----------------------------------------------------------------------------
    character(*), intent(in)  :: vector_type
    integer,      intent(in)  :: ih
    real(8),      intent(out) :: response_vector(*)
    integer,      intent(in)  :: length
    integer,      intent(in)  :: nr_freq
    integer,      intent(in)  :: ibtyp(*)
    integer,      intent(in)  :: ibtyp_pointer(*)
    real(8),      intent(in)  :: eigvec(*)
!   ----------------------------------------------------------------------------
    integer,      allocatable :: ivecs(:)
    real(8),      allocatable :: buffer(:)
    character(6)              :: file_name
    integer                   :: file_unit, nbtyp
!   ----------------------------------------------------------------------------

#include "dcbbas.h"
#include "dcborb.h"
#include "dcbxpr.h"
#include "dcbxrs.h"
#include "dgroup.h"

    select case (vector_type)
      case ('pp')
        nbtyp     = 1
        file_unit = luboe
        file_name = 'PAMBOE'
      case ('pn')
        nbtyp     = 2
        file_unit = lubop
        file_name = 'PAMBOP'
      case default
        call quit('error in construct_response_vector: unknown vector type')
    end select

    allocate(ivecs(  nr_freq))
    allocate(buffer( length*nz))

    open(file_unit,              &
         file   = file_name,     &
         form   = 'unformatted', &
         access = 'direct',      &
         recl   = 8*length*nz,   &
         status = 'old')

    call xrsxv1(ih,              &
                nbtyp,           &
                response_vector, &
                eigvec,          &
                nr_freq,         &
                ibtyp,           &
                ibtyp_pointer,   &
                ivecs,           &
                buffer)

    close(file_unit, status = 'keep')

    deallocate(ivecs)
    deallocate(buffer)

  end subroutine

  subroutine print_vector(vector, vector_dimension)

    real(8), intent(in) :: vector(:)
    integer,       intent(in) :: vector_dimension

    call prqmat(vector,           &
                vector_dimension, &
                1,                &
                vector_dimension, &
                1,                &
                1,                &
                (/1, 2, 3, 4/),   &
                6)

  end subroutine

  subroutine print_q_vector(vector, vector_dimension)

    real(8), intent(in) :: vector(*)
    integer,       intent(in) :: vector_dimension

#include "dgroup.h"

    call prqmat(vector,           &
                vector_dimension, &
                1,                &
                vector_dimension, &
                1,                &
                nz,               &
                (/1, 2, 3, 4/),   &
                6)

  end subroutine

  subroutine scatter_vector(length,                   &
                            orbital_rotation_indices, &
                            h,                        &
                            response_vector,          &
                            matrix,                   &
                            irep)

#include "dgroup.h"
#include "dcborb.h"

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: length
    integer, intent(in)  :: orbital_rotation_indices(2, length)
    real(8), intent(in)  :: h
    real(8), intent(in)  :: response_vector(length, *)
    real(8), intent(out) :: matrix(norbt, norbt, *)
    integer, intent(in)  :: irep
!   ----------------------------------------------------------------------------
    integer              :: i, s, is, iz
    real(8)              :: f
!   ----------------------------------------------------------------------------

    do is = 1, length

      i = orbital_rotation_indices(1, is)
      s = orbital_rotation_indices(2, is)

      do iz = 1, nz

        if (ipqtoq(iz, irep) > 1) then
          f = -1.0
        else
          f =  1.0
        end if

!              row
!              ¦
!              ¦  column
!              ¦  ¦
        matrix(s, i, iz) = matrix(s, i, iz) &
                         +     response_vector(is, iz)
        matrix(i, s, iz) = matrix(i, s, iz) &
                         - f*h*response_vector(is, iz)
      end do
    end do

  end subroutine
#endif /* ifdef PRG_DIRAC */

end module
