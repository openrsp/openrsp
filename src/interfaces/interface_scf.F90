module interface_scf

   use matrix_defop
   use interface_f77_memory
   use interface_pcm

   implicit none

   public interface_scf_init
   public interface_scf_finalize

   public get_is_restricted_scf_calculation

   public interface_scf_get_s
   public interface_scf_get_d
   public interface_scf_get_h1
   public interface_scf_get_g

   private

!  if false the interface will refuse to be accessed
   logical :: is_initialized = .false.

!  non-allocatables
   logical :: is_restricted_scf_calculation

contains

   subroutine interface_scf_init()

#include "inforb.h"

      is_restricted_scf_calculation = (nasht == 0)

      is_initialized = .true.

   end subroutine

   subroutine interface_scf_finalize()

      is_initialized = .false.

   end subroutine

   subroutine check_if_interface_is_initialized()
      if (.not. is_initialized) then
         print *, 'error: you try to access interface_scf'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine

   logical function get_is_restricted_scf_calculation()
      call check_if_interface_is_initialized()
      get_is_restricted_scf_calculation = is_restricted_scf_calculation
   end function

  !> \brief gets the AO density matrix
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \return D contains the AO density matrix
  subroutine interface_scf_get_d(D)

    type(matrix), intent(inout) :: D

#ifdef PRG_DALTON
    ! uses NCMOT, NASHT, NNASHX, N2BASX
#include "inforb.h"
    ! uses LUSIFC
#include "inftap.h"
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
    strt_cmo = get_f77_memory_next()
    ! start of active part of one-electron density matrix (MO)
    strt_dv = strt_cmo + NCMOT
    if ( get_is_restricted_scf_calculation() ) then
      ! start of active part of one-electron density matrix (AO)
      strt_dvao = strt_dv + 1
      ! start of left workspace
      call set_f77_memory_next(strt_dvao + 1)
      ! only calculates DCAO
      GETDC = .true.
      GETDV = .false.
    else
      ! start of active part of one-electron density matrix (AO)
      strt_dvao = strt_dv + NNASHX
      ! start of left workspace
      call set_f77_memory_next(strt_dvao + N2BASX)
      ! we calculate DCAO and DVAO
      GETDC = .true.
      GETDV = .true.
    end if
    ! checks if the memory is enough
    if ( get_f77_memory_left() < 0 ) call STOPIT( 'DALTON_IFC', 'di_get_dens', get_f77_memory_next()-1, get_f77_memory_total() )
    ! opens SIRIFC
    if ( LUSIFC <= 0 ) &
      call GPOPEN( LUSIFC, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED', idummy, .false. )
    rewind( LUSIFC )
    ! reads the molecular orbital coefficients
    call DZERO( f77_memory(strt_cmo), NCMOT )
    call rd_sirifc( 'CMO', found, f77_memory(strt_cmo), f77_memory(get_f77_memory_next()), get_f77_memory_left() )
    if ( .not. found ) call QUIT( 'CMO not found on SIRIFC!' )
    ! reads active part of one-electron density matrix (MO)
    if ( GETDV ) then
      call DZERO( f77_memory(strt_dv), NNASHX )
      call rd_sirifc( 'DV', found, f77_memory(strt_dv), f77_memory(get_f77_memory_next()), get_f77_memory_left() )
      if ( .not. found ) call QUIT( 'DV not found on SIRIFC!' )
      call DZERO( f77_memory(strt_dvao), N2BASX )
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
    call FCKDEN( GETDC, GETDV, D%elms, f77_memory(strt_dvao), f77_memory(strt_cmo), &
                 f77_memory(strt_dv), f77_memory(get_f77_memory_next()), get_f77_memory_left() )
    ! sums DCAO and DVAO
    if ( GETDV ) &
      D%elms = D%elms + reshape( f77_memory(strt_dvao : strt_dvao+N2BASX-1), &
                                 (/D%nrow, D%ncol/) )
    ! clean
    call set_f77_memory_next(strt_cmo)
#endif /* ifdef PRG_DALTON */

  end subroutine



  !> \brief gets the two electron contribution (G) to Fock matrix
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param D contains the AO density matrix
  !> \return G contains the two electron contribution
  subroutine interface_scf_get_g(D, G)

    type(matrix), intent(in)    :: D
    type(matrix), intent(inout) :: G

#ifdef PRG_DALTON
    ! uses NBAST, NNBAST, NNBASX, N2BASX
#include "inforb.h"
    ! uses NODC, NODV
#include "cbione.h"
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
    work_ao_dens = get_f77_memory_next()
    call set_f77_memory_next(work_ao_dens + N2BASX)
    if ( get_f77_memory_left() < 0 ) call STOPIT( 'DALTON_IFC', 'di_get_gmat', get_f77_memory_next()-1, get_f77_memory_total() )
    ! sets the total density matrix
    !> \todo this may fail for unrestricted calculations
    call DCOPY( N2BASX, D%elms, 1, f77_memory(work_ao_dens), 1 )
    call DSCAL( N2BASX, two, f77_memory(work_ao_dens), 1 )
    ! outputs the total density matrix to check
    ! only one density matrix
    NDMAT = 1
    ISYMDM = 1
    !> \todo determines IFCTYP run-time
    IFCTYP = 3
    ! calculates two electron contribution by calling SIRFCK
    call SIRFCK( G%elms, f77_memory(work_ao_dens), NDMAT, &
                 ISYMDM, IFCTYP, .true., f77_memory(get_f77_memory_next()), get_f77_memory_left() )
    ! PCM two-electron contributions
    if ( get_is_pcm_calculation() ) then
      ! assigns the work memory
      work_pcm = get_f77_memory_next()
      work_pcm2 = work_pcm + NNBASX
      call set_f77_memory_next(work_pcm2 + N2BASX)
      if ( get_f77_memory_left() < 0 ) call STOPIT( 'DALTON_IFC', 'di_pcmfck2', get_f77_memory_next()-1, get_f77_memory_total() )
      call WAVPCM_2EL( f77_memory(work_pcm), f77_memory(work_ao_dens), f77_memory(get_f77_memory_next()), get_f77_memory_left() )
!      call pcm_ao_rsp_2elfock( f77_memory(work_pcm), f77_memory(work_ao_dens) )

      ! transforms to square matrix
      call DZERO( f77_memory(work_pcm2), N2BASX )
      call DSPTSI( NBAST, f77_memory(work_pcm), f77_memory(work_pcm2) )

      ! adds to G
      call DAXPY( N2BASX, 1D0, f77_memory(work_pcm2), 1, G%elms, 1 )
    end if
    !N if ( .not. restrict_scf ) G%elmsb = G%elms
    ! cleans
    call set_f77_memory_next(work_ao_dens)
#endif /* ifdef PRG_DALTON */

  end subroutine


  !> \brief gets the overlap and one-electron Hamiltonian matrices
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \return S contains the overlap matrix
  subroutine interface_scf_get_s(S)
    type(matrix), intent(inout) :: S

#ifdef PRG_DALTON
    ! uses NBAST, NNBASX, N2BASX
#include "inforb.h"
    ! IO units, use LUPROP for file AOPROPER
#include "inftap.h"
    ! starts of the overlap and one electron Hamiltonian matrices in work memory
    integer work_ovlp
    ! integer constants
    real(8), parameter :: one = 1.0D+00
    ! dummy stuff
    integer idummy
    ! external DALTON function finding the corresponding label
    logical FNDLAB
    ! reads overlap and one electron Hamiltonian matrices by calling RDONEL
    work_ovlp = get_f77_memory_next()
    call set_f77_memory_next(work_ovlp + NNBASX)

    if ( get_f77_memory_left() < 0 ) call STOPIT( 'DALTON_IFC', 'di_get_SH1', get_f77_memory_next()-1, get_f77_memory_total() )
    call RDONEL( 'OVERLAP', .true., f77_memory(work_ovlp), NNBASX )
    ! fills the data into matrices S and H1
    !N N2BASX = NBAST * NBAST
    if ( get_f77_memory_left() < 0 ) call STOPIT( 'DALTON_IFC', 'DSPTSI', get_f77_memory_next()+N2BASX-1, get_f77_memory_total() )
    ! gets S
    call DSPTSI( NBAST, f77_memory(work_ovlp), S%elms )
    ! clean
    call set_f77_memory_next(work_ovlp)
#endif /* ifdef PRG_DALTON */

  end subroutine


  !> \brief gets the overlap and one-electron Hamiltonian matrices
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \return H1 contains the one-electron Hamiltonian matrix
  subroutine interface_scf_get_h1(H1)
    type(matrix), intent(inout) :: H1

#ifdef PRG_DALTON
    ! uses NBAST, NNBASX, N2BASX
#include "inforb.h"
    ! IO units, use LUPROP for file AOPROPER
#include "inftap.h"
    ! starts of the overlap and one electron Hamiltonian matrices in work memory
    integer work_ham1
    ! PCM one-electron contributions
    integer work_pcm
    ! integer constants
    real(8), parameter :: one = 1.0D+00
    ! dummy stuff
    integer idummy
    ! external DALTON function finding the corresponding label
    logical FNDLAB
    ! reads overlap and one electron Hamiltonian matrices by calling RDONEL
    work_ham1 = get_f77_memory_next()
    call set_f77_memory_next(work_ham1 + NNBASX)

    if ( get_f77_memory_left() < 0 ) call STOPIT( 'DALTON_IFC', 'di_get_SH1', get_f77_memory_next()-1, get_f77_memory_total() )
    call RDONEL( 'ONEHAMIL', .true., f77_memory(work_ham1), NNBASX )
    ! PCM one-electron contributions
    if ( get_is_pcm_calculation() ) then
      work_pcm = get_f77_memory_next()
      call set_f77_memory_next(work_pcm + NNBASX)
      if ( get_f77_memory_left() < 0 ) call STOPIT( 'DALTON_IFC', 'di_get_SH1', get_f77_memory_next()-1, get_f77_memory_total() )
#ifdef USE_WAVPCM
      call pcm_ao_rsp_1elfock( f77_memory(work_pcm) )
      ! adds to one-electron Hamiltonian
      f77_memory( work_ham1 : work_ham1 + NNBASX - 1 ) &
                    = f77_memory( work_ham1 : work_ham1 + NNBASX - 1 ) &
                    + f77_memory( work_pcm  : work_pcm  + NNBASX - 1 )
#endif
      ! cleans
      call set_f77_memory_next(work_pcm)
    end if
    ! fills the data into matrices S and H1
    !N N2BASX = NBAST * NBAST
    if ( get_f77_memory_left() < 0 ) call STOPIT( 'DALTON_IFC', 'DSPTSI', get_f77_memory_next()+N2BASX-1, get_f77_memory_total() )
    ! gets S
    ! gets H1
    call DSPTSI( NBAST, f77_memory(work_ham1), H1%elms )
    ! clean
    call set_f77_memory_next(work_ham1)
#endif /* ifdef PRG_DALTON */

  end subroutine

end module
