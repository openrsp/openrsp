module interface_scf

   use openrsp_const
   use matrix_defop
   use interface_f77_memory
   use interface_pcm
   use interface_io
   use interface_basis
   use interface_dirac_gen1int
   use eri_contractions, only: ctr_arg
   use eri_basis_loops,  only: unopt_geodiff_loop

   implicit none

   public interface_scf_init
   public interface_scf_finalize

   public get_is_restricted_scf_calculation

   public interface_scf_get_s
   public interface_scf_get_d
   public interface_scf_get_h1
   public interface_scf_get_g
#ifdef PRG_DIRAC
   public get_C
   public read_mo_coef
#endif

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
    call FCKDEN( GETDC, GETDV, D%elms_alpha, f77_memory(strt_dvao), f77_memory(strt_cmo), &
                 f77_memory(strt_dv), f77_memory(get_f77_memory_next()), get_f77_memory_left() )
    ! sums DCAO and DVAO
    if ( GETDV ) &
      D%elms_alpha(:, :, 1) = D%elms_alpha(:, :, 1) + reshape( f77_memory(strt_dvao : strt_dvao+N2BASX-1), &
                                 (/D%nrow, D%ncol/) )
    ! clean
    call set_f77_memory_next(strt_cmo)
#endif /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
      real(8), allocatable :: mo_coef(:)
      type(matrix)         :: C_i

! uses ntbas(0)
#include "dcbbas.h"
! uses nz
#include "dgroup.h"
! uses ncmotq and norbt
#include "dcborb.h"

      if (nz /= 4) then
         print *, 'adapt interface_scf_get_d for nz /= 4'
         stop 1
      end if

      allocate(mo_coef(ncmotq))

      call mat_init(C_i, nrow=ntbas(0), ncol=norbt, closed_shell=.true.)

      call read_mo_coef(mo_coef)

      call get_C(C_i, mo_coef, i=1.0d0, s=0.0d0, g=1.0d0, u=1.0d0)

      deallocate(mo_coef)

      D = C_i*(trps(C_i))
      D = 2.0d0*D

      C_i     = 0
#endif /* ifdef PRG_DIRAC */

   end subroutine



  !> \brief gets the two electron contribution (G) to Fock matrix
  !> \author Bin Gao
  !> \date 2009-12-08
  !> \param D contains the AO density matrix
  !> \return G contains the two electron contribution
  subroutine interface_scf_get_g(D, G)

    use interface_interest

    type(matrix), intent(in),    target :: D
    type(matrix), intent(inout), target :: G

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
    call DCOPY( N2BASX, D%elms_alpha, 1, f77_memory(work_ao_dens), 1 )
    call DSCAL( N2BASX, two, f77_memory(work_ao_dens), 1 )
    ! outputs the total density matrix to check
    ! only one density matrix
    NDMAT = 1
    ISYMDM = 1
    !> \todo determines IFCTYP run-time
    IFCTYP = 3
    ! calculates two electron contribution by calling SIRFCK
    call SIRFCK( G%elms_alpha, f77_memory(work_ao_dens), NDMAT, &
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
      call DAXPY( N2BASX, 1D0, f77_memory(work_pcm2), 1, G%elms_alpha, 1 )
    end if
    !N if ( .not. restrict_scf ) G%elms_alpha = G%elms_alpha
    ! cleans
    call set_f77_memory_next(work_ao_dens)
#endif /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
    G%elms_alpha = 0.0d0
    call interest_get_int(D%nrow, D%elms_alpha, G%elms_alpha, 0)
#endif /* ifdef PRG_DIRAC */

  end subroutine


   !> \brief gets the overlap
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
    call DSPTSI( NBAST, f77_memory(work_ovlp), S%elms_alpha )
    ! clean
    call set_f77_memory_next(work_ovlp)
#endif /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
      call get_1el_integrals(                            &
                             M=(/S/),                    &
                             prop_name="INT_OVERLAP",    &
                             num_ints=1,                 &
                             order_mom=0,                &
                             order_elec=1,               &
                             order_geo_total=0,          &
                             max_num_cent=0,             &
                             blocks=(/1, 1, 2, 2/),      &
                             print_unit=get_print_unit() &
                            )
!     call mat_print(S)
#endif /* ifdef PRG_DIRAC */

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
    call DSPTSI( NBAST, f77_memory(work_ham1), H1%elms_alpha )
    ! clean
    call set_f77_memory_next(work_ham1)
#endif /* ifdef PRG_DALTON */

#ifdef PRG_DIRAC
      type(matrix)       :: T, TX, TY, TZ

! uses nz
#include "dgroup.h"

      if (nz /= 4) then
         print *, 'adapt interface_scf_get_h1 for nz /= 4'
         stop 1
      end if

      T  = mat_alloc_like(H1)
      TX = mat_alloc_like(H1)
      TY = mat_alloc_like(H1)
      TZ = mat_alloc_like(H1)

      call get_1el_integrals(                                &
                             M=(/TX, TY, TZ/),               &
                             prop_name="INT_CART_MULTIPOLE", &
                             num_ints=3,                     &
                             order_mom=0,                    &
                             order_elec=1,                   &
                             order_geo_total=0,              &
                             max_num_cent=0,                 &
                             blocks=(/1, 2, 2, 1/),          &
                             print_unit=get_print_unit()     &
                            )

      T%elms_alpha = 0.0d0

      call dcopy(H1%nrow*H1%ncol, TX%elms_alpha, 1, T%elms_alpha(1, 1, 4), 1)
      call dcopy(H1%nrow*H1%ncol, TY%elms_alpha, 1, T%elms_alpha(1, 1, 3), 1)
      call dcopy(H1%nrow*H1%ncol, TZ%elms_alpha, 1, T%elms_alpha(1, 1, 2), 1)

      H1%elms_alpha = 0.0d0

      ! kinetic energy = c (\vec \alpha \cdot \vec p)
      H1 = H1 - openrsp_const_speed_of_light*T

      call get_1el_integrals(                            &
                             M=(/T/),                    &
                             prop_name="INT_OVERLAP",    &
                             num_ints=1,                 &
                             order_mom=0,                &
                             order_elec=1,               &
                             order_geo_total=0,          &
                             max_num_cent=0,             &
                             blocks=(/2, 2/),            &
                             print_unit=get_print_unit() &
                            )

      ! shifted rest mass energy = \beta' c*c
      H1 = H1 - 2.0d0*(openrsp_const_speed_of_light**2.0d0)*T

      call get_1el_integrals(                            &
                             M=(/T/),                    &
                             prop_name="INT_POT_ENERGY", &
                             num_ints=1,                 &
                             order_mom=0,                &
                             order_elec=1,               &
                             order_geo_total=0,          &
                             max_num_cent=0,             &
                             blocks=(/1, 1, 2, 2/),      &
                             print_unit=get_print_unit() &
                            )

      ! electron-nuclei attraction energy
      H1 = H1 + T

      T  = 0
      TX = 0
      TY = 0
      TZ = 0

#endif /* ifdef PRG_DIRAC */

  end subroutine

#ifdef PRG_DIRAC
   subroutine read_mo_coef(mo_coef)

!     --------------------------------------------------------------------------
      real(8)            :: mo_coef(*)
!     --------------------------------------------------------------------------
      integer, parameter :: file_unit = 66
      real(8)            :: energy
!     --------------------------------------------------------------------------

      call opnfil(file_unit, 'DFCOEF', 'OLD', 'read_mo_coef')

      call reacmo_no_work(file_unit, &
                          'DFCOEF',  &
                          mo_coef,   &
                          (/0.0d0/), &
                          (/0/),     &
                          energy,    &
                          2)

      close(file_unit, status = 'keep')

   end subroutine

   subroutine get_C(C, mo_coef, i, s, g, u)

!     --------------------------------------------------------------------------
      type(matrix)        :: C
      real(8), intent(in) :: mo_coef(*)
      real(8), intent(in) :: i, s, g, u
!     --------------------------------------------------------------------------
      integer             :: nr_g_mo_ns
      integer             :: nr_g_mo_pi
      integer             :: nr_g_mo_ps
      integer             :: nr_u_mo_ns
      integer             :: nr_u_mo_pi
      integer             :: nr_u_mo_ps
      integer             :: nr_g_mo
      integer             :: nr_u_mo
      integer             :: nr_g_ao
      integer             :: nr_u_ao
      integer             :: k, l, m, n, ir, iw1, iw2, iz
!     --------------------------------------------------------------------------

#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"

      nr_g_mo_ns = npsh(1)
      nr_g_mo_pi = nish(1)
      nr_g_mo_ps = nesh(1) - nr_g_mo_pi
      nr_g_mo    = nr_g_mo_ns + nr_g_mo_pi + nr_g_mo_ps
      nr_u_mo_ns = npsh(2)
      nr_u_mo_pi = nish(2)
      nr_u_mo_ps = nesh(2) - nr_u_mo_pi
      nr_u_mo    = nr_u_mo_ns + nr_u_mo_pi + nr_u_mo_ps
      nr_g_ao    = nfbas(1, 0)
      nr_u_ao    = nfbas(2, 0)

      ir = 0

      do iz = 1, nz
         k = 0
         do m = 1, nr_g_mo_ns
            k = k + 1
            do n = 1, nr_g_ao
               ir = ir + 1
               C%elms_alpha(n, k, iz) = mo_coef(ir)*g*s
            end do
         end do
         do m = 1, nr_g_mo_pi
            k = k + 1
            do n = 1, nr_g_ao
               ir = ir + 1
               C%elms_alpha(n, k, iz) = mo_coef(ir)*g*i
            end do
         end do
         do m = 1, nr_g_mo_ps
            k = k + 1
            do n = 1, nr_g_ao
               ir = ir + 1
               C%elms_alpha(n, k, iz) = mo_coef(ir)*g*s
            end do
         end do
      end do

      if (nfsym == 1) return

      do iz = 1, nz
         k = 0
         do m = 1, nr_u_mo_ns
            k = k + 1
            do n = nr_g_ao + 1, nr_g_ao + nr_u_ao
               ir = ir + 1
               C%elms_alpha(n, k, iz) = mo_coef(ir)*u*s
            end do
         end do
         do m = 1, nr_u_mo_pi
            k = k + 1
            do n = nr_g_ao + 1, nr_g_ao + nr_u_ao
               ir = ir + 1
               C%elms_alpha(n, k, iz) = mo_coef(ir)*u*i
            end do
         end do
         do m = 1, nr_u_mo_ps
            k = k + 1
            do n = nr_g_ao + 1, nr_g_ao + nr_u_ao
               ir = ir + 1
               C%elms_alpha(n, k, iz) = mo_coef(ir)*u*s
            end do
         end do
      end do

   end subroutine
#endif /* ifdef PRG_DIRAC */

end module
