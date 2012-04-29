! radovan: do not add anything here
!          this module is about to retire

module dalton_ifc

  use matrix_defop   !type matrix with operators
  use interface_f77_memory
  use interface_molecule
  use interface_io
  use interface_pcm
  use interface_scf

  implicit none

  public di_get_dens
  public di_get_gmat
  public di_twofck

  public NUCLEI_ifc
  public DIPNUC_ifc
  public QDRNUC_ifc
  public GRADNN_ifc
  public HESSNN_ifc
  public DPGNUC_ifc
  public AATNUC_ifc
  public VIBCTL_ifc

  contains



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
  end subroutine



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
  end subroutine



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
  end subroutine






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
  end subroutine


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
  end subroutine


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
    real(8), intent(out) :: Q(6)
    real(8), parameter :: zero = 0.0D+00
    call NUCQDR( CORD(:,1:NUCDEP), (/zero,zero,zero/), get_print_unit(), 0 )
    Q = (/QDRNUC(1,1),QDRNUC(1:2,2),QDRNUC(1:3,3)/)
  end subroutine


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
  end subroutine


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
  end subroutine


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
  end subroutine


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
  end subroutine


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
    call VIBCTL( f77_memory(get_f77_memory_next()), get_f77_memory_left() )
  end subroutine


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
    !radovan: fixme ugly dependency on basis_set
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
