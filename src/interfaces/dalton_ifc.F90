! radovan: do not add anything here
!          this module is about to retire

module dalton_ifc

  use interface_f77_memory
  use interface_io

  implicit none

  public vibctl_ifc
  public shells_nuclei_displace

  private

  contains

  !> \brief gets the
  !> \author Bin Gao
  !> \date 2009-12-12
  !> \param
  !> \return
  !> \note modified on linsca/linears/VIBCTL_interface.F
  subroutine VIBCTL_ifc( nc, w, ALPHA, GPRIME, THETA, DIPGRAD, &
                         dALPHAdR, dGPRIMEdR, dTHETAdR )
    implicit integer (i,m-n)
#include "implicit.h"
#include "mxcent.h"
#include "cbilnr.h"
#include "cbivib.h"
#include "abainf.h"
#include "moldip.h"
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
#ifdef PRG_DIRAC
    print *, 'fix vibini and vibctl call'
    stop 1
#else
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
#endif
  end subroutine


  !> Count the number of contracted Gaussian-type orbital shells
  !> \param ncgto, and number of exponents and contraction coefficents
  !> \param nectr, so that memory for data structures can be allocated
  subroutine SHELLS_find_sizes(ncgto, nectr)
    implicit integer (i,m-n)
#include "implicit.h"
    integer, intent(out) :: ncgto, nectr
    ! need MXSHEL
#include "maxorb.h"
    ! need NLRGSH NBCH NUCO NRCO
#include "shells.h"
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
#include "implicit.h"
    integer,          intent(in)  :: ncgto, nectr
    type(cgto),       intent(out) :: bas(ncgto)
    real(8), target, intent(out) :: ectr(nectr)
    ! need MXSHEL
#include "maxorb.h"
    ! need NLRGSH NBCH NUCO NRCO NUMCF CENT NHKT NSTRT
#include "shells.h"
    ! need MXCONT
#include "aovec.h"
    ! need PRIEXP PRICCF
#include "primit.h"
    ! need MAXQNM
#include "maxmom.h"
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
#include "implicit.h"
    integer,  intent(in) :: ic
    real(8), intent(in) :: dc !step
    ! need MXSHEL
#include "maxorb.h"
    ! need CENT
#include "shells.h"
    ! need MXCOOR
#include "mxcent.h"
    ! need CORD
#include "nuclei.h"
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
