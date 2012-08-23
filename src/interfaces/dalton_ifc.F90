! radovan: do not add anything here
!          this module is about to retire

module dalton_ifc

  use interface_f77_memory
  use interface_io

  implicit none

  public vibctl_ifc

  private

  contains

#ifdef VAR_LSDALTON
  !> \brief unknown function
  !> \author Thomas Kjaergaard
  !> \date 2012-08-23
  subroutine VIBCTL_ifc( nc, w, ALPHA, GPRIME, THETA, DIPGRAD, &
                         dALPHAdR, dGPRIMEdR, dTHETAdR )
    implicit none
    integer :: nc
    real(8) :: w
    real(8) :: ALPHA(3,3), GPRIME(3,3), THETA(3,6)
    real(8) :: DIPGRAD(3,nc), dALPHAdR(3,3,nc)
    real(8) :: dGPRIMEdR(3,3,nc), dTHETAdR(3,6,nc)
    stop 'VIBCTL_ifc not implemented in LSDALTON interface. TK'
  end subroutine VIBCTL_ifc
#else
  !> \brief gets the
  !> \author Bin Gao
  !> \date 2009-12-12
  !> \param
  !> \return
  !> \note modified on linsca/linears/VIBCTL_interface.F
  subroutine VIBCTL_ifc( nc, w, ALPHA, GPRIME, THETA, DIPGRAD, &
                         dALPHAdR, dGPRIMEdR, dTHETAdR )
#include "mxcent.h"
#include "cbilnr.h"
#include "cbivib.h"
#include "abainf.h"
#include "moldip.h"
    integer :: nc
    real(8) :: w
    real(8) :: ALPHA(3,3), GPRIME(3,3), THETA(3,6)
    real(8) :: DIPGRAD(3,nc), dALPHAdR(3,3,nc)
    real(8) :: dGPRIMEdR(3,3,nc), dTHETAdR(3,6,nc)
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
#endif

end module
