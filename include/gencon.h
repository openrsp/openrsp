!
! VERSION : $Revision: 10155 $
! DATE    : $Date: 2010-01-10 12:53:17 +0100 (Sun, 10 Jan 2010) $
! FILE    : gencon.h
!
!   - Parameter for maximum number of functions in a block:
!     MXCON_PARAM = NPRIM * l(l+1)/2 * NCONTR
!   - BSET_is_GENCON(i) true if contracted basis set is used for basis set "i"
!
      INTEGER MXCON_PARAM
      PARAMETER (MXCON_PARAM = 300)

      LOGICAL BSET_is_GENCON(-MXBSETS:MXBSETS)
      COMMON /GENCON/ BSET_is_GENCON
