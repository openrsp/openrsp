!     FILE: orgcom.h
!
! Common block with origins:
!
!     TOTMAS, TNUCMAS : total mass, total nuclear mass
!     CMXYZ : center of mass
!     ORIGIN: center of origin of most operators
!     DIPORG: center of origin of dippole operators
!     GAGORG: center of origin of gauge for magnetic operators
!     CAVORG: center of origin of cavity in continuum models
!     FMMORI: use local origins for moment integrals (used for FMM)
!
      REAL*8  TOTMAS, TNUCMAS
      REAL*8  CMXYZ(3), ORIGIN(3), DIPORG(3), GAGORG(3), CAVORG(3)
      LOGICAL FMMORI
      COMMON /ORGCOM/ TOTMAS, TNUCMAS, CMXYZ, ORIGIN, DIPORG, GAGORG,   &
     &                CAVORG,          FMMORI
! -- end of orgcom.h --
