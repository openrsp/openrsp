!
! FILE: cbilnr.h
!
! Used to control linear response solver in ABACUS
!
      LOGICAL         ALFA, ROAA, ROAG, STATIC
      PARAMETER       (MXFR  = 101, MAXLN = 80)
      CHARACTER*8     LABALN
      COMMON /LNLBL / LABALN(MAXLN)
      COMMON /CBILNR/ THCLNR, FRVAL(MXFR), NFRVAL, IPRLNR, NABALN,       &
     &                ALFA, ROAA, ROAG, STATIC
! -- end of cbilnr.h --
