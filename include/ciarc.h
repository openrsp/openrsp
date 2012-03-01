C ==========================================
C Configuration of CC_CIA
C See subroutine CC_CIA for more information
C ==========================================
C IPRLVL: Print level threshold.
      PARAMETER (IPRLVL = 3, NTMCIA = 6)
      LOGICAL SYMTRZ, LIDEN, CIAMIO, GETMNM, INXINT
      CHARACTER*9 CIAROU
C The following 6 common blocks should be initialized
C by the calling routine.
      COMMON / CIARCA / IOFA1(8),LVIRA(8),NX1AMA(8),IX1AMA(8,8)
      COMMON / CIARCB / IOFB1(8),LVIRB(8),NX1AMB(8),IX1AMB(8,8)
      COMMON / CIARCI / IX2SQ(8,8),NX2SQ
      COMMON / CIARCC / NTOVEC(8),ISYCH1,ISYCH2,ITYP1,ITYP2
      COMMON / CIARCP / IPRCIA
      COMMON / CIARCL / LIDEN, SYMTRZ,CIAMIO,GETMNM,INXINT
C The following 3 common blocks are internal to CC_CIA 
      COMMON / CIATIM / TIMCIA(NTMCIA)
      COMMON / CIABAT / NCHBAT(8)
      COMMON / CIAIDI / CIAROU
C ===========================
C End configuration of CC_CIA
C ===========================

