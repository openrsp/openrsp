C
C     For adding two-electron Darwin terms to electron repulsion integrals
C     (or one-electron four-center overlap integrals)
C
      INTEGER LTUV34
#if defined (SYS_CRAY)
      REAL DPTFAC,DARFAC 
#else
      DOUBLE PRECISION DPTFAC,DARFAC
#endif
      LOGICAL DO2DAR, AD2DAR, S4CENT, S4CDUM, DERNAI, FINDPT, NO2DPT,
     &        DPTINT, BPH2OO, NOPICH
      COMMON /DARTWO/ DARFAC, DO2DAR, AD2DAR, S4CENT, S4CDUM, DERNAI, 
     &                FINDPT, DPTFAC, NO2DPT, LTUV34, DPTINT,
     &                BPH2OO, NOPICH
