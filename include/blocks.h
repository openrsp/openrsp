      INTEGER  MAXSHL, NLRGBL, NSMLBL, NHKTSH, KHKTSH, KCKTSH,
     &         ISTBSH, NUCOSH, NORBSH, NSTRSH, NCNTSH, NSETSH,
     &         IORBSB, NRCSH,  LCLASH, NO2INT, NLRBL,  ISYMBL, NSYMBL,
     &         MBSISH
#if defined (SYS_CRAY)
      REAL CENTSH
#else
      DOUBLE PRECISION CENTSH
#endif
      LOGICAL BIGVEC, SEGMEN, SEGMSH, SPHRSH
      COMMON /BLOCKS/ CENTSH(MXSHEL,3),
     &                MAXSHL, BIGVEC, SEGMEN,NLRGBL,NSMLBL,
     &                NHKTSH(MXSHEL), KHKTSH(MXSHEL), KCKTSH(MXSHEL),
     &                ISTBSH(MXSHEL), NUCOSH(MXSHEL), NORBSH(MXSHEL),
     &                NSTRSH(MXSHEL), NCNTSH(MXSHEL), NSETSH(MXSHEL,2),
     &                IORBSB(0:MXCORB-1), NRCSH(MXSHEL), SEGMSH(MXSHEL),
     &                LCLASH(MXSHEL), SPHRSH(MXCORB),
     &                NLRBL,ISYMBL(MXSHEL,8),NSYMBL,
     &                MBSISH(MXSHEL)
C     MBSISH has been added for multiple basis sets (WK/UniKA/31-10-2002).
