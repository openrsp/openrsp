!  FILE : ccdeco.h
      INTEGER MXCHVC, LRDTOT, LREDU, NREDUC, ISYSCR, INAOSH, NUMCHO,
     &        LENCHO, IDNTCH, NDIAG, IDIAG, MAXDIA, MAXDI1, IATPRI,
     &        NUMCEN, IALBET, CHOVER, IPRCHO
      REAL*8  DIASCR, THRCOM, THRDEF, THINDI, THSUDI, SPAN

      PARAMETER (MXCHVC = 80 000)
C
      LOGICAL CHOINT,COMP,RSTDIA,RSTCHO,SCDIAG,REDUCE,DIACAL,NEWSCF
C
      COMMON /CHOINT/ DIASCR(MXSHEL,MXSHEL),
     &                THRCOM, THRDEF, THINDI, THSUDI, SPAN,
     &                CHOVER, IPRCHO,
     &                LRDTOT,LREDU,NREDUC(8),
     &                ISYSCR(MXSHEL,MXSHEL),INAOSH(MXCORB),
     &                NUMCHO(8), LENCHO(MXCHVC,8),IDNTCH(MXCHVC,8),
     &                NDIAG,IDIAG(8),MAXDIA,MAXDI1,
     &                IATPRI(MXCORB),NUMCEN,
     &                IALBET(MXCORB*(MXCORB+1)/2,2), 
     &                CHOINT,COMP,RSTDIA,RSTCHO,SCDIAG,REDUCE,DIACAL,
     &                NEWSCF
C
C     Old common with some stuff to have several files.
C
C     COMMON /CHOINT/ DIASCR(MXSHEL,MXSHEL),
C    *                THRCOM,THRDEF,THINDI,THSUDI,SPAN,
C    *                ISYSCR(MXSHEL,MXSHEL),
C    *                NUMCHO(8),IADRTO(8),IDNTLU(0:9,1:8),
C    *                INDCHO(MXCHVC,8),LENCHO(MXCHVC,8),
C    *                IDNTCH(MXCHVC,8),ISTRFI(0:9,1:8),
C    *                IDIAG(8),IIBST(8),
C    *                NDIAG,INAOSH(MXCORB),
C    *                MAXDIA,MAXDI1,NUMFIL(8),
C    *                IATPRI(MXCORB),NUMCEN,
C    *                IALBET(MXCORB*(MXCORB+1)/2,2), 
C    *                CHOINT,COMP,RSTDIA,RSTCHO,SCDIAG
C
!  -- end of ccdeco.h --
