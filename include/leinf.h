C     LEINF : space nedded in LELIN = LLEWA + nsim*LLEWB
C             NSIDE=0, symmetric; =1 from left side; =2 from right side
#if defined (SYS_CRAY)
      REAL THRLE
#else
      DOUBLE PRECISION THRLE
#endif
      INTEGER MAXLE, IPRLE, LETYPA, LETYPB, LETOT, NSIDE, NCREF,
     *        LULEA, LULEB, LULEC, LLEWA, LLEWB, KLELIN, NLELIN,LLELIN,
     *        LETYP2, LETYP1
C             NCREF = # of reference vectors to othogonalize against
      COMMON /LEINF / THRLE, MAXLE, IPRLE, LETYPA,LETYPB,LETOT, NSIDE,
     *                NCREF, LULEA, LULEB, LULEC,
     *                LLEWA, LLEWB, KLELIN(20),   NLELIN,LLELIN
      EQUIVALENCE (LETYP2,LETYPA), (LETYP1,LETYPB)
