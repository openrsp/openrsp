C
C     Parameters NINT? must be updated after changes (for parallelization)
C
      PARAMETER (NINTI = 14,
     &           NINTL = 7)
      LOGICAL JNODV, JNOPV, JNOCNT, JRETUR, JTKTIM, JSOLVN, JRELCL
      COMMON /PARINT/ JATOM,  JLUDAS, JLUINT, JLUONE, JLUSOL, JLUSUP,
     &                JMXDIF, JMXREP, JNFMAT, JTASK,  JTYPE,  J2TYP,
     &                JCEDIF, JFTHRS,
     &                JNODV,  JNOPV,  JNOCNT, JRETUR, JTKTIM, JSOLVN,
     &                JRELCL
