C File : cbieri.h
C
C     Parameters NCBI? must be updated after changes (for parallelization)
C
C     NOTE: New logicals should appear after the last logical (NCLERI)
C           New integers should appear after the last integer (IFITDM)
C           Reals should appear at the end.
C
      INTEGER NCBII,  NCBIL
      PARAMETER (NCBII = 18,
     &           NCBIL = 30)
      INTEGER IPRERI, IAOBCH, IPRNT1, IPRNT2, IPROD1, IPROD2
      INTEGER LBFINP, MAXDST, NDMAT,  IANGMO, NSPMAX, MAXDSD,
     &        MXBCH,  MXBCH0, IFITDM
      LOGICAL RUNERI, TIMERI, PMSAB,  PMSCD,  PMS12,  RTNERI, OFFCNT,
     &        DIASRT, OLDCR1, CANIND, WRTSCR, NOWRIT, EXTPRI,
     &        INTPRI, DODIST, INTSKP, DISTST, WRTINT, FCKINT, PROFIL,
     &        NOLOCS, GRDZER, OLDDER, EXPERI, DOERIP, ERITWO, CCRUN,
     &        COMPRS, GENCON_ERI, NCLERI, DGABAB
      DIMENSION IANGMO(4)
      COMMON /CBIERI/ RUNERI, TIMERI, PMSAB,  PMSCD,  PMS12,  RTNERI,
     &                OFFCNT, DIASRT, OLDCR1, CANIND, WRTSCR, NOWRIT,
     &                EXTPRI, INTPRI, DODIST, INTSKP, DISTST, WRTINT,
     &                FCKINT, PROFIL, NOLOCS, GRDZER, OLDDER, EXPERI,
     &                DOERIP, ERITWO, CCRUN,  COMPRS, GENCON_ERI,
     &                NCLERI, DGABAB,
     &                NDMAT,  IPROD1, IPROD2, IAOBCH, IPRNT1, IPRNT2,
     &                IPRERI, LBFINP, MAXDST, IANGMO, NSPMAX, MAXDSD,
     &                MXBCH,  MXBCH0, IFITDM
C -- end of cbieri.h --
