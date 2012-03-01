      LOGICAL SKIP,   DIRECT, DIRGRD, CCRSTR, 
     &        FROIMP, FROEXP, NOCCIT,
     &        CCSAVE, STOLD,  JACEXP, LHTR,
     &        DEBUG,  CCSTST, ANAAOD,
     &        HERDIR, FREEZE, KEEPAOIN, NOEONL, NOSORT,
CSPAS 15.11.2009 implementing AO-SOPPA
C    &        SIRSOP, LVVVV, ONLYMO
     &        SIRSOP, LVVVV, ONLYMO, AOSOPPA
CKeinSPASmehr

      LOGICAL ETACCPT, DIRKAPB

      INTEGER MXDIIS, MXLRV,
     &        ICHANG, IPRINT, KEEPAOTWO

      LOGICAL CCS,  CIS,  MP2,  CC2,   CC1A,  CC1B,
     &        MCC2, CCP2, CCD,  CCSD,  CC3,   CCPT,
     &        CCP3, CCRT, CCR3, CCR1A, CCR1B, CCT,
     &        CCSDT,CCR12, MTRIP, CHOPT

      COMMON /CCSDGNINP/ SKIP,   DIRECT, DIRGRD, CCRSTR,
     &                   FROIMP, FROEXP, NOCCIT,
     &                   CCSAVE, STOLD,  JACEXP, LHTR,
     &                   DEBUG,  CCSTST, ANAAOD,
     &                   MXDIIS, MXLRV,
     &                   ICHANG, IPRINT, KEEPAOTWO, HERDIR,
     &                   ETACCPT, DIRKAPB,
     &                   FREEZE, KEEPAOIN, NOEONL, NOSORT,
CSPAS 15.11.2009 implementing AO-SOPPA
C    &                   SIRSOP, LVVVV, ONLYMO
     &                   SIRSOP, LVVVV, ONLYMO, AOSOPPA
CKeinSPASmehr

      COMMON /CCMODELS/ CCS,  CIS,  MP2,  CC2,   CC1A,  CC1B,
     &                  MCC2, CCP2, CCD,  CCSD,  CC3,   CCPT,
     &                  CCP3, CCRT, CCR3, CCR1A, CCR1B, CCT,
     &                  CCSDT,CCR12,MTRIP, CHOPT

