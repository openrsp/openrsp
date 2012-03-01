C The PL1 list: 
C Common blocks for the projected L1=Zeta^X response 
C Lagrangian multipliers. List type: 'PL1'.....sonia
C Vector and operator have the same symmetry (ISYPL1)
C
      INTEGER MAXPL1LBL
      PARAMETER ( MAXPL1LBL = 60 )
      LOGICAL LORXPL1, LPRPL1, LPL1OPN
      INTEGER NPL1LBL, ISYPL1, ISYOFPL1
      INTEGER ISYSPL1, ISTPL1
      CHARACTER*8 LBLPL1

#if defined (SYS_CRAY)
      REAL FRQPL1, EIGPL1
#else
      DOUBLE PRECISION FRQPL1, EIGPL1
#endif
      !integer common block
      COMMON/IPL1RSP/ NPL1LBL,ISYPL1(MAXPL1LBL),ISYOFPL1(8),
     &                LPRPL1(MAXPL1LBL), LORXPL1(MAXPL1LBL), 
     &                LPL1OPN,
     &                ISYSPL1(MAXPL1LBL), ISTPL1(MAXPL1LBL)
      !character common block
      COMMON/CPL1RSP/ LBLPL1(MAXPL1LBL)
      !real common block
      COMMON/RPL1RSP/ FRQPL1(MAXPL1LBL), EIGPL1(MAXPL1LBL)
