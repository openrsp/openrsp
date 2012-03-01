C
C QRSHG  - SHG process                     : w_2= w_3
C QROPRF - optical refractivity calculation: w_2=-w_3
C
      PARAMETER ( MXQROP = MAXLBL , MBQRFR = MAXLBL , MCQRFR = MAXLBL )
      LOGICAL  HYPCAL,AQROP,BQROP,CQROP,REFCHK,QRSPEC,QRSHG,QRPOCK,
     *         QROPRF,SOCOLL,SSCOLL
      CHARACTER*8 AQRLB,BQRLB,CQRLB
      COMMON /INFHYP/ HYPCAL,IPRHYP, NBQRFR,
     *                NCQRFR, BQRFR(MBQRFR),
     *                CQRFR(MCQRFR),NAQROP(8), NBQROP(8), NCQROP(8) ,
     *                AQROP(MXQROP), BQROP(MXQROP), CQROP(MXQROP),
     *                REFCHK,IAABB,QRSPEC,QRSHG,QRPOCK,QROPRF,SOCOLL,
     *                SSCOLL
      COMMON /CHRHYP/ AQRLB(8,MXQROP), BQRLB(8,MXQROP),
     *                CQRLB(8,MXQROP)
