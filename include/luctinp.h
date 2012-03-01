! --- luctinp.h ---
!
!     Common block for LUCITA input information,
!     set in lucita/lucidal_input.F:LUCITAINP
!     Information is transferred to LUCITA routine dirluct.F
!     for further processing.
!
!     CRDINA : Card specifying inactive orbitals
!     CRDGAS : INGASD cards specifying GAS orbitals
!     CRDGOC : INGASD cards specifying min and max GAS occupation
!     CRDFRO : Card specifying frozen orbitals
!     CRDRS1/2/3: Cards specifying RAS1/2/3 orbitals
!
      INTEGER NTABLE
      PARAMETER (NTABLE = 30)

      INTEGER MXNGAS
      PARAMETER (MXNGAS = 16)

      CHARACTER*72 WAFFCD, CALCTP, SZCALD, TITLUC, CRDINA,              &
     &             CRDGAS(MXNGAS), CRDGOC(MXNGAS), CRDFRO,              &
     &             CRDRS1, CRDRS2, CRDRS3

      COMMON /LUCTINFC/ CRDINA, CRDGAS, CRDGOC, CRDFRO, CRDRS1, CRDRS2, &
     &                  CRDRS3, WAFFCD, CALCTP, SZCALD, TITLUC

      INTEGER IMOKW, NROOTD, ISSYMD, NACTED, IMULTD, IPRNGD, IPRNLD,    &
     &        IDENSD, MXHL1D, MXEL3D, INGASD, NSEQCD, IRSTLT, MXCIVE,   &
     &        IPARMODEL, ICIMAXITER, IMAXBLKSIZE,                       &
     &        I_USE_DIST_ROUTE, IN_MEMFAC
      COMMON /LUCTINFI/IMOKW(NTABLE), NROOTD, ISSYMD, NACTED, IMULTD,   &
     &                 IPRNGD,IPRNLD,IDENSD,MXHL1D,MXEL3D,INGASD,NSEQCD,&
     &                 IRSTLT, MXCIVE, IPARMODEL, ICIMAXITER,           &
     &                 IMAXBLKSIZE, I_USE_DIST_ROUTE, IN_MEMFAC
      REAL*8 CTRUNC_FAC
	COMMON /LUCTINFR/ CTRUNC_FAC
! --- end of luctinp.h ---
