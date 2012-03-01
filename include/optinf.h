!
! File: optinf.h
!
! Information for geometry optimization
! controlled in abaopt.F, abaop2.F, and abarint.F
!
      PARAMETER (MAXPRE = 10)
      CHARACTER*80    SPBSTX, PREBTX
      LOGICAL GECONV, NOTRST, NOBRKS, BRKSYM, NWSYMM, DOSPE,  DOPRE,
     &                FINPRE, VRML,   VRBOND, VREIGV, VRCORD, VRVIBA,
     &                VRML_SYM, VISUAL, BOFILL, INITHS, HSFILE, BFGSR1,
     &                STEEPD, RANKON, PSB,    DFP,    BFGS,   SCHLEG, 
     &                NEWTON, QUADSD, KEEPHE, BAKER,  REDINT, CARTCO,
     &                INRDHS, FSTORD, SNDORD, REJINI, GRDINI, MULTI,
     &                CHGRDT, CONOPT, MODHES, INMDHS, FINDRE, TRSTRG,
     &                RATFUN, GDIIS,  DELINT, RSTARR, LNSRCH, SADDLE,
     &                REBILD, CMBMOD, HFPROP, CONFRM, NOAUX,  NODIHE,
     &                LINDHD, ADDCRD, REDRED, NOADDA, NATNRM, NOHSWR
      COMMON /OPTINF/ TRSTRA, TRSTIN, TRSTDE, RTENBD, RTENGD, RTRJMN,
     &                RTRJMX, ENERGY, ERGOLD, ERGPRD, ERGPRO, STPNRM,
     &                STPNRO, GRADNM, THRERG, GRDTHR, THRSTP, THRSYM,
     &                THGRMX, THSTMX, EVLINI, DISPLA, PRVRMS, PRVMAX,
     &                STPDIA(8*MXCENT), 
     &                STPSYM(8*MXCENT), STPINT(8*MXCENT), 
     &                GRDDIA(8*MXCENT), EVAL(8*MXCENT), 
     &                EVALOL(8*MXCENT), GRDINT(8*MXCENT),
     &                CRDIN1(8*MXCENT), CRDINT(8*MXCENT),
     &                CNDHES(0:7), INDHES(0:7), INTCRD(8*MXCENT,6),
     &                ICONF(0:5), ICNSTR(8*MXCENT), IADDCR(0:10,1:4),
     &                IFREEZ(0:10),
     &                ISTBLZ, IAUXRD, ITOTRJ, KEPTIT, NSPMOD, NCNSTP,
     &                INDTOT, ITRNMR, ITRMAX, MAXREJ, IPRINT, NCRTOT,
     &                NCART,  NPROJ,  NTMAT,  IINTCR, IREDIC, ICRTCR,
     &                ICONDI, ITRBRK, NUMPRE, IPRE,   ITRFRZ,
     &                PREBTX(MAXPRE),
     &                SPBSTX, GECONV, NOTRST, NOBRKS, BRKSYM, NWSYMM,
     &                DOSPE,  DOPRE,  FINPRE, VRML,   VRBOND, VREIGV,
     &                VRCORD, VRVIBA, VRML_SYM, VISUAL, INITHS, HSFILE,
     &                BFGSR1, STEEPD, RANKON, PSB   , DFP,    BFGS,
     &                SCHLEG, NEWTON, QUADSD, KEEPHE, BAKER,  REDINT,
     &                CARTCO, INRDHS, FSTORD, SNDORD, REJINI, GRDINI,
     &                MULTI,  CHGRDT, CONOPT, MODHES, INMDHS, FINDRE,
     &                TRSTRG, RATFUN, GDIIS,  DELINT, RSTARR, LNSRCH,
     &                SADDLE, REBILD, BOFILL, CMBMOD, HFPROP, CONFRM,
     &                NOAUX,  NODIHE, LINDHD, ADDCRD, REDRED, NOADDA,
     &                NATNRM, NOHSWR
! --- end of optinf.h ---
