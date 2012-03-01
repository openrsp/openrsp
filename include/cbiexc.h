!
!     cbiexc.h - Control common block for abacus/abaexc.F
!
      LOGICAL         SKIP, CUT, DIPSTR, ROTSTR, ROTVEL, FNAC, EXCTRP,
     &                SUMRUL, OOTV
      PARAMETER       (MAXPP = 200)
      CHARACTER*8     LABAPP
      COMMON /PPLBL / LABAPP(MAXPP), LABSYM(MAXPP)
      COMMON /CBIEXC/ THREXC,
     &                NEXCIT(8), MAXITE, MXNEXI, MXRM,
     &                MXPHP, NABAPP, IPREXC, IPR1IN,
     &                SKIP, CUT, DIPSTR, ROTSTR, ROTVEL, FNAC, EXCTRP,
     &                SUMRUL, OOTV
! -- end of abaexc.h --
