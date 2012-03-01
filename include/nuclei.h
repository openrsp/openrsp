!
!     FILE: nuclei.h
!
      REAL*8  CHARGE, CORD, GNUEXP
      INTEGER NUCPRE, NUCNUM, NUCDEG, ISTBNU, NCTOT,                    &
     &        NUCIND, NUCDEP, NTRACO, ITRACO, NATOMS, NFLOAT,           &
     &        NBASIS, NLARGE, NSMALL, NPBAS,  NPLRG,  NPSML,            &
     &        NCHTOT, INCENT, INUNIQ, NDEGNM, ISOTOP, IZATOM,           &
     &        NBASISAUX, NPBASAUX,    NAUX,   NPAUX
      INTEGER MULBSI
!     MULBSI has been added for multiple basis sets (WK/UniKA/31-10-2002).
      LOGICAL NOORBT, GAUNUC
      COMMON /NUCLEI/ CHARGE(MXCENT), CORD(3,MXCENT), GNUEXP(MXCENT),   &
     &                NUCPRE(MXCENT), NUCNUM(MXCENT,8), NUCDEG(MXCENT), &
     &                ISTBNU(MXCENT), NDEGNM(MXCENT), NUCIND, NUCDEP,   &
     &                NTRACO, ITRACO(3),                                &
     &                NATOMS, NFLOAT, NBASIS, NLARGE, NSMALL, NPBAS,    &
     &                NPLRG, NPSML, NCHTOT, INCENT(MXCENT), NCTOT,      &
     &                INUNIQ(MXCENT), ISOTOP(MXCENT),IZATOM(MXCENT),    &
     &                NBASISAUX, NPBASAUX, NAUX, NPAUX,                 &
     &                MULBSI(MXCENT),                                   &
     &                GAUNUC, NOORBT(MXCENT)
      CHARACTER       NAMN*4, NAMEX*6, NAMDEP*6, NAMDPX*8
      COMMON /NUCLEC/ NAMN(MXCENT),   NAMEX(MXCOOR),                    &
     &                NAMDEP(MXCENT), NAMDPX(MXCOOR)
!  -- end of nuclei.h --
