! File : inftap.h
!
!     Most file variables used in Dalton
!     ("TAP" for tapes here for historical reasons)
!
      INTEGER         LUCME,  LUMOL,  LUPROP, LUSOL,  LUINTA,           &
     &                LUONEL, LUSUPM, LUTLM,  LUDA1,  LUITMP,           &
     &                LU2DER, LUDASP, LURDR,  LURDI,  LUGDR,  LUGDI,    &
     &                LUGDT,  LURDT,  LUDFCK, LUSFDA, LUFDC,  LUWLK,    &
     &                LUPAO,  LUPAS,  LUNR1,  LUNR3,  LUNR5,            &
     &                LUINTR, LUMOLDEN,LUR12, LUORDA, LUMINT, LUSRINT
      INTEGER         LUAHSO, LUCRV1, LUCRV2, LUXYVE, LUCRVE,           &
     &                LURSP3, LURSP4, LURSP5, LUMHSO, LURSP
      INTEGER         LUINTM, LUIT1,  LUIT2,  LUIT3,  LUIT5,  LUINF,    &
     &                LUH2AC, LUSIFC, LBINTM, LBINTD, LBONEL, LBINFO
      INTEGER         LUPMOM, LUMOM,  LUEIND, LUENUC, LUESITE,LUEOBAR,  &
     &                LUVDWSE,LUENSA, LUQM3E, LUQM3P, LUOSCR, LUMMOV,   &
     &                LUOVER, LUNDPF, LUNMPF
      CHARACTER ABATLM*10, FNSOL*8,   ABARDR*9,  ABARDI*10, ABAGDR*9,   &
     &          ABAGDI*10, ABAGDT*10, ABARDT*10, ABADFK*10, ABASF*9,    &
     &          ABAWLK*10, ABATRJ*10, ABAIRC*10, ABANR1*10, ABANR3*10,  &
     &          ABANR5*10, FNONEL*8,  FNINTM*8,  FNSUPM*8,  FNSIFC*6,   &
     &          LBSIFC*8,  CLIND*8,   QMIND*8,   EINDFILE*8,ENUCFILE*8, &
     &          ESITFILE*8,OBARFILE*8,VDWSIEP*8, ENSAFILE*8,            &
     &          ELFDMM*8,  POTMM*8,   QM3CORD*8, MMOVER*6,MMOVER_INFO*11
      COMMON /INFTAP/ LUCME,  LUMOL,  LUPROP, LUSOL, LUINTA,            &
     &                LUONEL, LUSUPM, LUTLM,  LUDA1, LUITMP,            &
     &                LU2DER, LUDASP, LURDR,  LURDI, LUGDR, LUGDI,      &
     &                LUGDT,  LURDT,  LUDFCK, LUSFDA,LUFDC, LUWLK,      &
     &                LUPAO,  LUPAS,  LUNR1,  LUNR3, LUNR5,             &
     &                LUINTR, LUMOLDEN, LUR12(20),                      &
     &                LUPMOM, LUMOM,  LUEIND, LUENUC, LUESITE, LUEOBAR, &
     &                LUVDWSE,LUENSA, LUQM3E, LUQM3P, LUOSCR,LUMMOV,    &
     &                LUOVER, LUNDPF, LUSRINT,LUNMPF
      COMMON /RSPTAP/ LUAHSO, LUCRV1, LUCRV2, LUXYVE, LUCRVE,           &
     &                LURSP3, LURSP4, LURSP5, LUMHSO, LURSP
      COMMON /SIRTAP/ LUINTM, LUIT1,  LUIT2,  LUIT3,  LUIT5,  LUINF,    &
     &                LUH2AC, LUSIFC, LUORDA, LUMINT,                   &
     &                LBINTM, LBINTD, LBONEL, LBINFO, LBSIFC
      COMMON /CHRTAP/ ABATLM, FNSOL,  ABARDR, ABARDI, ABAGDR, ABAGDI,   &
     &                ABAGDT, ABARDT, ABADFK, ABASF,  ABAWLK, ABATRJ,   &
     &                ABAIRC, ABANR1, ABANR3, ABANR5, FNINTM, FNONEL,   &
     &                FNSUPM, FNSIFC, CLIND,  QMIND,  EINDFILE,         &
     &                ENUCFILE, ESITFILE, OBARFILE, VDWSIEP,            &
     &                ENSAFILE,ELFDMM,POTMM,QM3CORD,MMOVER,MMOVER_INFO
! -- end of inftap.h --
