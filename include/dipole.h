#if defined (SYS_CRAY)
      REAL DIPMN, DIPME, DDIPN, DDIPE, DDIPS, DDIPR
      REAL DIPME2, DIPME3, DLFN, DLFE
#else
      DOUBLE PRECISION DIPMN, DIPME, DDIPN, DDIPE, DDIPS, DDIPR
      DOUBLE PRECISION DIPME2, DIPME3, DLFN, DLFE
#endif
      COMMON /DIPOLE/ DIPMN(3), DIPME(3), DLFN(3), DLFE(3),              &
     &                DDIPN(3,MXCOOR), DDIPE(3,MXCOOR), DDIPS(3,MXCOOR), &
     &                DDIPR(3,MXCOOR),                                   &
     &                DIPME2(3), DIPME3(3)
