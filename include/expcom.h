#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      REAL
#else
      DOUBLE PRECISION
#endif
     &        SIGN1X, SIGN1Y, SIGN1Z, SIGN2X, SIGN2Y, SIGN2Z,
     &        SIGN3X, SIGN3Y, SIGN3Z, SIGN4X, SIGN4Y, SIGN4Z
      LOGICAL TWOCEN, THRCEN, FOUCEN, DERONE, DERTWO
      COMMON /EXPCOM/ SIGN1X, SIGN1Y, SIGN1Z, SIGN2X, SIGN2Y, SIGN2Z,
     &                SIGN3X, SIGN3Y, SIGN3Z, SIGN4X, SIGN4Y, SIGN4Z,
     &                NCENT1, NCENT2, NCENT3, NCENT4,
     &                ICENT1, ICENT2, ICENT3, ICENT4,
     &                ISO1,   ISO2,   ISO3,   ISO4,
     &                DERONE, DERTWO, TWOCEN, THRCEN, FOUCEN
