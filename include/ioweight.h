#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
      DATA DAWEIG /1./, SQWEIG /1./
#else
      DATA DAWEIG /1.D0/, SQWEIG /1.D0/
#endif
