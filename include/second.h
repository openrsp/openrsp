C
C     type declaration for the function SECOND()
C
C     note that SECOND() is an internal function for the g77 compiler
C     and should then NOT be declared explicitly...
C
#if defined (SYS_CRAY)
      REAL SECOND
#else
#if !defined (VAR_G77)
      DOUBLE PRECISION SECOND
#endif
#endif
