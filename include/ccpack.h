C 
C------------------------------------------------------------
C information for the packing of the AO integrals on the
C presorted CCAOIN_* files
C------------------------------------------------------------
C
      LOGICAL LPACKINT 
      INTEGER IPCKTABINT(0:255)
      INTEGER IOFFINT(MXCORB)
      INTEGER NPCKINT(MXCORB)
      INTEGER NTOTINT, NTOTPCK

#if defined (SYS_CRAY)
      REAL THRPCKINT, PCKRATIO, PCKTIME
#else
      DOUBLE PRECISION THRPCKINT, PCKRATIO, PCKTIME
#endif
      COMMON /CCPACK/ THRPCKINT, PCKRATIO, PCKTIME,
     &                IPCKTABINT, IOFFINT, NPCKINT, 
     &                NTOTINT, NTOTPCK, LPACKINT

