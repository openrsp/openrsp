      PARAMETER (NPSORT = 9)
      LOGICAL DOSORT
      COMMON /ODBTCH/ NODBCH, NODCLS, NODPPR, DOSORT(NPSORT)
C     NPSORT has been increased from 8 to 9 for basis-set identifiers 
C     (WK/UniKA/31-10-2002).
