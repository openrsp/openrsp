      INTEGER CHANNEL_MAX
      PARAMETER (CHANNEL_MAX = 150)
      INTEGER CHANNEL_ORB(CHANNEL_MAX), CHANNEL_NORB,
     &     CHANNEL_VIRT(8)
C     CHANNEL_VIRT(ISYM) is the number of virtual orbitals to
C     include in a response calculation in each symmetry. 
C     Set CHANNEL_VIRT(ISYM) to -1 to include all orbitals in this 
C     symmetry
      LOGICAL CHANNEL_CALC,CHANNEL_VCALC
      COMMON /CHANNEL/ CHANNEL_ORB, CHANNEL_NORB, CHANNEL_VIRT,
     &     CHANNEL_CALC,CHANNEL_VCALC
