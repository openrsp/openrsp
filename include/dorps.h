C
C     File: dorps.h
C     Purpose: control which nuclear coordinates to treat in ABACUS.
C
C     Parameter NDORL must be updated after changes (for parallelization)
C     NDORL is length of DOREPS and DOCOOR together
C     NOTE: There must be no variables between DOREPS and DOCOOR.
C
      INTEGER NDORL, NDCORD
      PARAMETER (NDORL = 8 + 3*MXCENT)
      LOGICAL DOREPS, DOCOOR, DCORD, DCORGD, DOPERT
      COMMON /DORPS/ NDCORD(2),
     &               DOREPS(0:7), DOCOOR(3,MXCENT),
     &               DCORD(MXCENT,3,2), DCORGD(0:MXCENT,3,2),
     &               DOPERT(0:3*MXCENT,2)
