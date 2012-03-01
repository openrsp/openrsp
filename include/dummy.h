! --- dummy.h ---
      REAL*8  DUMMY
      REAL*8  DUMMY2(2)
      INTEGER IDUMMY
!     make DUMMY write protected (on most computers)
      PARAMETER ( DUMMY = 1.0D20 , IDUMMY = - 9999999 )
! --- end of dummy.h ---
