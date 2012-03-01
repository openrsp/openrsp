! File: ibtfun.h
!
! Purpose: redefine bit functions as statement functions,
!          because some systems had chosen different names
!
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = SHIFTL(I,J)
      IBTSHR(I,J) = SHIFTR(I,J)
      IBTXOR(I,J) = XOR(I,J)
#else
      IBTAND(I,J) = IAND(I,J)
      IBTOR(I,J)  = IOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = IEOR(I,J)
#endif
! -- end of ibtfun.h --
