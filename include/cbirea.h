!
!     File: cbirea.h - common block for reading of MOLECULE.INP (herrdn.F)
!     Last revision: Apr 2010 /hjaaj
!
!     MAXPRD = default value for MAXPRI
      INTEGER MAXPRD
      PARAMETER ( MAXPRD = 35 )
!     MAXFAMEXP = maximum number of exponents in family basis sets
      INTEGER MXFAMEXP
      PARAMETER (MXFAMEXP = 100)

      REAL*8  ZCMVAL, TOL_SYMADD, FAMEXP(MXFAMEXP, 2), FAMPAR(4)
      INTEGER IPREAD, MAXPRI, LUMLCL, LCMMAX, NCMSTR, NCMEND, NFAMEXP(2)
      LOGICAL BIGVC,  GENCON, INPTST, DIRAC,  BASIS,  PRIBAS, ATOMBA,   &
     &        UNCONT, WRTLIN, ANGS,   ATOMDF, CNTMAT, NODDYDF
      COMMON /CBIREA/ ZCMVAL, TOL_SYMADD,     FAMEXP, FAMPAR,           &
     &        IPREAD, MAXPRI, LUMLCL, LCMMAX, NCMSTR, NCMEND, NFAMEXP,  &
     &        BIGVC,  GENCON, INPTST, DIRAC,  BASIS,  PRIBAS, ATOMBA,   &
     &        UNCONT, WRTLIN, ANGS,   ATOMDF, CNTMAT, NODDYDF
!
!     Info for multiple basis sets (p.t. used with r12int.h in Dalton)

      INTEGER    MXMULB
      PARAMETER (MXMULB = 2)
      CHARACTER*80 MULNAM(MXMULB)
      COMMON /CBIREA_C/ MULNAM

      INTEGER NMULBS
      LOGICAL LMULBS
      COMMON /CMMBAS/ NMULBS, LMULBS
! -- end of cbirea.h --
