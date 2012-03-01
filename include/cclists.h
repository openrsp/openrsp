*======================================================================*
*    
*     Purpose: centralize info about the leading dimensions of the
*              transformation lists
*
*
*======================================================================*

      INTEGER MXDIM_FTRAN  ! F matrix transformation, CC_FMATRIX
      PARAMETER (MXDIM_FTRAN  = 3)

      INTEGER MXDIM_JTRAN  ! J matrix transformation, CC_JMATRIX
      PARAMETER (MXDIM_JTRAN  = 3)

      INTEGER MXDIM_GTRAN  ! G matrix tranformation, CC_GMAT
      PARAMETER (MXDIM_GTRAN  = 4)

      INTEGER MXDIM_HTRAN  ! H matrix tranformation, CC_HMAT
      PARAMETER (MXDIM_HTRAN  = 5)


      INTEGER MXDIM_BTRAN  ! B matrix tranformation, CC_BMATRIX
      PARAMETER (MXDIM_BTRAN  = 3)

      INTEGER MXDIM_KTRAN  ! K matrix tranformation, CC_KMATRIX
      PARAMETER (MXDIM_KTRAN  = 3)

      INTEGER MXDIM_CTRAN  ! C matrix tranformation, CC_CMAT
      PARAMETER (MXDIM_CTRAN  = 4)

      INTEGER MXDIM_DTRAN  ! D matrix tranformation, CC_DMAT
      PARAMETER (MXDIM_DTRAN  = 5)


      INTEGER MXDIM_XEVEC  ! Xksi{O} and Eta{O} vectors, CC_XIETA
      PARAMETER (MXDIM_XEVEC  = 8)

      INTEGER MXDIM_FATRAN ! F{O} matrix tranformation, CCQR_FADRV
      PARAMETER (MXDIM_FATRAN = 5)

      INTEGER MXDIM_AATRAN ! A{O} matrix transformation, CC_AAMAT
      PARAMETER (MXDIM_AATRAN = 7)

      INTEGER MXDIM_BATRAN ! B{O} matrix tranformation, CC_BAMAT
      PARAMETER (MXDIM_BATRAN = 8)

      INTEGER MXDIM_CATRAN ! C{O} matrix tranformation, CC_CAMAT
      PARAMETER (MXDIM_CATRAN = 9)

*======================================================================*
