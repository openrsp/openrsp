*----------------------------------------------------------------------*
* flags for test (noddy) implementation of CC3 response: 
*----------------------------------------------------------------------*
      CHARACTER*12 FILNODFOCK, FILNODINTS, FILNODINTT, FILNODINTV,
     &             FILNODINTX, FILNODT30,  FILNODL30
      PARAMETER (FILNODFOCK = 'CCNODDYFOCK0',
     &           FILNODINTS = 'CCNODDYINTS0',
     &           FILNODINTT = 'CCNODDYINTT0',
     &           FILNODINTV = 'CCNODDYINTV0',
     &           FILNODINTX = 'CCNODDYINTXY',
     &           FILNODT30  = 'CCNODDYT30AM',
     &           FILNODL30  = 'CCNODDYL30AM')

      LOGICAL NODDY_XI, NODDY_ETA, NODDY_XI_ALTER, NODDY_OVLP,
     &        NODDY_ETA_ALTER, NODDY_FMAT, CCSDT_F_ALTER, NODDY_OMEGA,
     &        NODDY_LHTR, NODDY_RHTR, NODDY_DEN, NODDY_FAMAT, 
     &        NODDY_GMAT, NODDY_BMAT, NODDY_AAMAT, NODDY_INIT,
     &        NODDY_HMAT, NODDY_FA_ALTER

      COMMON /CC_NODDY/ NODDY_XI, NODDY_ETA, NODDY_XI_ALTER, NODDY_OVLP,
     &        NODDY_ETA_ALTER, NODDY_FMAT, CCSDT_F_ALTER, NODDY_OMEGA,
     &        NODDY_LHTR, NODDY_RHTR, NODDY_DEN, NODDY_FAMAT,
     &        NODDY_GMAT, NODDY_BMAT, NODDY_AAMAT, NODDY_INIT,
     &        NODDY_HMAT, NODDY_FA_ALTER
*----------------------------------------------------------------------*
