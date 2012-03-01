*---------------------------------------------------------------------*
*
*   Purpose: store input specific for excitation energies/eigenvector
*            (CC_EXCI section, read in CC_EXCINP)
*
*   NCCEXCI : # of requested states per symmetry & multiplicity
*   NCC     : # of ??
*
*   STVEC   : flag that start vectors are given by input
*   ISTVEC  : start vector information
*
*   OMEINP  : flag that omegas for triples calculations have been
*             given in the input
*   NOMINP  : # omegas given
*   IOMINP  : indeces?
*   EOMINP  : omegas...
*
*   FDEXCI  : flag for diagonlization of fin.-diff. jacobian
*   JACEXP  : flag for explicite construction of the jacobian
*   JACTST  : flag for test of the jacobian
*
*   STSD    : flag for .STSD option in triples calculation
*   OMESC   : flag for .NOSCOM option in triples calculation
*
* Christof Haettig, Januar 1999
*---------------------------------------------------------------------*
      INTEGER MAXOME
      PARAMETER (MAXOME = 50)

      LOGICAL FDEXCI, FDJAC, JACTST, STVEC, OMEINP, STSD, OMESC
      LOGICAL CCSPIC,CC2PIC,CCSDPI,MARGIN,SQROVLP,CCSDTRENRM
      LOGICAL EXCI_CONT

      INTEGER NCCEXCI(8,3), NOMINP(8,3)
      INTEGER IOMINP(MAXOME,8,3), ISTVEC(MAXOME,8)
C     INTEGER NCC(8,3)

#if defined (SYS_CRAY)
      REAL EOMINP(MAXOME,8,3), TOLSC
      REAL OMPCCS,OMPCC2,OMPCCSD,XMARGIN
#else
      DOUBLE PRECISION EOMINP(MAXOME,8,3), TOLSC
      DOUBLE PRECISION OMPCCS,OMPCC2,OMPCCSD,XMARGIN
#endif

      COMMON /CCEXCIT/ EOMINP, TOLSC,
     &                 NCCEXCI, NOMINP, IOMINP, ISTVEC,
     &                 FDEXCI, FDJAC, JACTST, STVEC, 
     &                 OMEINP, STSD, OMESC,
     &                 CCSPIC,CC2PIC,CCSDPI,MARGIN,SQROVLP,CCSDTRENRM,
     &                 EXCI_CONT
      COMMON /OMEPIC/ OMPCCS,OMPCC2,OMPCCSD,XMARGIN
*---------------------------------------------------------------------*
