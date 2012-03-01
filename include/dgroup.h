! File : dgroup.h
!
!     NZ - number of matrices (1,2,4)
!     NFSYM - number of fermion irreps
!     NBSYM - number of boson irreps
!     JSPINR(IM,IC,IFRP) - boson irrep for a given spinor position
!     JBTOF(IBRP,IC) - points to fermion irrep for given boson 
!                      irrep/component
!     JFSYM(NZ,IFRP) - points to boson symmetries associated with 
!                      given fermion
!                      ircop
!     IQDEF(4)
!     JQBAS(IBRP,IC) - points to quaternionic position
!     IPQTOQ(4,IBRP) - pointer from packed quaternion to quaternion unit
!     IQTOPQ(4,IBRP) - pointer from quaternion to packed quaternion unit
!     IHQMAT(4,IH, ) - matrix symmetry(symmetric/antisymmetric)
!     IRQMAT(4,IBRP)    - irrep of matrix
!
      CHARACTER*3 FREP(2)
      COMMON /DNAME/ FREP

      INTEGER NDGROUP
      PARAMETER (NDGROUP = 207)

      INTEGER I1_DGROUP, NZ, NFSYM, NBSYM, IFDIRP(4, 2), JBTOF(0:7, 2), &
     &        JMROI(4), JSPINR(4, 2, 2), JQBAS(0:7, 2), JFSYM(4, 2),    &
     &        IPQTOQ(4, 0:7), IQTOPQ(4, 0:7), IQPH(4, 2),               &
     &        IHQMAT(4,-1:1), IRQMAT(4, 0:7), IQDEF(4), ITOT, NZC1 ,    &
     &        I2_DGROUP
      COMMON /DGROUP/ I1_DGROUP, NZ, NFSYM, NBSYM, IFDIRP, JBTOF, JMROI,&
     &                JSPINR, JQBAS, JFSYM, IPQTOQ, IQTOPQ, IQPH,       &
     &                IHQMAT, IRQMAT, IQDEF, ITOT, NZC1, I2_DGROUP

!     Do NOT add any variables after the end tag !
      LOGICAL LINEAR
      COMMON /DGRPL/ LINEAR

!     ... NOPTYP is number of different operator types in DIRAC
      INTEGER NOPTYP
      PARAMETER (NOPTYP = 20)

! *** 4-component operators : Matrix part ***
!
!     M = I_4,iA_z,iA_y,iA_x
!
!     JM4REP(0:7) - bosonirrep of matrix operator
!     JM4POS(0:7) - quaternionic position of matrix operator 
!                   (1,2(i),3(j),4(k))
!     JM4FAS(0:7) - sign of matrix operator (+/- 1)
!
!     OPERATORS (NOPTYP types):
!
!     MCMP  - number of components for a given operator
!     JM4   - pointer to matrix operator for a given component
!     JCOM  - coefficient of component
!
      INTEGER MCMP(NOPTYP), JM4(3,NOPTYP), JCOM(3,NOPTYP),              &
     &         JM4REP(0:7),   JM4POS(0:7),     JM4FAS(0:7)
      COMMON /ONEOPS/ MCMP, JM4, JCOM, JM4REP, JM4POS, JM4FAS
! -- end of dgroup.h --
