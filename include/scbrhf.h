C scbrhf.h: Sirius Common Block for RHF (in fact for SCF - both RHF and DFT)
C
C Common block for SCF information, defined in SCFINP input module /hjaaj Oct 2003
C (called scbrhf.h because it was made before DFT was implemented)
C
C THRRHF: convergence threshold for SCF gradient
C NFRRHF: frozen SCF orbitals
C NMVO* : see input option .FC MVO
C IOPRHF: for 1 open shell RHF
C MXHFMA, MXHFMI: max macro and micro iterations in QC-SCF
C MAXFCK: max Roothaan Fock interations
C MXDIIS: max DIIS iterations
C MAXEVC: max error vectors in DIIS
C NRHFEL: number of electrons in SCF
C INIOCC: how is initial orbital occupation determined
C RHFCAN: generate canonical SCF orbitalsl, diagonalizing Fock matrix
C AUTOCC: use automatic orbital occupation
C BCKSTP: allow for backstepping to previous Fock matrix if DIIS stalls
C in infinp.h:
C HSROHF: High spin open-shell RHF
C
      REAL*8  THRRHF
      LOGICAL RHFCAN, AUTOCC, BCKSTP
      COMMON /SCBRHF/ THRRHF,
     &                NFRRHF(8),NMVO(8),NMVOT,IOPRHF,MXHFMA,MXHFMI,
     &                MAXFCK, MXDIIS, MAXEVC, NRHFEL,INIOCC,
     &                RHFCAN, AUTOCC, BCKSTP
C .. end of scbrhf.h ...
