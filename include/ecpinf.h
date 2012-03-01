C ==========================================================================
C >>>>>>>>>>>              ECP specifics                          <<<<<<<<<<
C ==========================================================================
C MXECP maximum number of ECP atoms
C MXEXP maximum number of exponents per ECP atoms
C MXANG maximum angular momentum quantum number
C NTYECP number of ECP atom types
C INDECP(I) pointer to first ECP atom of I
C NECP(I) number of ECP atoms of type I 
C 
C Vpp = sum_i ( r^(ncr(i)-2) ccr(i) exp(-zcr(i)*r^2) )      (local part)
C     + sum_l (sum_i ( r^(ncr(i)-2) ccr(i) exp(-zcr(i)*r^2) ) |l><l| ) 
C                                                        (nonlocal part)
C ncr(i) exponents n for r**(n-2)
C zcr(i) exponents
C ccr(i) coefficients
C lcr(k) maximal l value + 1 of the pseudopotential
C nkcrl(l,k) lower index of the pseudopotential k for projection operator l
C nkcru(l,k) upper index of the pseudopotential k for projection operator l
C
C KCRU pointer to ECP atom number
C LCRU is LMAX value for current ECP atom
C DFAC hardcoded numbers from subroutine CORTAB
C ==========================================================================
C
      LOGICAL ECP
      INTEGER MXECP, MXNONT, MXEXP, MXANG
      INTEGER NTYECP, NECP, INDECP
      INTEGER LCR, NCR, NKCRL, NKCRU
      DOUBLE PRECISION ZCR, CCR
C
      PARAMETER ( MXECP = 5, MXNONT = 80, MXEXP = 20, MXANG = 6 )
C
      COMMON /ECP_logical/ ECP
      COMMON /ECP_integer/ NTYECP, NECP(MXECP), INDECP(MXECP*MXNONT)
      COMMON /TRQ99/  LCR(MXECP), NCR(MXECP*MXEXP), NKCRL(MXANG,MXECP),
     &                NKCRU(MXANG,MXECP) 
      COMMON /TRQ98/  ZCR(MXECP*MXEXP), CCR(MXECP*MXEXP)
