*======================================================================*
*    
*     Purpose: store information about effective Fock matrices
*              calculated from (single) one-index transformed
*              integrals:
*              1) label list maintained via function I1DXFCK
*              2) the start addresses on file for the fock matrices
*
*         L1DXFCKOPN    --  flag to allow/disallow extension for Fock m.
*         MAX1DXFCKLBL  --  maximum dimension of the list
*         N1DXFCKLBL    --  list length
*         LBL1DXFCK(*)  --  labels of the operators 
*         LSTRLX1DX(*)  --  list type of relaxation vector
*         IRELAX1DX(*)  --  list index of relaxation vector
*         LBLRLX(*)     --  label for relaxation vector
*         FRQRLX(*)     --  frequency for relaxation vector
*         ISYMRLX(*)    --  symmetry of relaxation vector
*         FIL1DXFCK     --  file name for fock matrix results
*         IADR1DXF(2,*) --  addresses of the fock matrix results on file
*                           1) CC effective Fock matrix (kappa+R,F)
*                           2) Fock matrix needed for kappa rhs vector
*
*
*======================================================================*
      INTEGER MAX1DXFCKLBL
      PARAMETER ( MAX1DXFCKLBL = 120 )
      CHARACTER*8 FIL1DXFCK
      PARAMETER (FIL1DXFCK = 'CC1DXFCK')

      LOGICAL L1DXFOPN
      INTEGER N1DXFLBL

      INTEGER IADR1DXF(2,MAX1DXFCKLBL)
      INTEGER IRELAX1DX(MAX1DXFCKLBL), ISYMRLX(MAX1DXFCKLBL)
      CHARACTER*8 LBL1DXFCK(MAX1DXFCKLBL), LABRLX(MAX1DXFCKLBL)
      CHARACTER*3 LST1DXFCK(MAX1DXFCKLBL)

#if defined (SYS_CRAY)
      REAL FRQRLX(MAX1DXFCKLBL)
#else
      DOUBLE PRECISION FRQRLX(MAX1DXFCKLBL)
#endif

      COMMON/I1DXFOCK/ IADR1DXF, IRELAX1DX, N1DXFLBL, ISYMRLX, L1DXFOPN
      COMMON/C1DXFOCK/ LBL1DXFCK, LST1DXFCK, LABRLX
      COMMON/R1DXFOCK/ FRQRLX

*======================================================================*
