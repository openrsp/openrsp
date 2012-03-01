*======================================================================*
*    
*     Purpose: store information about expectation values and 
*              effective Fock matrices:
*              1) label list initialized by CC_EXPFCK_INIT and
*                 maintained via functions IEXPECT/IEFFFOCK
*              2) results for expectation values and the file name
*                 and the start addresses on file for the fock matrices
*
*         LEFCKOPN      --  flag to allow/disallow extension for Fock m.
*         LEXPTOPN      --  flag to allow/disallow extension for expect.
*         MAXEXPFCKLBL  --  maximum dimension of the list
*         NEXPFCKLBL    --  list length
*         LBLEXPFCK(*)  --  labels of the operators 
*         ORDEXPFCK(*)  --  derivative order of the operator
*         LEXPFCK(2,*)  --  flags for expt. values and fock matrices
*         ISYEFCK(*)    --  symmetry of the results
*         EXPVALUE(4,*) --  results for the expectation values:
*                           1) one-electron contrib. to CC exp. value
*                           2) two-electron contrib. to CC exp. value
*                           3) one-electron contrib. to SCF exp. value
*                           4) two-electron contrib. to SCF exp. value
*         FILFCKEFF     --  file name for fock matrix results
*         IADRFCK(2,*)  --  addresses of the fock matrix results on file
*                           1) CC effective Fock matrix
*                           2) Fock matrix needed for kappa rhs vector
*
*
*======================================================================*
      INTEGER MAXEXPFCKLBL
      PARAMETER ( MAXEXPFCKLBL = 120 )
      CHARACTER*8 FILFCKEFF
      PARAMETER (FILFCKEFF = 'CCEFFOCK')

      LOGICAL LEFCKOPN, LEXPTOPN
      INTEGER NEXPFCKLBL

      LOGICAL LEXPFCK(2,MAXEXPFCKLBL)
      INTEGER ORDEXPFCK(MAXEXPFCKLBL)
      INTEGER ISYEXPFCK(MAXEXPFCKLBL)
      INTEGER IADRFCK(2,MAXEXPFCKLBL)
      CHARACTER*8 LBLEXPFCK(MAXEXPFCKLBL)

#if defined (SYS_CRAY)
      REAL EXPVALUE(2,MAXEXPFCKLBL)
#else
      DOUBLE PRECISION EXPVALUE(4,MAXEXPFCKLBL)
#endif

      COMMON/IEXPFCK/ ISYEXPFCK, ORDEXPFCK, IADRFCK, NEXPFCKLBL, 
     &                LEXPFCK, LEFCKOPN, LEXPTOPN
      COMMON/CEXPFCK/ LBLEXPFCK
      COMMON/REXPFCK/ EXPVALUE

*======================================================================*
