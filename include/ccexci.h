*---------------------------------------------------------------------*
*
*     Purpose: storage of the excited states information
*
*     MAXEXCI : maximum dimension
*     NEXCI   : total number of states requested by input
*
*     the excited states energies and vectors are sorted after
*           1) symmetry classes
*           2) multiplicities
*           3) energy (lowest eigenvalues first)
*     (see CC_EXCI for more information)
*
*     ISYEXC  : the symmetries of the states
*     IMULT   : the multiplicities (singlett/triplett) of the states
*     EIGVAL  : eigenvalues (excitation energies)
*     XNORM   : overlap between left and right vector
*               (why is this stored on a common block? should be
*                1.00 after convergence of the left vectors...)
*       
*     ISYOFE  : offsets for the symmetry classes (singlett states)
*     ITROFE  : offsets for the triplett states for each symmetry
*
* Christof Haettig, Januar 1999
*---------------------------------------------------------------------*
      INTEGER MAXEXCI
      PARAMETER ( MAXEXCI = 100)

      INTEGER NEXCI
      INTEGER ISYEXC(MAXEXCI),IMULTE(MAXEXCI)
      INTEGER ISYOFE(8), ITROFE(8)

#if defined (SYS_CRAY)
      REAL EIGVAL(MAXEXCI),XNORM(MAXEXCI)
#else
      DOUBLE PRECISION EIGVAL(MAXEXCI),XNORM(MAXEXCI)
#endif

      COMMON/IEXCI  / NEXCI,ISYEXC,IMULTE,ISYOFE,ITROFE
      COMMON/EEXCI  / EIGVAL,XNORM
*---------------------------------------------------------------------*
