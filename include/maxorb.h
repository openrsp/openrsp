!
!     File: maxorb.h
!
!     MXSHEL = maximum number of shells (insert shell definition here).
!              (if modified: also  change MXSHEL for __CVERSION__ in infpar.h)
!     MXPRIM = maximum number of primitives.
!     MXCORB = maximum number of orbitals (possibly contracted).
!     MAXOCC = maximum number of occupied orbitals
!
!     IF you change any of these parameters you should do a "make depend"
!     and then rebuild the program using the command "make".
!
      INTEGER MXSHEL, MXPRIM, MXCORB, MXORBT, MAXOCC
      PARAMETER (MXSHEL = 1000, MXPRIM = 8000, MXCORB = 1200,           &
     &           MAXOCC = 400, MXORBT = MXCORB*(MXCORB + 1)/2)
