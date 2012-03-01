!
!     File: aovec.h
!
!     MXAOVC = maximum number of AO input blocks for a given center
!              (in "MOLECULE.INP" input file)
!
!     IF you change MXAOVC you should do a "make depend"
!     and then rebuild the program using the command "make".
!
      INTEGER MXAOVC, MXCONT
      PARAMETER (MXAOVC = 25, MXCONT = 35)
!     hjaaj: MXCONT needs to be 35 for ANO-DK3 uncontracted for DKH Hamiltonian
