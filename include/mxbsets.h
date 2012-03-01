! File: mxbsets.h
!
! MXBSETS = max number of regular basis sets
!   1) in Dalton typically 2 : normal, auxiliary (for r12)
!   2) in Dirac  typically 2 : large component, small component
!
! Basis set indices will go from -MXBSETS..+MXBSETS
!   Basis set -i will be density fitting basis set for basis +i
!   Basis set 0 (zero) is for Huckel basis set
!
!     IF you change any of these parameters you should do a "make depend"
!     and then rebuild the program using the command "make".
!
      INTEGER MXBSETS, MXBSETS_TOT
      PARAMETER (MXBSETS = 2, MXBSETS_TOT = 2*MXBSETS+1)
! -- end of mxbsets.h --
