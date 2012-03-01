!  FILE : incore.h
!
!  Purpose : storage and control variables for "AO integrals in core" feature
!
C     PARAMETER (MXTSK=(MXSHEL*(MXSHEL+1)/2), IEIR=0, MMWORK=0,
C    &     MXSHL_KEEP=MAX(MXTSK,IEIR), MNSHL_KEEP=MIN(MXTSK,IEIR) )
Chj 11-May-2005: changes needed for compilation with g77:
C     MAX and MIN not allowed by "g77"; zero size arrays not allowed
      PARAMETER (MXTSK=(MXSHEL*(MXSHEL + 1)/2),
     &           MMWORK=INSTALL_MMWORK )

      REAL*8  AOINTSCORE
      COMMON /AOINTS_INCORE/ AOINTSCORE(MMWORK)

      LOGICAL AOSAVE, LINTSV, INITX, LINTMP, MSAVE
      COMMON /AOINTS_INCORE_2/
     &     ISCORE, MMCORE, LMCORE, JSCORE, ITOTNT,
     &     INDX_SHL1, INDX_SHL2, INDX_SHL3, INDX_SHL4,
     &     I_SHL, N_SHL,
     &     IPREVA, IPREVB, IPREVC, IPREVD,
     &     INDX_SHL(MXTSK), INDX_C(3,MXTSK),
     &     AOSAVE, LINTSV, INITX, LINTMP, MSAVE
#if 0
#if defined (VAR_MPI)
      LOGICAL INDMPI
      COMMON /MPISAV/ INDMPI(MXTSK)
#endif

/* MMWORK and IEIR must be set to nonzero values here if the 
  AOSAVE-feature of DALTON is to be used. Suggested values: 
  MMWORK=80000000 and IEIR=3000000, 
  if your hardware can take it.*/
#endif
! -- end of incore.h --
