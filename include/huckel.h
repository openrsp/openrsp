C
C     File: huckel.h
C
C  HUCCNT    : Huckel factor, normally 1.75
C  HUCEXC(i) : atomic orbital energies for Huckel SOs organized after
C              symmetry , i=1..NHUCBA
C  NHUCAO(8) : number of Huckel SOs per symmetry
C  NHUCBA    : total number of Huckel SOs (e.g. 3 p-functions)
C  NHUCOCC(8): Huckel orbital occupation for NSYM symmetries 
C  IHUCPT(j) : pointer to Huckel AO "i" from total basis fu list
C              (used in her1*.F when Huckel matrix is constructed)
C              j=1..KMAX (KMAX is total number of symmetry independent
C                         AO shells, )
C  ISETHUCKEL : which basis set number is Huckel basis set
C  EWMO      : do EWMO instead of EHT
C
      REAL*8  HUCCNT, HUCEXC(MXSHEL)
      INTEGER NHUCAO(8), NHUCBA, IHUCPT(MXSHEL), IPRHUC, ISETHUCKEL
      LOGICAL DOHUCKEL, EWMO, HUCPROJCMO
C
      COMMON /CBHUCKEL/ HUCCNT, HUCEXC,
     &       NHUCAO, NHUCBA, IHUCPT, IPRHUC, ISETHUCKEL,
     &                DOHUCKEL, EWMO,HUCPROJCMO
C --- end of huckel.h ---
