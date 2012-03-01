! FILE : aosotr.h
!     AO atomic orbital, SO = symmetry orbital
!     CTRAN, ITRAN, JTRAN: transformation from SO to AO
!     IPTYP, IPCEN: gto type and sym.indep. center index for each SO
!     IAOAO, JAOAO: AO to AO transformation (see AOTOAO subroutine in herrdn.F)
!     IAOINFO: information about each atomic orbital
!        IAOINFO(:,1) =   sym. indep. AO index
!        IAOINFO(:,2) =   sym. dep. atom index (e.g. to NAMDEP)
!        IAOINFO(:,3) =   orbital type index

      REAL*8  CTRAN(MXCORB,8)
      INTEGER ITRAN(MXCORB,8),JTRAN(MXCORB), IPTYP(MXCORB),IPCEN(MXCORB)&
     &      , IAOAO(MXCORB), JAOAO(MXCORB,2), IAOINFO(MXCORB,3)
      COMMON /AOSOTR/ CTRAN, ITRAN, JTRAN, IPTYP, IPCEN, IAOAO, JAOAO,  &
     &        IAOINFO
! -- end of aosotr.h
