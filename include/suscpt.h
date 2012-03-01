C     used for transfer of magn. susceptibility between routines
      LOGICAL DONETWOLOP
      REAL*8  SUSDIA, SUS2EL, SUSFS,  SUSREL,
     &        SUSTOT, SUSFSY, SUSCOM, SUSDZD,
     &        SUSDFT
      COMMON /SUSCPT/ SUSDIA(3,3), SUS2EL(3,3), SUSFS(3,3), SUSREL(3,3),
     &                SUSTOT(3,3), SUSFSY(3,3), SUSCOM(3,3),SUSDZD(3,3),
     &                SUSDFT(3,3), DONETWOLOP
