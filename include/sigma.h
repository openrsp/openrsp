C     used for transfer of NMR shieldings between routines (SIGMA)
      LOGICAL         DONS 
      REAL*8          SIGMAD, SIGMAS, SIGMAR, SIGMAT,
     &                SIGMADZ, SIGMASFTP, SIGMASFTM
      COMMON /SIGMA/  SIGMAD(3,MXCOOR), SIGMAS(3,MXCOOR),
     &                SIGMAR(3,MXCOOR), SIGMAT(3,MXCOOR),
     &                SIGMADZ(3,MXCOOR),SIGMASFTP(3,MXCOOR),
     &                SIGMASFTM(3,MXCOOR), DONS
