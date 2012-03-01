!
!     file: mxcent.h
!
!     MXCENT_QM = max number of QM nuclei + point charges + ghost orbital centers
!     MXCENT = MXCENT_QM + max number of MM nuclei
!
!     IF you change MXCENT you should do a "make depend"
!     and then rebuild the program using the command "make".
!
!     In case of QM3 MXCENT_QM is used to allocate memory in herrdn.F.
!     To run a QM3 calculation in most cases MXCENT will have to be
!     around 2000 - 3000!!! Remember to set MXQM3 = MXCENT_QM in qm3.h!!!
!
      INTEGER MXCENT_QM, MXCENT, MXCOOR
      PARAMETER (MXCENT_QM =120, MXCENT = 120, MXCOOR = 3*MXCENT)
! -- end of mxcent.h --
