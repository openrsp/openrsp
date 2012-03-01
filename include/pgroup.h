!  FILE : pgroup.h
!  Point group information
!
!     SYMOP - Symmetry operations
!     ===========================
!
!         zyx
!     0   000    E   identity
!     1   001   Oyz  reflection in the yz-plane
!     2   010   Oxz  reflection in the xz-plane
!     3   011   C2z  rotation about the z-axis
!     4   100   Oxy  reflection in the xy-plane
!     5   101   C2y  rotation about the y-axis
!     6   110   C2x  rotation about the x-axis
!     7   111    i   inversion centre
!
      CHARACTER*3 GROUP, REP(0:7)
      COMMON /PGROUP/ GROUP, REP
!
      CHARACTER*3 SYMOP(0:7), GROUPS(0:7)
      SAVE SYMOP, GROUPS
      DATA SYMOP / ' E ','Oyz','Oxz','C2z',                             &
     &             'Oxy','C2y','C2x',' i '/
      DATA GROUPS/ 'C1 ','C2 ','Cs ','Ci ',                             &
     &             'D2 ','C2v','C2h','D2h'/
! -- end of pgroup.h --
