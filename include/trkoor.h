C
C     NCOOR  - total number of Cartesian coordinates
C     NCOORS - number of Cartesian coordinates of each symmetry
C     DEPEND - true for dependent symmetry coordinates
C     NCDEP  - total number of dependent/trarot coordinates
C     NCDEPS - number of dependent/trarot coordinates of each symmetry
C     IPTTRO - identifies symmetry-ordered trarot coordinate as
C              Tx, Ty, Tz, Rx, Ry, or Rz
C
      LOGICAL DEPEND
      COMMON /TRKOOR/ NCOOR, NCOORS(8), DEPEND(MXCOOR),
     *                NCDEP, NCDEPS(8), NCIND, NCINDS(8),
     *                IPTTRO(6,8), NTRREP(0:7)
