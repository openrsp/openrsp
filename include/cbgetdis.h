C     -----------------------------------------------------------------
C     Purpose of /CBGETD/: transfer information to GETINT about
C     how the "H2AC" integrals are available.  The distribution
C     type DISTYP tells about packing and index symmetry.
C     IADINT .lt. 0 means integrals are in-core.
C     IADINT .ge. 0 is offset for reading an active-active distribution
C     from disk.  IADH2, IADH2X, IADH2D are used for off-sets for
C     H2AC, H2XAC, and H2DAC, resp. (H2DAC from ABACUS).
      INTEGER DISTYP
      COMMON /CBGETD/ DISTYP, IADINT,IADH2,IADH2X,IADH2D
C     -----------------------------------------------------------------
