C
C     Variables named ?1/2/3??? are used for triplet properties
C     calculations. Others are for singlet properties.
C     To reduce static memory allocations, all large matrices are allocated
C     when needed in rspsoppa.F, K.Ruud-96
C
C      INTEGER AB1ADR,AB2ADR,AB3ADR
      LOGICAL A2EXIST
      COMMON /INFSOP/ NSOO(8), NTOO(8), NSVV(8), NTVV(8),
     *                NS2P2H(8), NT2P2H(8), N2P2HT(8), N2P2HS(8),
     *                N12P2H(8), N22P2H(8), N32P2H(8), A2EXIST,
     *                KAB2,KABSAD,KABTAD,KIJSAD,KIJTAD,KIJ1AD,
     *                KIJ2AD,KIJ3AD,KIADR1
C
