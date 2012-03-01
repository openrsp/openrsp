      INTEGER A, B, C, D, E
      PARAMETER (A=1, B=2, C=3, D=4, E=5)

      INTEGER AB, AC, AD, AE, BC, BD, BE, CD, CE, DE
      PARAMETER (AB=1, AC=2, AD=4, AE=7,
     &                 BC=3, BD=5, BE=8, 
     &                       CD=6, CE=9,
     &                             DE=10 )
 
      INTEGER ABC, ABD, ACD, BCD, ABE, ACE, ADE, BCE, BDE, CDE
      PARAMETER (ABC=1, ABD=2, ABE=5,
     &                  ACD=3, ACE=6,
     &                         BCE=7,
     &                  BCD=4, ADE=8,
     &                         BDE=9,
     &                         CDE=10  )


      INTEGER I1(15), I2(15), I3(15), I4(15), I5(15), IP1(15), IP2(15)

      DATA I1  / A, A, A,  B, B, B,  C, C, C,  D, D, D,  E, E, E  /

      DATA I2  / B, B, B,  A, A, A,  A, A, A,  A, A, A,  A, A, A  /
      DATA I3  / C, D, E,  C, D, E,  B, D, E,  B, C, E,  B, C, D  /
      DATA IP1 / BC,BD,BE, AC,AD,AE, AB,AD,AE, AB,AC,AE, AB,AC,AD /

      DATA I4  / D, C, C,  D, C, C,  D, B, B,  C, B, B,  C, B, B  /
      DATA I5  / E, E, D,  E, E, D,  E, E, D,  E, E, C,  D, D, C  /
      DATA IP2 / DE,CE,CD, DE,CE,CD, DE,BE,BD, CE,BE,BC, CD,BD,BC /


      INTEGER J1(5), J2(5), J3(5), J4(5), J5(5), JP1(5), JP2(5)

      DATA J1  / A, B, C, D, E  /
      DATA J2  / B, C, D, E, A  /
      DATA J3  / C, D, E, A, B  /
      DATA J4  / D, E, A, B, C  /
      DATA J5  / E, A, B, C, D  /
      DATA JP1 / AB,BC,CD,DE,AE /
      DATA JP2 / AC,BD,CE,AD,BE /

