      INTEGER A, B, C, D, E, F
      PARAMETER (A=1, B=2, C=3, D=4, E=5, F=6)

      INTEGER AB, AC, AD, AE, AF, BC, BD, BE, BF, CD, CE, CF, DE, DF,EF
      PARAMETER (AB=1, AC=2, AD=4, AE=7,  AF=11,
     &                 BC=3, BD=5, BE=8,  BF=12,
     &                       CD=6, CE=9,  CF=13, 
     &                             DE=10, DF=14,
     &                                    EF=15 )

      INTEGER ABC, ABD, ABE, ABF, ACD, ACE, ACF, ADE, ADF, AEF,
     &                            BCD, BCE, BCF, BDE, BDF, BEF,
     &                                           CDE, CDF, CEF,
     &                                                     DEF   
      PARAMETER (ABC=1, ABD=2, ABE=5,  ABF=11, 
     &                  ACD=3, ACE=6,  ACF=12, 
     &                         BCE=7,  BCF=13, 
     &                                 ADF=14,
     &                  BCD=4, ADE=8,  BDF=15, 
     &                         BDE=9,  CDF=16, 
     &                                 AEF=17,
     &                         CDE=10, BEF=18, 
     &                                 CEF=19,
     &                                 DEF=20  )

      INTEGER I1(15), I2(15), I3(15), I4(15), I5(15), I6(15)
      INTEGER IP1(15), IP2(15), IP3(15)

      DATA I1  / A, A, A,  A, A, A,  A, A, A,  A, A, A,  A, A, A  /
      DATA I2  / B, C, D,  B, C, E,  B, E, D,  E, C, D,  F, F, F  /
      DATA IP1 / AB,AC,AD, AB,AC,AE, AB,AE,AD, AE,AC,AD, AF,AF,AF /

      DATA I3  / C, B, B,  C, B, B,  C, B, B,  B, B, B,  B, B, B  /
      DATA I4  / D, D, C,  E, E, C,  F, D, E,  F, F, F,  E, D, C  /
      DATA IP2 / CD,BD,BC, CE,BE,BC, CF,BD,BE, BF,BF,BF, BE,BD,BC /

      DATA I5  / E, E, E,  D, D, D,  D, C, C,  C, D, C,  C, C, D  /
      DATA I6  / F, F, F,  F, F, F,  E, F, F,  D, E, E,  D, E, E  /
      DATA IP3 / EF,EF,EF, DF,DF,DF, DE,CF,CF, CD,DE,CE, CD,CE,DE /



      INTEGER K1(45), K2(45), KP1(45), KP2(45), KP3(45)

      DATA K1  / A, A, A,  A, A, A,  A, A, A,  A, A, A,  A, A, A,   !I1
     &           C, B, B,  C, B, B,  C, B, B,  B, B, B,  B, B, B,   !I3
     &           E, E, E,  D, D, D,  D, C, C,  C, D, C,  C, C, D  / !I5

      DATA K2  / B, C, D,  B, C, E,  B, E, D,  E, C, D,  F, F, F,   !I2
     &           D, D, C,  E, E, C,  F, D, E,  F, F, F,  E, D, C,   !I4
     &           F, F, F,  F, F, F,  E, F, F,  D, E, E,  D, E, E  / !I6

      DATA KP1 / AB,AC,AD, AB,AC,AE, AB,AE,AD, AE,AC,AD, AF,AF,AF, !IP1
     &           CD,BD,BC, CE,BE,BC, CF,BD,BE, BF,BF,BF, BE,BD,BC, !IP2
     &           EF,EF,EF, DF,DF,DF, DE,CF,CF, CD,DE,CE, CD,CE,DE /!IP3

      DATA KP2 / CD,BD,BC, CE,BE,BC, CF,BD,BE, BF,BF,BF, BE,BD,BC, !IP2
     &           EF,EF,EF, DF,DF,DF, DE,CF,CF, CD,DE,CE, CD,CE,DE, !IP3
     &           AB,AC,AD, AB,AC,AE, AB,AE,AD, AE,AC,AD, AF,AF,AF /!IP1

      DATA KP3 / EF,EF,EF, DF,DF,DF, DE,CF,CF, CD,DE,CE, CD,CE,DE, !IP3
     &           AB,AC,AD, AB,AC,AE, AB,AE,AD, AE,AC,AD, AF,AF,AF, !IP1
     &           CD,BD,BC, CE,BE,BC, CF,BD,BE, BF,BF,BF, BE,BD,BC /!IP2


      INTEGER IT1(10), IT2(10)
      DATA IT1 / ABC, ABD, ACD, BCD, ABE, ACE, BCE, ADE, BDE, CDE /
      DATA IT2 / DEF, CEF, BEF, AEF, CDF, BDF, ADF, BCF, ACF, ABF /


      INTEGER JP(15), J1(15), J2(15), J3(15), J4(15), J5(15)
      DATA JP /AB,AC,BC,AD,BD,CD,AE,BE,CE,DE,AF,BF,CF,DF,EF/
      DATA J1 / C, B, A, B, A, A, B, A, A, A, B, A, A, A, A/
      DATA J2 / D, D, D, C, C, B, C, C, B, B, C, C, B, B, B/
      DATA J3 / E, E, E, E, E, E, D, D, D, C, D, D, D, C, C/
      DATA J4 / F, F, F, F, F, F, F, F, F, F, E, E, E, E, D/


