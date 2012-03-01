      INTEGER A, B, C, D
      PARAMETER (A=1, B=2, C=3, D=4)
      INTEGER AB, AC, AD, BC, BD, CD
      PARAMETER (AB=1, AC=2, AD=3, BC=4, BD=5, CD=6)
      INTEGER ABC, ABD, ACD, BCD
      PARAMETER (ABC=1, ABD = 2, ACD = 3, BCD = 4)

      INTEGER I1(6), I2(6), I3(6), I4(6), IP(6), IP1(3), IP2(3)
      DATA I1  / A, A, A, B, B, C /
      DATA I2  / B, C, D, C, D, D /
      DATA I3  / C, B, B, A, A, A /
      DATA I4  / D, D, C, D, C, B /
      DATA IP  / CD, BD, BC, AD, AC, AB /
      DATA IP1 / AB, AC, AD /
      DATA IP2 / CD, BD, BC /
