      DIMENSION   IDBTYP(4,4)
      CHARACTER*9 COBTYP(4)
      SAVE      IDBTYP, COBTYP
      DATA      IDBTYP/-1,-2,-3,-4,
     &                 -2, 1, 2, 4,
     &                 -3, 2, 3, 5,
     &                 -4, 4, 5, 6/
      DATA      COBTYP/'frozen   ','inactive ','active   ','secondary'/
