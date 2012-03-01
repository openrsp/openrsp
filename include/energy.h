#if defined (SYS_CRAY)
      REAL ENERKE, ENERNA, ENEREE, ENERNN
      REAL GRADKE, GRADNA, GRADEE, GRADNN, GRADFS
      REAL GRADFT
#else
      DOUBLE PRECISION ENERKE, ENERNA, ENEREE, ENERNN
      DOUBLE PRECISION GRADKE, GRADNA, GRADEE, GRADNN, GRADFS
      DOUBLE PRECISION GRADFT
#endif
      COMMON /ENERGY/ ENERKE, GRADKE(MXCOOR),                            &
     &                ENERNA, GRADNA(MXCOOR),                            &
     &                ENEREE, GRADEE(MXCOOR),                            &
     &                ENERNN, GRADNN(MXCOOR),                            &
     &                        GRADFS(MXCOOR),                            &
     &                        GRADFT(MXCOOR)
