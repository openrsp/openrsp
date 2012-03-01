#include <iratdef.h>
#include <maxorb.h>
      PARAMETER (LRECAO = 64*1024, LRINT=LRECAO)
      INTEGER*8 NSOINT
      COMMON /N_SORT/ RINT(LRINT),  NSOINT(106),
     &             LASTAD(106),NREC(106), NSOBAT(106),
     &             IT(MXCORB),IS(MXCORB),IBATCH(666),
     &             KEEP(8),NBAS(8),IBAS(8),NSYM,
     &             IDATA(18),LINT,LBATCH, INTSYM,IPRINT,IPRFIO,IW
      INTEGER IINT(IRAT*LRINT)
      EQUIVALENCE (RINT, IINT)
