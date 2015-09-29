#include <stdarg.h>

/* QcMatrix library */
#include "qcmatrix.h"

/* abstract OpenRSP context */
/*typedef struct _p_OpenRSP* OpenRSP;*/

/* tests struct*[] */
typedef struct {
    QInt value;
} TestStruct;

/* abstract OpenRSP context */
typedef struct OpenRSP OpenRSP;

struct OpenRSP {
    /* implementation-specific data */
    QVoid *data;
    /* implementation-specific OpenRSP APIs */
    QErrorCode (*openrspaddoneoper)(OpenRSP*,
                                    const QInt,
                                    const QInt*,
                                    const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                    QVoid*,
#endif
                                    va_list);
    QErrorCode (*openrspdestroy)(OpenRSP*);
};

/* OpenRSP APIs */
extern QErrorCode OpenRSPCreate(OpenRSP*);
extern QErrorCode OpenRSPAddOneOper(OpenRSP*,
                                    const QInt,
                                    const QInt*,
                                    const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                    QVoid*,
#endif
                                    ...);
extern QErrorCode OpenRSPDestroy(OpenRSP*);

