#include "OpenRSP.h"

/* abstract OpenRSP context */
typedef struct _p_OpenRSP _p_OpenRSP;

struct _p_OpenRSP {
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

