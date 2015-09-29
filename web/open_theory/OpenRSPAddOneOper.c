#include "OpenRSP.h"

QErrorCode OpenRSPAddOneOper(OpenRSP *open_rsp,
                             const QInt num_pert,
                             const QInt *pert_labels,
                             const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             QVoid *user_ctx,
#endif
                             ...)
{
    va_list callback_one_oper;
#if defined(OPENRSP_C_USER_CONTEXT)
    va_start(callback_one_oper, user_ctx)
#else
    va_start(callback_one_oper, pert_max_orders);
#endif
    open_rsp->openrspaddoneoper(open_rsp,
                                num_pert,
                                pert_labels,
                                pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                user_ctx,
#endif
                                callback_one_oper);
    va_end(callback_one_oper);
    return QSUCCESS;
}
