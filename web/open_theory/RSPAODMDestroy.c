#include "OpenRSPAODM.h"

QErrorCode RSPAODMDestroy(OpenRSP *open_rsp)
{
    RSPAODM *rsp_ao_dm;

    QcMat *val_int[2];
    QReal *val_exp;

    rsp_ao_dm = (RSPAODM *)(open_rsp->data);
    if (rsp_ao_dm!=NULL) {
        if (rsp_ao_dm->one_oper!=NULL) {


    rsp_ao_dm->one_oper->get_one_oper_mat(rsp_ao_dm->one_oper->num_pert,
                                          rsp_ao_dm->one_oper->pert_labels,
                                          rsp_ao_dm->one_oper->pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                          rsp_ao_dm->one_oper->user_ctx,
#endif
                                          2,
                                          val_int);

    rsp_ao_dm->one_oper->get_one_oper_exp(rsp_ao_dm->one_oper->num_pert,
                                          rsp_ao_dm->one_oper->pert_labels,
                                          rsp_ao_dm->one_oper->pert_max_orders,
                                          2,
                                          val_int,
#if defined(OPENRSP_C_USER_CONTEXT)
                                          rsp_ao_dm->one_oper->user_ctx,
#endif
                                          2,
                                          val_exp);


            free(rsp_ao_dm->one_oper->pert_labels);
            rsp_ao_dm->one_oper->pert_labels = NULL;
            free(rsp_ao_dm->one_oper->pert_max_orders);
            rsp_ao_dm->one_oper->pert_max_orders = NULL;
#if defined(OPENRSP_C_USER_CONTEXT)
            rsp_ao_dm->one_oper->user_ctx = NULL;
#endif
            rsp_ao_dm->one_oper->get_one_oper_mat = NULL;
            rsp_ao_dm->one_oper->get_one_oper_exp = NULL;
            rsp_ao_dm->one_oper->next_oper = NULL;
            free(rsp_ao_dm->one_oper);
            rsp_ao_dm->one_oper = NULL;
        }
    }

    return QSUCCESS;
}
