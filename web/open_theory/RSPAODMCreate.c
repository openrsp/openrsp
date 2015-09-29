#include "OpenRSPAODM.h"

QErrorCode RSPAODMCreate(OpenRSP *open_rsp)
{
    RSPAODM *rsp_ao_dm;

    rsp_ao_dm = (RSPAODM *)malloc(sizeof(RSPAODM));
    if (rsp_ao_dm==NULL) {
        QErrorExit(FILE_AND_LINE, "malloc rsp_ao_dm");
    }
    rsp_ao_dm->one_oper = NULL;
    open_rsp->data = (QVoid *)rsp_ao_dm;
    open_rsp->openrspaddoneoper = RSPAODMAddOneOper;
    open_rsp->openrspdestroy = RSPAODMDestroy;

    return QSUCCESS;
}
