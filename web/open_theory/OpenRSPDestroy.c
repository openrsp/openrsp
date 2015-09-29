#include "OpenRSP.h"

QErrorCode OpenRSPDestroy(OpenRSP *open_rsp)
{
    QErrorCode ierr;

    ierr = open_rsp->openrspdestroy(open_rsp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling openrspdestroy");

    if (open_rsp->data!=NULL) {
        free(open_rsp->data);
        open_rsp->data = NULL;
    }

    open_rsp->openrspaddoneoper = NULL;
    open_rsp->openrspdestroy = NULL;

    return QSUCCESS;
}

