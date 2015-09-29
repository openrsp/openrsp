#include "OpenRSPAODM.h"

QErrorCode OpenRSPCreate(OpenRSP *open_rsp)
{
    QErrorCode ierr;

    ierr = RSPAODMCreate(open_rsp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPAODMCreate");

    return QSUCCESS;
}

