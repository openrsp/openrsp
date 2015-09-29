#include "OpenRSPAODM.h"

QErrorCode RSPAODMAddOneOper(OpenRSP *open_rsp,
                             const QInt num_pert,
                             const QInt *pert_labels,
                             const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             QVoid *user_ctx,
#endif
                             va_list callback_one_oper)
{
    RSPAODM *rsp_ao_dm;
    va_list ao_dm_callback;
    QInt ipert;
    QInt jpert;

    QcMat *val_int[2];
    QReal *val_exp;
/* test TestStruct*[] */
    QInt num_ptr;
    TestStruct **ptr;
    TestStruct **tmp;

    rsp_ao_dm = (RSPAODM *)(open_rsp->data);
    rsp_ao_dm->one_oper = (AODMOneOper *)malloc(sizeof(AODMOneOper));
    if (rsp_ao_dm->one_oper==NULL) {
        QErrorExit(FILE_AND_LINE, "malloc rsp_ao_dm->one_oper");
    }

    if (num_pert>0) {
        rsp_ao_dm->one_oper->num_pert = num_pert;
    }
    else {
        printf("RSPAODMAddOneOper>> number of perturbations %"QINT_FMT"\n", num_pert);
        QErrorExit(FILE_AND_LINE, "invalid number of perturbations");
    }
    rsp_ao_dm->one_oper->pert_labels = (QInt *)malloc(num_pert*sizeof(QInt));
    if (rsp_ao_dm->one_oper->pert_labels==NULL) {
        printf("RSPAODMAddOneOper>> number of perturbations %"QINT_FMT"\n", num_pert);
        QErrorExit(FILE_AND_LINE, "malloc rsp_ao_dm->one_oper->pert_labels");
    }
    rsp_ao_dm->one_oper->pert_max_orders = (QInt *)malloc(num_pert*sizeof(QInt));
    if (rsp_ao_dm->one_oper->pert_max_orders==NULL) {
        printf("RSPAODMAddOneOper>> number of perturbations %"QINT_FMT"\n", num_pert);
        QErrorExit(FILE_AND_LINE, "malloc rsp_ao_dm->one_oper->pert_max_orders");
    }
    for (ipert=0; ipert<num_pert; ipert++) {
        /* each element of \var{pert_labels} should be unique */
        for (jpert=0; jpert<ipert; jpert++) {
            if (pert_labels[jpert]==pert_labels[ipert]) {
                printf("RSPAODMAddOneOper>> perturbation %"QINT_FMT" is %"QINT_FMT"\n",
                       jpert,
                       pert_labels[jpert]);
                printf("RSPAODMAddOneOper>> perturbation %"QINT_FMT" is %"QINT_FMT"\n",
                       ipert,
                       pert_labels[ipert]);
                QErrorExit(FILE_AND_LINE, "same perturbation not allowed");
            }
        }
        rsp_ao_dm->one_oper->pert_labels[ipert] = pert_labels[ipert];
        if (pert_max_orders[ipert]<1) {
            printf("RSPAODMAddOneOper>> order of %"QINT_FMT"-th perturbation (%"QINT_FMT") is %"QINT_FMT"\n",
                   ipert,
                   pert_labels[ipert],
                   pert_max_orders[ipert]);
            QErrorExit(FILE_AND_LINE, "only positive order allowed");
        }
        rsp_ao_dm->one_oper->pert_max_orders[ipert] = pert_max_orders[ipert];
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    rsp_ao_dm->one_oper->user_ctx = user_ctx;
#endif

/*
https://github.com/dotnet/coreclr/issues/329
src/pal/src/cruntime/silent_printf.cpp: Silent_PAL_vfprintf()

As per C11 standard (http://www.open-std.org/JTC1/SC22/WG14/www/docs/n1570.pdf):

If access to the varying arguments is desired, the called function shall
declare an object (generally referred to as ap in this subclause) having type
va_list. The object ap may be passed as an argument to another function; if
that function invokes the va_arg macro with parameter ap, the value of ap in
the calling function is indeterminate and shall be passed to the va_end macro
prior to any further reference to ap.

It is permitted to create a pointer to a va_list and pass that pointer to
another function, in which case the original function may make further use of
the original list after the other function returns."
*/
    va_copy(ao_dm_callback, callback_one_oper);
    rsp_ao_dm->one_oper->get_one_oper_mat = va_arg(ao_dm_callback, AODMOneOperMat);
    rsp_ao_dm->one_oper->get_one_oper_exp = va_arg(ao_dm_callback, AODMOneOperExp);
/* test TestStruct*[] */
    num_ptr = va_arg(ao_dm_callback, QInt);
    ptr = (TestStruct **)malloc(num_ptr*sizeof(TestStruct*));
    if (ptr==NULL) {
        QErrorExit(FILE_AND_LINE, "malloc ptr");
    }
    ptr[0] = va_arg(ao_dm_callback, TestStruct*);
    printf("address is %p\n", ptr[0]);
    printf("value is %d\n", ptr[0]->value);
    ptr[0] = NULL;
    tmp = va_arg(ao_dm_callback, TestStruct**);
    printf("address is %p\n", tmp);
    for (ipert=0; ipert<num_ptr; ipert++,tmp++) {
        ptr[ipert] = *tmp;
        printf("address is %p, %p (%d)\n", ptr[ipert], &ptr[ipert], ipert);
        printf("value is %d (%d)\n", ptr[ipert]->value, ipert);
    }
    ptr[0]->value = 33;
    ptr[1]->value = 44;
    free(ptr);
    ptr=NULL;

    va_end(ao_dm_callback);
    rsp_ao_dm->one_oper->next_oper = NULL;


    rsp_ao_dm->one_oper->get_one_oper_mat(num_pert,
                                          pert_labels,
                                          pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                          user_ctx,
#endif
                                          2,
                                          val_int);

    rsp_ao_dm->one_oper->get_one_oper_exp(num_pert,
                                          pert_labels,
                                          pert_max_orders,
                                          2,
                                          val_int,
#if defined(OPENRSP_C_USER_CONTEXT)
                                          user_ctx,
#endif
                                          2,
                                          val_exp);

    return QSUCCESS;
}
