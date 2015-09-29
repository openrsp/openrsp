#include "OpenRSP.h"

#define NUM_ALL_PERT 3
#define PERT_GEOMETRIC 1
#define PERT_DIPOLE 2
#define PERT_MAGNETIC 5
#define MAX_ORDER_GEOMETRIC 7
#define MAX_ORDER_DIPOLE 1
#define MAX_ORDER_MAGNETIC 7

QVoid get_one_oper_mat(const QInt num_pert,
                       const QInt *pert_labels,
                       const QInt *pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                       QVoid *user_ctx,
#endif
                       const QInt num_int,
                       QcMat *val_int[])
{
    printf("get_one_oper_mat %d\n", num_pert);
}

QVoid get_one_oper_exp(const QInt num_pert,
                       const QInt *pert_labels,
                       const QInt *pert_orders,
                       const QInt num_dmat,
                       QcMat *dens_mat[],
#if defined(OPENRSP_C_USER_CONTEXT)
                       QVoid *user_ctx,
#endif
                       const QInt num_exp,
                       QReal *val_exp)
{
    printf("get_one_oper_exp %d\n", num_pert);
}

QErrorCode main()
{
    OpenRSP open_rsp;
    QErrorCode ierr;                                  /* error information */

    QInt oneham_num_pert = 2;
    QInt oneham_pert_labels[2] = {PERT_GEOMETRIC,PERT_MAGNETIC};
    QInt oneham_pert_orders[2] = {MAX_ORDER_GEOMETRIC,MAX_ORDER_MAGNETIC};
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *oneham_context = "ONEHAM";
#endif

/* test TestStruct*[] */
    QInt num_ptr=2;
    TestStruct **ptr;
    TestStruct val[2];
    val[0].value = 12;
    val[1].value = 20;
    ptr = (TestStruct **)malloc(num_ptr*sizeof(TestStruct*));
    if (ptr==NULL) {
        QErrorExit(FILE_AND_LINE, "malloc ptr");
    }
    ptr[0] = &val[0];
    ptr[1] = &val[1];
    printf("value is %d (%d)\n", ptr[0]->value, 0);
    printf("value is %d (%d)\n", ptr[1]->value, 1);

printf("address %p\n", (QVoid *)ptr);
printf("address %p\n", (QVoid *)ptr[0]);
printf("address %p\n\n", (QVoid *)ptr[1]);
printf("address %p\n", (QVoid *)&ptr);
printf("address %p\n", (QVoid *)&ptr[0]);
printf("address %p\n", (QVoid *)&ptr[1]);

    ierr = OpenRSPCreate(&open_rsp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPCreate");

    /* adds one-electron Hamiltonian */
    ierr = OpenRSPAddOneOper(&open_rsp,
                             oneham_num_pert,
                             oneham_pert_labels,
                             oneham_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             (QVoid *)oneham_context,
#endif
                             &get_one_oper_mat,
                             &get_one_oper_exp,
/* test TestStruct*[] */
                             num_ptr,
                             ptr[0],
                             ptr);
    free(ptr);
    ptr=NULL;
    printf("value is %d (%d)\n", val[0].value, 0);
    printf("value is %d (%d)\n", val[1].value, 1);

    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPAddOneOper(h)");

    ierr = OpenRSPDestroy(&open_rsp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling OpenRSPDestroy");

    return QSUCCESS;
}

