/*
  OpenRSP: open-ended library for response theory
  Copyright 2015 Radovan Bast,
                 Daniel H. Friese,
                 Bin Gao,
                 Dan J. Jonsson,
                 Magnus Ringholm,
                 Kenneth Ruud,
                 Andreas Thorvaldsen
  This source code form is subject to the terms of the
  GNU Lesser General Public License, version 2.1.
  If a copy of the GNU LGPL v2.1 was not distributed with this
  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

*/

#include "RSPZeroOper.h"

/* <function name='RSPZeroOperCreate'
             attr='private'
             author='Bin Gao'
             date='2015-02-12'>
     Create a node of a linked list for a given zero-electron operator, should
     be called at first
     <param name='zero_oper' direction='inout'>
       The linked list of zero-electron operators
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act as
       perturbations on the zero-electron operator
     </param>
     <param name='pert_labels' direction='in'>
       All the different perturbation labels
     </param>
     <param name='pert_max_orders' direction='in'>
       Allowed maximal order of a perturbation described by exactly one of
       the above different labels
     </param>
     <param name='user_ctx' direction='in'>
       User-defined callback-function context
     </param>
     <param name='get_zero_oper_contrib' direction='in'>
       User-specified function for calculating contribution of the
       zero-electron operator and its derivatives
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPZeroOperCreate(RSPZeroOper **zero_oper,
                             const QInt num_pert_lab,
                             const QcPertInt *pert_labels,
                             const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             void *user_ctx,
#endif
                             const GetZeroOperContrib get_zero_oper_contrib)
{
    RSPZeroOper *new_oper;  /* new operator */
    QInt ilab;  /* incremental recorders over perturbation labels */
    QInt jlab;
    new_oper = (RSPZeroOper *)malloc(sizeof(RSPZeroOper));
    if (new_oper==NULL) {
        QErrorExit(FILE_AND_LINE, "allocates memory for zero-electron operator");
    }
    if (num_pert_lab<0) {
        printf("RSPZeroOperCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        QErrorExit(FILE_AND_LINE, "invalid number of perturbation labels");
    }
    else if (num_pert_lab>OPENRSP_PERT_LABEL_MAX) {
        printf("RSPZeroOperCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        printf("RSPZeroOperCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
               OPENRSP_PERT_LABEL_MAX);
        QErrorExit(FILE_AND_LINE, "too many perturbation labels");
    }
    new_oper->num_pert_lab = num_pert_lab;
    if (new_oper->num_pert_lab>0) {
        new_oper->pert_max_orders = (QInt *)malloc(num_pert_lab*sizeof(QInt));
        if (new_oper->pert_max_orders==NULL) {
            printf("RSPZeroOperCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for allowed maximal orders");
        }
        new_oper->oper_pert_orders = (QInt *)malloc(num_pert_lab*sizeof(QInt));
        if (new_oper->oper_pert_orders==NULL) {
            printf("RSPZeroOperCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for pert. orders on zero-electron operator");
        }
        new_oper->pert_labels = (QcPertInt *)malloc(num_pert_lab*sizeof(QcPertInt));
        if (new_oper->pert_labels==NULL) {
            printf("RSPZeroOperCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for perturbation labels");
        }
        new_oper->oper_pert_labels = (QcPertInt *)malloc(num_pert_lab*sizeof(QcPertInt));
        if (new_oper->oper_pert_labels==NULL) {
            printf("RSPZeroOperCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for pert. labels on zero-electron operator");
        }
        for (ilab=0; ilab<num_pert_lab; ilab++) {
            if (pert_labels[ilab]>OPENRSP_PERT_LABEL_MAX) {
                printf("RSPZeroOperCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ilab,
                       pert_labels[ilab]);
                printf("RSPZeroOperCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
                       OPENRSP_PERT_LABEL_MAX);
                QErrorExit(FILE_AND_LINE, "invalid perturbation label");
            }
            /* each element of <pert_labels> should be unique */
            for (jlab=0; jlab<ilab; jlab++) {
                if (pert_labels[jlab]==pert_labels[ilab]) {
                    printf("RSPZeroOperCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                           jlab,
                           pert_labels[jlab]);
                    printf("RSPZeroOperCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                           ilab,
                           pert_labels[ilab]);
                    QErrorExit(FILE_AND_LINE, "repeated perturbation labels not allowed");
                }
            }
            new_oper->pert_labels[ilab] = pert_labels[ilab];
            if (pert_max_orders[ilab]<1) {
                printf("RSPZeroOperCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ilab,
                       pert_labels[ilab]);
                printf("RSPZeroOperCreate>> allowed maximal order is %"QINT_FMT"\n",
                       pert_max_orders[ilab]);
                QErrorExit(FILE_AND_LINE, "only positive order allowed");
            }
            new_oper->pert_max_orders[ilab] = pert_max_orders[ilab];
        }
    }
    else {
        new_oper->pert_max_orders = NULL;
        new_oper->oper_pert_orders = NULL;
        new_oper->pert_labels = NULL;
        new_oper->oper_pert_labels = NULL;
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    new_oper->user_ctx = user_ctx;
#endif
    new_oper->get_zero_oper_contrib = get_zero_oper_contrib;
    new_oper->next_oper = NULL;
    *zero_oper = new_oper;
    return QSUCCESS;
}

/* <function name='RSPZeroOperAdd'
             attr='private'
             author='Bin Gao'
             date='2015-02-12'>
     Add a given zero-electron operator to the linked list
     <param name='zero_oper' direction='inout'>
       The linked list of zero-electron operators
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act as
       perturbations on the zero-electron operator
     </param>
     <param name='pert_labels' direction='in'>
       All the different perturbation labels
     </param>
     <param name='pert_max_orders' direction='in'>
       Allowed maximal order of a perturbation described by exactly one of
       the above different labels
     </param>
     <param name='user_ctx' direction='in'>
       User-defined callback-function context
     </param>
     <param name='get_zero_oper_contrib' direction='in'>
       User-specified function for calculating contribution of the
       zero-electron operator and its derivatives
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPZeroOperAdd(RSPZeroOper *zero_oper,
                          const QInt num_pert_lab,
                          const QcPertInt *pert_labels,
                          const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                          void *user_ctx,
#endif
                          const GetZeroOperContrib get_zero_oper_contrib)
{
    RSPZeroOper *new_oper;  /* new operator */
    RSPZeroOper *cur_oper;  /* current operator */
    QErrorCode ierr;        /* error information */
    /* creates the new operator */
    ierr = RSPZeroOperCreate(&new_oper,
                             num_pert_lab,
                             pert_labels,
                             pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             user_ctx,
#endif
                             get_zero_oper_contrib);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPZeroOperCreate()");
    /* walks to the last operator */
    cur_oper = zero_oper;
    while (cur_oper->next_oper!=NULL) {
        cur_oper = cur_oper->next_oper;
    }
    /* inserts the new operator to the tail of the linked list */
    cur_oper->next_oper = new_oper;
    return QSUCCESS;
}

/* <function name='RSPZeroOperAssemble'
             attr='private'
             author='Bin Gao'
             date='2015-02-12'>
     Assembles the context of zero-electron operator
     <param name='zero_oper' direction='inout'>
       The context of zero-electron operator
     </param>
     <param name='rsp_pert' direction='in'>
       The context of perturbations
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPZeroOperAssemble(RSPZeroOper *zero_oper,
                                  const RSPPert *rsp_pert)
{
    QInt ioper;             /* incremental recorder over operators */
    RSPZeroOper *cur_oper;  /* current operator */
    QErrorCode ierr;        /* error information */
    ioper = 0;
    cur_oper = zero_oper;
    do {
        if (cur_oper->num_pert_lab>0 &&
            (cur_oper->pert_labels==NULL || cur_oper->pert_max_orders==NULL)) {
            printf("RSPZeroOperAssemble>> %"QINT_FMT"-th zero-electron operator\n",
                   ioper);
            QErrorExit(FILE_AND_LINE, "perturbations of zero-electron operator not set");
        }
        if (zero_oper->get_zero_oper_contrib==NULL) {
            printf("RSPZeroOperAssemble>> %"QINT_FMT"-th zero-electron operator\n",
                   ioper);
            QErrorExit(FILE_AND_LINE, "callback function of zero-electron operator not set");
        }
        /* checks perturbation labels and allowed maximal orders against
           all known perturbations */
        ierr = RSPPertValidateLabelOrder(rsp_pert,
                                         cur_oper->num_pert_lab,
                                         cur_oper->pert_labels,
                                         cur_oper->pert_max_orders);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertValidateLabelOrder()");
        /* moves to the next operator */
        ioper++;
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}

/* <function name='RSPZeroOperWrite'
             attr='private'
             author='Bin Gao'
             date='2015-02-12'>
     Writes the linked list of zero-electron operators
     <param name='zero_oper' direction='in'>
       The linked list of zero-electron operators
     </param>
     <param name='fp_oper' direction='inout'>File pointer</param>
     <return>Error information</return>
   </function> */
QErrorCode RSPZeroOperWrite(RSPZeroOper *zero_oper, FILE *fp_oper)
{
    QInt ioper;             /* incremental recorder over operators */
    RSPZeroOper *cur_oper;  /* current operator */
    QInt ilab;              /* incremental recorder over perturbation labels */
    ioper = 0;
    cur_oper = zero_oper;
    do {
        fprintf(fp_oper, "RSPZeroOperWrite>> operator %"QINT_FMT"\n", ioper);
        fprintf(fp_oper,
                "RSPZeroOperWrite>> number of pert. labels that zero-electron operator depends on %"QINT_FMT"\n",
                cur_oper->num_pert_lab);
        fprintf(fp_oper, "RSPZeroOperWrite>> label           maximum-order\n");
        for (ilab=0; ilab<cur_oper->num_pert_lab; ilab++) {
            fprintf(fp_oper,
                    "RSPZeroOperWrite>>       %"QCPERTINT_FMT"                  %"QINT_FMT"\n",
                    cur_oper->pert_labels[ilab],
                    cur_oper->pert_max_orders[ilab]);
        }
#if defined(OPENRSP_C_USER_CONTEXT)
        if (cur_oper->user_ctx!=NULL) {
            fprintf(fp_oper, "RSPZeroOperWrite>> user-defined function context given\n");
        }
#endif
        /* moves to the next operator */
        ioper++;
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}

/* <function name='RSPZeroOperGetContribution'
             attr='private'
             author='Bin Gao'
             date='2015-10-15'>
     Calculates contribution of the linked list of zero-electron operators
     <param name='zero_oper' direction='inout'>
       The linked list of zero-electron operators
     </param>
     <param name='oper_len_tuple' direction='in'>
       Length of the perturbation tuple on the linked list of zero-electron
       operators
     </param>
     <param name='oper_pert_tuple' direction='in'>
       Perturbation tuple on the linked list of zero-electron operators
     </param>
     <param name='size_pert' direction='in'>
       Size of the perturbations on the linked list of zero-electron operators
     </param>
     <param name='val_oper' direction='inout'>
       The contribution of the linked list of zero-electron operators
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPZeroOperGetContribution(RSPZeroOper *zero_oper,
                                      const QInt oper_len_tuple,
                                      const QcPertInt *oper_pert_tuple,
                                      const QInt size_pert,
                                      QReal *val_oper)
{
    RSPZeroOper *cur_oper;  /* current operator */
    QErrorCode ierr;        /* error information */
    cur_oper = zero_oper;
    do {
        /* gets perturbation labels and corresponding orders out of the internal
           perturbation tuple on the zero-electron operator */
        ierr = RSPPertInternTupleToHostLabelOrder(oper_len_tuple,
                                                  oper_pert_tuple,
                                                  cur_oper->num_pert_lab,
                                                  cur_oper->pert_labels,
                                                  cur_oper->pert_max_orders,
                                                  &cur_oper->oper_num_pert,
                                                  cur_oper->oper_pert_labels,
                                                  cur_oper->oper_pert_orders);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertInternTupleToHostLabelOrder()");
        /* checks if the perturbations on the zero-electron operator
           result in zero values */
        if (cur_oper->oper_num_pert>=0) {
            /* calculates contribution of zero-electron operator using the
               callback function */
            cur_oper->get_zero_oper_contrib(cur_oper->oper_num_pert,
                                            cur_oper->oper_pert_labels,
                                            cur_oper->oper_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                            cur_oper->user_ctx,
#endif
                                            size_pert,
                                            val_oper);
        }
        /* moves to the next operator */
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}

/* <function name='RSPZeroOperDestroy'
             attr='private'
             author='Bin Gao'
             date='2015-02-12'>
     Destroys the linked list of zero-electron operators, should be called
     at the end
     <param name='zero_oper' direction='inout'>
       The linked list of zero-electron operators
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPZeroOperDestroy(RSPZeroOper **zero_oper)
{
    RSPZeroOper *cur_oper;   /* current operator */
    RSPZeroOper *next_oper;  /* next operator */
    /* walks to the last operator */
    cur_oper = *zero_oper;
    while (cur_oper!=NULL) {
        if (cur_oper->pert_max_orders!=NULL) {
            free(cur_oper->pert_max_orders);
            cur_oper->pert_max_orders = NULL;
        }
        if (cur_oper->oper_pert_orders!=NULL) {
            free(cur_oper->oper_pert_orders);
            cur_oper->oper_pert_orders = NULL;
        }
        if (cur_oper->pert_labels!=NULL) {
            free(cur_oper->pert_labels);
            cur_oper->pert_labels = NULL;
        }
        if (cur_oper->oper_pert_labels!=NULL) {
            free(cur_oper->oper_pert_labels);
            cur_oper->oper_pert_labels = NULL;
        }
#if defined(OPENRSP_C_USER_CONTEXT)
        cur_oper->user_ctx = NULL;
#endif
        cur_oper->get_zero_oper_contrib = NULL;
        next_oper = cur_oper->next_oper;
        free(cur_oper);
        cur_oper = NULL;
        cur_oper = next_oper;
    }
    return QSUCCESS;
}

