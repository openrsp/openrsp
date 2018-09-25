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

#include "RSPOneOper.h"

/* <function name='RSPOneOperCreate'
             attr='private'
             author='Bin Gao'
             date='2014-07-30'>
     Create a node of a linked list for a given one-electron operator, should
     be called at first
     <param name='one_oper' direction='inout'>
       The linked list of one-electron operators
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act as
       perturbations on the one-electron operator
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
     <param name='get_one_oper_mat' direction='in'>
       User-specified function for calculating integral matrices of the
       one-electron operator and its derivatives
     </param>
     <param name='get_one_oper_exp' direction='in'>
       User-specified function for calculating expectation values of the
       one-electron operator and its derivatives
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPOneOperCreate(RSPOneOper **one_oper,
                            const QInt num_pert_lab,
                            const QcPertInt *pert_labels,
                            const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                            void *user_ctx,
#endif
                            const GetOneOperMat get_one_oper_mat,
                            const GetOneOperExp get_one_oper_exp)
{
    RSPOneOper *new_oper;  /* new operator */
    QInt ilab;             /* incremental recorders over perturbation labels */
    QInt jlab;
    new_oper = (RSPOneOper *)malloc(sizeof(RSPOneOper));
    if (new_oper==NULL) {
        QErrorExit(FILE_AND_LINE, "allocates memory for one-electron operator");
    }
    if (num_pert_lab<0) {
        printf("RSPOneOperCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        QErrorExit(FILE_AND_LINE, "invalid number of perturbation labels");
    }
    else if (num_pert_lab>OPENRSP_PERT_LABEL_MAX) {
        printf("RSPOneOperCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        printf("RSPOneOperCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
               OPENRSP_PERT_LABEL_MAX);
        QErrorExit(FILE_AND_LINE, "too many perturbation labels");
    }
    new_oper->num_pert_lab = num_pert_lab;
    if (new_oper->num_pert_lab>0) {
        new_oper->pert_max_orders = (QInt *)malloc(num_pert_lab*sizeof(QInt));
        if (new_oper->pert_max_orders==NULL) {
            printf("RSPOneOperCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for allowed maximal orders");
        }
        new_oper->oper_pert_orders = (QInt *)malloc(num_pert_lab*sizeof(QInt));
        if (new_oper->oper_pert_orders==NULL) {
            printf("RSPOneOperCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for pert. orders on 1el operator");
        }
        new_oper->pert_labels = (QcPertInt *)malloc(num_pert_lab*sizeof(QcPertInt));
        if (new_oper->pert_labels==NULL) {
            printf("RSPOneOperCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for perturbation labels");
        }
        new_oper->oper_pert_labels = (QcPertInt *)malloc(num_pert_lab*sizeof(QcPertInt));
        if (new_oper->oper_pert_labels==NULL) {
            printf("RSPOneOperCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for pert. labels on 1el operator");
        }
        for (ilab=0; ilab<num_pert_lab; ilab++) {
            if (pert_labels[ilab]>OPENRSP_PERT_LABEL_MAX) {
                printf("RSPOneOperCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ilab,
                       pert_labels[ilab]);
                printf("RSPOneOperCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
                       OPENRSP_PERT_LABEL_MAX);
                QErrorExit(FILE_AND_LINE, "invalid perturbation label");
            }
            /* each element of <pert_labels> should be unique */
            for (jlab=0; jlab<ilab; jlab++) {
                if (pert_labels[jlab]==pert_labels[ilab]) {
                    printf("RSPOneOperCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                           jlab,
                           pert_labels[jlab]);
                    printf("RSPOneOperCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                           ilab,
                           pert_labels[ilab]);
                    QErrorExit(FILE_AND_LINE, "repeated perturbation labels not allowed");
                }
            }
            new_oper->pert_labels[ilab] = pert_labels[ilab];
            if (pert_max_orders[ilab]<1) {
                printf("RSPOneOperCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ilab,
                       pert_labels[ilab]);
                printf("RSPOneOperCreate>> allowed maximal order is %"QINT_FMT"\n",
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
    new_oper->get_one_oper_mat = get_one_oper_mat;
    new_oper->get_one_oper_exp = get_one_oper_exp;
    new_oper->next_oper = NULL;
    *one_oper = new_oper;
    return QSUCCESS;
}

/* <function name='RSPOneOperAdd'
             attr='private'
             author='Bin Gao'
             date='2014-07-30'>
     Add a given one-electron operator to the linked list
     <param name='one_oper' direction='inout'>
       The linked list of one-electron operators
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act as
       perturbations on the one-electron operator
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
     <param name='get_one_oper_mat' direction='in'>
       User-specified function for calculating integral matrices of the
       one-electron operator and its derivatives
     </param>
     <param name='get_one_oper_exp' direction='in'>
       User-specified function for calculating expectation values of the
       one-electron operator and its derivatives
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPOneOperAdd(RSPOneOper *one_oper,
                         const QInt num_pert_lab,
                         const QcPertInt *pert_labels,
                         const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                         void *user_ctx,
#endif
                         const GetOneOperMat get_one_oper_mat,
                         const GetOneOperExp get_one_oper_exp)
{
    RSPOneOper *new_oper;  /* new operator */
    RSPOneOper *cur_oper;  /* current operator */
    QErrorCode ierr;       /* error information */
    /* creates the new operator */
    ierr = RSPOneOperCreate(&new_oper,
                            num_pert_lab,
                            pert_labels,
                            pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                            user_ctx,
#endif
                            get_one_oper_mat,
                            get_one_oper_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOneOperCreate()");
    /* walks to the last operator */
    cur_oper = one_oper;
    while (cur_oper->next_oper!=NULL) {
        cur_oper = cur_oper->next_oper;
    }
    /* inserts the new operator to the tail of the linked list */
    cur_oper->next_oper = new_oper;
    return QSUCCESS;
}

/* <function name='RSPOneOperAssemble'
             attr='private'
             author='Bin Gao'
             date='2014-07-30'>
     Assembles the linked list of one-electron operators
     <param name='one_oper' direction='inout'>
       The linked list of one-electron operators
     </param>
     <param name='rsp_pert' direction='in'>
       The context of perturbations
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPOneOperAssemble(RSPOneOper *one_oper, const RSPPert *rsp_pert)
{
    QInt ioper;            /* incremental recorder over operators */
    RSPOneOper *cur_oper;  /* current operator */
    QErrorCode ierr;       /* error information */
    ioper = 0;
    cur_oper = one_oper;
    do {
        if (cur_oper->num_pert_lab>0 &&
            (cur_oper->pert_labels==NULL || cur_oper->pert_max_orders==NULL)) {
            printf("RSPOneOperAssemble>> %"QINT_FMT"-th one-electron operator\n",
                   ioper);
            QErrorExit(FILE_AND_LINE, "perturbations of one-electron operator not set");
        }
        if (cur_oper->get_one_oper_mat==NULL || cur_oper->get_one_oper_exp==NULL) {
            printf("RSPOneOperAssemble>> %"QINT_FMT"-th one-electron operator\n",
                   ioper);
            QErrorExit(FILE_AND_LINE, "callback functions of one-electron operator not set");
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

/* <function name='RSPOneOperWrite'
             attr='private'
             author='Bin Gao'
             date='2014-07-30'>
     Writes the linked list of one-electron operators
     <param name='one_oper' direction='in'>
       The linked list of one-electron operators
     </param>
     <param name='fp_oper' direction='inout'>File pointer</param>
     <return>Error information</return>
   </function> */
QErrorCode RSPOneOperWrite(RSPOneOper *one_oper, FILE *fp_oper)
{
    QInt ioper;            /* incremental recorder over operators */
    RSPOneOper *cur_oper;  /* current operator */
    QInt ilab;             /* incremental recorder over perturbation labels */
    ioper = 0;
    cur_oper = one_oper;
    do {
        fprintf(fp_oper, "RSPOneOperWrite>> operator %"QINT_FMT"\n", ioper);
        fprintf(fp_oper,
                "RSPOneOperWrite>> number of pert. labels that one-electron operator depends on %"QINT_FMT"\n",
                cur_oper->num_pert_lab);
        fprintf(fp_oper, "RSPOneOperWrite>> label           maximum-order\n");
        for (ilab=0; ilab<cur_oper->num_pert_lab; ilab++) {
            fprintf(fp_oper,
                    "RSPOneOperWrite>>       %"QCPERTINT_FMT"                  %"QINT_FMT"\n",
                    cur_oper->pert_labels[ilab],
                    cur_oper->pert_max_orders[ilab]);
        }
#if defined(OPENRSP_C_USER_CONTEXT)
        if (cur_oper->user_ctx!=NULL) {
            fprintf(fp_oper, "RSPOneOperWrite>> user-defined function context given\n");
        }
#endif
        /* moves to the next operator */
        ioper++;
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}

/* <function name='RSPOneOperGetMat'
             attr='private'
             author='Bin Gao'
             date='2015-10-15'>
     Calculates integral matrices of the linked list of one-electron operators
     <param name='one_oper' direction='inout'>
       The linked list of one-electron operators
     </param>
     <param name='oper_len_tuple' direction='in'>
       Length of the perturbation tuple on the linked list of one-electron
       operators
     </param>
     <param name='oper_pert_tuple' direction='in'>
       Perturbation tuple on the linked list of one-electron operators
     </param>
     <param name='num_int' direction='in'>
       Number of the integral matrices
     </param>
     <param name='val_int' direction='inout'>
       The integral matrices
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPOneOperGetMat(RSPOneOper *one_oper,
                            const QInt oper_len_tuple,
                            const QcPertInt *oper_pert_tuple,
                            const QInt num_int,
                            QcMat *val_int[])
{
    RSPOneOper *cur_oper;  /* current operator */
    QErrorCode ierr;       /* error information */
    cur_oper = one_oper;
    do {
        /* gets perturbation labels and corresponding orders out of the internal
           perturbation tuple on the one-electron operator */
        ierr = RSPPertInternTupleToHostLabelOrder(oper_len_tuple,
                                                  oper_pert_tuple,
                                                  cur_oper->num_pert_lab,
                                                  cur_oper->pert_labels,
                                                  cur_oper->pert_max_orders,
                                                  &cur_oper->oper_num_pert,
                                                  cur_oper->oper_pert_labels,
                                                  cur_oper->oper_pert_orders);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertInternTupleToHostLabelOrder()");
        /* checks if the perturbations on the one-electron operator
           result in zero values */
        if (cur_oper->oper_num_pert>=0) {
            /* calculates integral matrices using the callback function */
            cur_oper->get_one_oper_mat(cur_oper->oper_num_pert,
                                       cur_oper->oper_pert_labels,
                                       cur_oper->oper_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                       cur_oper->user_ctx,
#endif
                                       num_int,
                                       val_int);
        }
        /* moves to the next operator */
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}

/* <function name='RSPOneOperGetExp'
             attr='private'
             author='Bin Gao'
             date='2015-10-15'>
     Calculates expectation values of the linked list of one-electron operators
     <param name='one_oper' direction='inout'>
       The linked list of one-electron operators
     </param>
     <param name='oper_len_tuple' direction='in'>
       Length of the perturbation tuple on the linked list of one-electron
       operators
     </param>
     <param name='oper_pert_tuple' direction='in'>
       Perturbation tuple on the linked list of one-electron operators
     </param>
     <param name='num_dmat' direction='in'>
       Number of atomic orbital (AO) based density matrices
     </param>
     <param name='dens_mat' direction='in'>
       The AO based density matrices
     </param>
     <param name='num_exp' direction='in'>
       Number of the expectation values
     </param>
     <param name='val_exp' direction='inout'>
       The expectation values
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPOneOperGetExp(RSPOneOper *one_oper,
                            const QInt oper_len_tuple,
                            const QcPertInt *oper_pert_tuple,
                            const QInt num_dmat,
                            QcMat *dens_mat[],
                            const QInt num_exp,
                            QReal *val_exp)
{
    RSPOneOper *cur_oper;  /* current operator */
    QErrorCode ierr;       /* error information */
    cur_oper = one_oper;
    do {
        /* gets perturbation labels and corresponding orders out of the internal
           perturbation tuple on the one-electron operator */
        ierr = RSPPertInternTupleToHostLabelOrder(oper_len_tuple,
                                                  oper_pert_tuple,
                                                  cur_oper->num_pert_lab,
                                                  cur_oper->pert_labels,
                                                  cur_oper->pert_max_orders,
                                                  &cur_oper->oper_num_pert,
                                                  cur_oper->oper_pert_labels,
                                                  cur_oper->oper_pert_orders);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertInternTupleToHostLabelOrder()");
        /* checks if the perturbations on the one-electron operator
           result in zero values */
        if (cur_oper->oper_num_pert>=0) {
            /* calculates expectation values using the callback function */
            cur_oper->get_one_oper_exp(cur_oper->oper_num_pert,
                                       cur_oper->oper_pert_labels,
                                       cur_oper->oper_pert_orders,
                                       num_dmat,
                                       dens_mat,
#if defined(OPENRSP_C_USER_CONTEXT)
                                       cur_oper->user_ctx,
#endif
                                       num_exp,
                                       val_exp);
        }
        /* moves to the next operator */
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}

/* <function name='RSPOneOperDestroy'
             attr='private'
             author='Bin Gao'
             date='2014-07-30'>
     Destroys the linked list of one-electron operators, should be called
     at the end
     <param name='one_oper' direction='inout'>
       The linked list of one-electron operators
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPOneOperDestroy(RSPOneOper **one_oper)
{
    RSPOneOper *cur_oper;   /* current operator */
    RSPOneOper *next_oper;  /* next operator */
    /* walks to the last operator */
    cur_oper = *one_oper;
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
        cur_oper->get_one_oper_mat = NULL;
        cur_oper->get_one_oper_exp = NULL;
        next_oper = cur_oper->next_oper;
        free(cur_oper);
        cur_oper = NULL;
        cur_oper = next_oper;
    }
    return QSUCCESS;
}

