/*
  OpenRSP: open-ended library for response theory
  Copyright 2015 Radovan Bast,
                 Daniel H. Friese,
                 Bin Gao,
                 Dan J. Jonsson,
                 Magnus Ringholm,
                 Kenneth Ruud,
                 Andreas Thorvaldsen

  OpenRSP is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.

  OpenRSP is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with OpenRSP. If not, see <http://www.gnu.org/licenses/>.

*/

#include "RSPTwoOper.h"

/* <function name='RSPTwoOperCreate'
             attr='private'
             author='Bin Gao'
             date='2014-08-06'>
     Create a node of a linked list for a given two-electron operator, should
     be called at first
     <param name='two_oper' direction='inout'>
       The linked list of two-electron operators
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act as
       perturbations on the two-electron operator
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
     <param name='get_two_oper_mat' direction='in'>
       User-specified function for calculating integral matrices of the
       two-electron operator and its derivatives
     </param>
     <param name='get_two_oper_exp' direction='in'>
       User-specified function for calculating expectation values of the
       two-electron operator and its derivatives
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPTwoOperCreate(RSPTwoOper **two_oper,
                            const QInt num_pert_lab,
                            const QcPertInt *pert_labels,
                            const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                            void *user_ctx,
#endif
                            const GetTwoOperMat get_two_oper_mat,
                            const GetTwoOperExp get_two_oper_exp)
{
    RSPTwoOper *new_oper;  /* new operator */
    QInt ilab;             /* incremental recorders over perturbation labels */
    QInt jlab;
    new_oper = (RSPTwoOper *)malloc(sizeof(RSPTwoOper));
    if (new_oper==NULL) {
        QErrorExit(FILE_AND_LINE, "allocates memory for two-electron operator");
    }
    if (num_pert_lab<0) {
        printf("RSPTwoOperCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        QErrorExit(FILE_AND_LINE, "invalid number of perturbation labels");
    }
    else if (num_pert_lab>OPENRSP_PERT_LABEL_MAX) {
        printf("RSPTwoOperCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        printf("RSPTwoOperCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
               OPENRSP_PERT_LABEL_MAX);
        QErrorExit(FILE_AND_LINE, "too many perturbation labels");
    }
    new_oper->num_pert_lab = num_pert_lab;
    if (new_oper->num_pert_lab>0) {
        new_oper->pert_max_orders = (QInt *)malloc(num_pert_lab*sizeof(QInt));
        if (new_oper->pert_max_orders==NULL) {
            printf("RSPTwoOperCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for allowed maximal orders");
        }
        new_oper->oper_pert_orders = (QInt *)malloc(num_pert_lab*sizeof(QInt));
        if (new_oper->oper_pert_orders==NULL) {
            printf("RSPTwoOperCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for pert. orders on 2el operator");
        }
        new_oper->pert_labels = (QcPertInt *)malloc(num_pert_lab*sizeof(QcPertInt));
        if (new_oper->pert_labels==NULL) {
            printf("RSPTwoOperCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for perturbation labels");
        }
        new_oper->oper_pert_labels = (QcPertInt *)malloc(num_pert_lab*sizeof(QcPertInt));
        if (new_oper->oper_pert_labels==NULL) {
            printf("RSPTwoOperCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for pert. labels on 2el operator");
        }
        for (ilab=0; ilab<num_pert_lab; ilab++) {
            if (pert_labels[ilab]>OPENRSP_PERT_LABEL_MAX) {
                printf("RSPTwoOperCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ilab,
                       pert_labels[ilab]);
                printf("RSPTwoOperCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
                       OPENRSP_PERT_LABEL_MAX);
                QErrorExit(FILE_AND_LINE, "invalid perturbation label");
            }
            /* each element of <pert_labels> should be unique */
            for (jlab=0; jlab<ilab; jlab++) {
                if (pert_labels[jlab]==pert_labels[ilab]) {
                    printf("RSPTwoOperCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                           jlab,
                           pert_labels[jlab]);
                    printf("RSPTwoOperCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                           ilab,
                           pert_labels[ilab]);
                    QErrorExit(FILE_AND_LINE, "repeated perturbation labels not allowed");
                }
            }
            new_oper->pert_labels[ilab] = pert_labels[ilab];
            if (pert_max_orders[ilab]<1) {
                printf("RSPTwoOperCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ilab,
                       pert_labels[ilab]);
                printf("RSPTwoOperCreate>> allowed maximal order is %"QINT_FMT"\n",
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
    new_oper->get_two_oper_mat = get_two_oper_mat;
    new_oper->get_two_oper_exp = get_two_oper_exp;
    new_oper->next_oper = NULL;
    *two_oper = new_oper;
    return QSUCCESS;
}

/* <function name='RSPTwoOperAdd'
             attr='private'
             author='Bin Gao'
             date='2014-08-06'>
     Add a given two-electron operator to the linked list
     <param name='two_oper' direction='inout'>
       The linked list of two-electron operators
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act as
       perturbations on the two-electron operator
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
     <param name='get_two_oper_mat' direction='in'>
       User-specified function for calculating integral matrices of the
       two-electron operator and its derivatives
     </param>
     <param name='get_two_oper_exp' direction='in'>
       User-specified function for calculating expectation values of the
       two-electron operator and its derivatives
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPTwoOperAdd(RSPTwoOper *two_oper,
                         const QInt num_pert_lab,
                         const QcPertInt *pert_labels,
                         const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                         void *user_ctx,
#endif
                         const GetTwoOperMat get_two_oper_mat,
                         const GetTwoOperExp get_two_oper_exp)
{
    RSPTwoOper *new_oper;  /* new operator */
    RSPTwoOper *cur_oper;  /* current operator */
    QErrorCode ierr;       /* error information */
    /* creates the new operator */
    ierr = RSPTwoOperCreate(&new_oper,
                            num_pert_lab,
                            pert_labels,
                            pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                            user_ctx,
#endif
                            get_two_oper_mat,
                            get_two_oper_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPTwoOperCreate()");
    /* walks to the last operator */
    cur_oper = two_oper;
    while (cur_oper->next_oper!=NULL) {
        cur_oper = cur_oper->next_oper;
    }
    /* inserts the new operator to the tail of the linked list */
    cur_oper->next_oper = new_oper;
    return QSUCCESS;
}

/* <function name='RSPTwoOperAssemble'
             attr='private'
             author='Bin Gao'
             date='2014-08-06'>
     Assembles the linked list of two-electron operators
     <param name='two_oper' direction='inout'>
       The linked list of two-electron operators
     </param>
     <param name='rsp_pert' direction='in'>
       The context of perturbations
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPTwoOperAssemble(RSPTwoOper *two_oper, const RSPPert *rsp_pert)
{
    QInt ioper;            /* incremental recorder over operators */
    RSPTwoOper *cur_oper;  /* current operator */
    QErrorCode ierr;       /* error information */
    ioper = 0;
    cur_oper = two_oper;
    do {
        if (cur_oper->num_pert_lab>0 &&
            (cur_oper->pert_labels==NULL || cur_oper->pert_max_orders==NULL)) {
            printf("RSPTwoOperAssemble>> %"QINT_FMT"-th two-electron operator\n",
                   ioper);
            QErrorExit(FILE_AND_LINE, "perturbations of two-electron operator not set");
        }
        if (cur_oper->get_two_oper_mat==NULL || cur_oper->get_two_oper_exp==NULL) {
            printf("RSPTwoOperAssemble>> %"QINT_FMT"-th two-electron operator\n",
                   ioper);
            QErrorExit(FILE_AND_LINE, "callback functions of two-electron operator not set");
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

/* <function name='RSPTwoOperWrite'
             attr='private'
             author='Bin Gao'
             date='2014-08-06'>
     Writes the linked list of two-electron operators
     <param name='two_oper' direction='in'>
       The linked list of two-electron operators
     </param>
     <param name='fp_oper' direction='inout'>File pointer</param>
     <return>Error information</return>
   </function> */
QErrorCode RSPTwoOperWrite(RSPTwoOper *two_oper, FILE *fp_oper)
{
    QInt ioper;            /* incremental recorder over operators */
    RSPTwoOper *cur_oper;  /* current operator */
    QInt ilab;             /* incremental recorder over perturbation labels */
    ioper = 0;
    cur_oper = two_oper;
    do {
        fprintf(fp_oper, "RSPTwoOperWrite>> operator %"QINT_FMT"\n", ioper);
        fprintf(fp_oper,
                "RSPTwoOperWrite>> number of pert. labels that two-electron operator depends on %"QINT_FMT"\n",
                cur_oper->num_pert_lab);
        fprintf(fp_oper, "RSPTwoOperWrite>> label           maximum-order\n");
        for (ilab=0; ilab<cur_oper->num_pert_lab; ilab++) {
            fprintf(fp_oper,
                    "RSPTwoOperWrite>>       %"QCPERTINT_FMT"                  %"QINT_FMT"\n",
                    cur_oper->pert_labels[ilab],
                    cur_oper->pert_max_orders[ilab]);
        }
#if defined(OPENRSP_C_USER_CONTEXT)
        if (cur_oper->user_ctx!=NULL) {
            fprintf(fp_oper, "RSPTwoOperWrite>> user-defined function context given\n");
        }
#endif
        /* moves to the next operator */
        ioper++;
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}

/* <function name='RSPTwoOperGetMat'
             attr='private'
             author='Bin Gao'
             date='2015-10-15'>
     Calculates integral matrices of the linked list of two-electron operators
     <param name='two_oper' direction='inout'>
       The linked list of two-electron operators
     </param>
     <param name='oper_len_tuple' direction='in'>
       Length of the perturbation tuple on the linked list of two-electron
       operators
     </param>
     <param name='oper_pert_tuple' direction='in'>
       Perturbation tuple on the linked list of two-electron operators
     </param>
     <param name='num_int' direction='in'>
       Number of the integral matrices
     </param>
     <param name='val_int' direction='inout'>
       The integral matrices
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPTwoOperGetMat(RSPTwoOper *two_oper,
                            const QInt oper_len_tuple,
                            const QcPertInt *oper_pert_tuple,
                            const QInt num_dmat,
                            QcMat *dens_mat[],
                            const QInt num_int,
                            QcMat *val_int[])
{
    RSPTwoOper *cur_oper;  /* current operator */
    QErrorCode ierr;       /* error information */
    cur_oper = two_oper;
    do {
        /* gets perturbation labels and corresponding orders out of the internal
           perturbation tuple on the two-electron operator */
        ierr = RSPPertInternTupleToHostLabelOrder(oper_len_tuple,
                                                  oper_pert_tuple,
                                                  cur_oper->num_pert_lab,
                                                  cur_oper->pert_labels,
                                                  cur_oper->pert_max_orders,
                                                  &cur_oper->oper_num_pert,
                                                  cur_oper->oper_pert_labels,
                                                  cur_oper->oper_pert_orders);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertInternTupleToHostLabelOrder()");
        /* checks if the perturbations on the two-electron operator
           result in zero values */
        if (cur_oper->oper_num_pert>=0) {
            /* calculates integral matrices using the callback function */
            cur_oper->get_two_oper_mat(cur_oper->oper_num_pert,
                                       cur_oper->oper_pert_labels,
                                       cur_oper->oper_pert_orders,
                                       num_dmat,
                                       dens_mat,
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

/* <function name='RSPTwoOperGetExp'
             attr='private'
             author='Bin Gao'
             date='2015-10-15'>
     Calculates expectation values of the linked list of two-electron operators
     <param name='two_oper' direction='inout'>
       The linked list of two-electron operators
     </param>
     <param name='oper_len_tuple' direction='in'>
       Length of the perturbation tuple on the linked list of two-electron
       operators
     </param>
     <param name='oper_pert_tuple' direction='in'>
       Perturbation tuple on the linked list of two-electron operators
     </param>
     <param name='dmat_len_tuple' direction='in'>
       Length of different perturbation tuples of the left-hand-side (LHS) and
       right-hand-side (RHS) atomic orbital (AO) based density matrices passed
     </param>
     <param name='num_LHS_dmat' direction='in'>
       Number of LHS AO based density matrices passed for each LHS density
       matrix perturbation tuple
     </param>
     <param name='LHS_dens_mat' direction='in'>
       The LHS AO based density matrices
     </param>
     <param name='num_RHS_dmat' direction='in'>
       Number of RHS AO based density matrices passed for each RHS density
       matrix perturbation tuple
     </param>
     <param name='RHS_dens_mat' direction='in'>
       The RHS AO based density matrices
     </param>
     <param name='num_exp' direction='in'>
       Number of the expectation values
     </param>
     <param name='val_exp' direction='inout'>
       The expectation values
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPTwoOperGetExp(RSPTwoOper *two_oper,
                            const QInt oper_len_tuple,
                            const QcPertInt *oper_pert_tuple,
                            const QInt dmat_len_tuple,
                            const QInt *num_LHS_dmat,
                            QcMat *LHS_dens_mat[],
                            const QInt *num_RHS_dmat,
                            QcMat *RHS_dens_mat[],
                            const QInt num_exp,
                            QReal *val_exp)
{
    RSPTwoOper *cur_oper;  /* current operator */
    QErrorCode ierr;       /* error information */
    cur_oper = two_oper;
    do {
        /* gets perturbation labels and corresponding orders out of the internal
           perturbation tuple on the two-electron operator */
        ierr = RSPPertInternTupleToHostLabelOrder(oper_len_tuple,
                                                  oper_pert_tuple,
                                                  cur_oper->num_pert_lab,
                                                  cur_oper->pert_labels,
                                                  cur_oper->pert_max_orders,
                                                  &cur_oper->oper_num_pert,
                                                  cur_oper->oper_pert_labels,
                                                  cur_oper->oper_pert_orders);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertInternTupleToHostLabelOrder()");
        /* checks if the perturbations on the two-electron operator
           result in zero values */
        if (cur_oper->oper_num_pert>=0) {
            /* calculates expectation values using the callback function */
            cur_oper->get_two_oper_exp(cur_oper->oper_num_pert,
                                       cur_oper->oper_pert_labels,
                                       cur_oper->oper_pert_orders,
                                       dmat_len_tuple,
                                       num_LHS_dmat,
                                       LHS_dens_mat,
                                       num_RHS_dmat,
                                       RHS_dens_mat,
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

/* <function name='RSPTwoOperDestroy'
             attr='private'
             author='Bin Gao'
             date='2014-08-06'>
     Destroys the linked list of two-electron operators, should be called
     at the end
     <param name='two_oper' direction='inout'>
       The linked list of two-electron operators
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPTwoOperDestroy(RSPTwoOper **two_oper)
{
    RSPTwoOper *cur_oper;   /* current operator */
    RSPTwoOper *next_oper;  /* next operator */
    /* walks to the last operator */
    cur_oper = *two_oper;
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
        cur_oper->get_two_oper_mat = NULL;
        cur_oper->get_two_oper_exp = NULL;
        next_oper = cur_oper->next_oper;
        free(cur_oper);
        cur_oper = NULL;
        cur_oper = next_oper;
    }
    return QSUCCESS;
}

