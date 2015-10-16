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

#include "RSPXCFun.h"

/* <function name='RSPXCFunCreate'
             attr='private'
             author='Bin Gao'
             date='2015-06-23'>
     Create a node of a linked list for a given XC functional, should
     be called at first
     <param name='xc_fun' direction='inout'>
       The linked list of XC functionals
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act as
       perturbations on the XC functional
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
     <param name='get_xc_fun_mat' direction='in'>
       User-specified function for calculating integral matrices of the
       XC functional and its derivatives
     </param>
     <param name='get_xc_fun_exp' direction='in'>
       User-specified function for calculating expectation values of the
       XC functional and its derivatives
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPXCFunCreate(RSPXCFun **xc_fun,
                          const QInt num_pert_lab,
                          const QcPertInt *pert_labels,
                          const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                          void *user_ctx,
#endif
                          const GetXCFunMat get_xc_fun_mat,
                          const GetXCFunExp get_xc_fun_exp)
{
    RSPXCFun *new_xc;  /* new XC functional */
    QInt ilab;         /* incremental recorders over perturbation labels */
    QInt jlab;
    new_xc = (RSPXCFun *)malloc(sizeof(RSPXCFun));
    if (new_xc==NULL) {
        QErrorExit(FILE_AND_LINE, "allocates memory for XC functional");
    }
    if (num_pert_lab<0) {
        printf("RSPXCFunCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        QErrorExit(FILE_AND_LINE, "invalid number of perturbation labels");
    }
    else if (num_pert_lab>OPENRSP_PERT_LABEL_MAX) {
        printf("RSPXCFunCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        printf("RSPXCFunCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
               OPENRSP_PERT_LABEL_MAX);
        QErrorExit(FILE_AND_LINE, "too many perturbation labels");
    }
    new_xc->num_pert_lab = num_pert_lab;
    if (new_xc->num_pert_lab>0) {
        new_xc->pert_max_orders = (QInt *)malloc(num_pert_lab*sizeof(QInt));
        if (new_xc->pert_max_orders==NULL) {
            printf("RSPXCFunCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for allowed maximal orders");
        }
        new_xc->pert_labels = (QcPertInt *)malloc(num_pert_lab*sizeof(QcPertInt));
        if (new_xc->pert_labels==NULL) {
            printf("RSPXCFunCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for perturbation labels");
        }
        new_xc->xc_len_tuple = 0;
        for (ilab=0; ilab<num_pert_lab; ilab++) {
            if (pert_labels[ilab]>OPENRSP_PERT_LABEL_MAX) {
                printf("RSPXCFunCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ilab,
                       pert_labels[ilab]);
                printf("RSPXCFunCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
                       OPENRSP_PERT_LABEL_MAX);
                QErrorExit(FILE_AND_LINE, "invalid perturbation label");
            }
            /* each element of <pert_labels> should be unique */
            for (jlab=0; jlab<ilab; jlab++) {
                if (pert_labels[jlab]==pert_labels[ilab]) {
                    printf("RSPXCFunCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                           jlab,
                           pert_labels[jlab]);
                    printf("RSPXCFunCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                           ilab,
                           pert_labels[ilab]);
                    QErrorExit(FILE_AND_LINE, "repeated perturbation labels not allowed");
                }
            }
            new_xc->pert_labels[ilab] = pert_labels[ilab];
            if (pert_max_orders[ilab]<1) {
                printf("RSPXCFunCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ilab,
                       pert_labels[ilab]);
                printf("RSPXCFunCreate>> allowed maximal order is %"QINT_FMT"\n",
                       pert_max_orders[ilab]);
                QErrorExit(FILE_AND_LINE, "only positive order allowed");
            }
            new_xc->pert_max_orders[ilab] = pert_max_orders[ilab];
            new_xc->xc_len_tuple += pert_max_orders[ilab];
        }
        new_xc->xc_pert_tuple = (QcPertInt *)malloc(new_xc->xc_len_tuple*sizeof(QcPertInt));
        if (new_xc->xc_pert_tuple==NULL) {
            printf("RSPXCFunCreate>> length of perturbation tuple %"QINT_FMT"\n",
                   new_xc->xc_len_tuple);
            QErrorExit(FILE_AND_LINE, "allocates memory for pert. tuple on XC functional");
        }
    }
    else {
        new_xc->pert_max_orders = NULL;
        new_xc->pert_labels = NULL;
        new_xc->xc_pert_tuple = NULL;
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    new_xc->user_ctx = user_ctx;
#endif
    new_xc->get_xc_fun_mat = get_xc_fun_mat;
    new_xc->get_xc_fun_exp = get_xc_fun_exp;
    new_xc->next_xc = NULL;
    *xc_fun = new_xc;
    return QSUCCESS;
}

/* <function name='RSPXCFunAdd'
             attr='private'
             author='Bin Gao'
             date='2015-06-23'>
     Add a given XC functional to the linked list
     <param name='xc_fun' direction='inout'>
       The linked list of XC functionals
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act as
       perturbations on the XC functional
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
     <param name='get_xc_fun_mat' direction='in'>
       User-specified function for calculating integral matrices of the
       XC functional and its derivatives
     </param>
     <param name='get_xc_fun_exp' direction='in'>
       User-specified function for calculating expectation values of the
       XC functional and its derivatives
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPXCFunAdd(RSPXCFun *xc_fun,
                       const QInt num_pert_lab,
                       const QcPertInt *pert_labels,
                       const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                       void *user_ctx,
#endif
                       const GetXCFunMat get_xc_fun_mat,
                       const GetXCFunExp get_xc_fun_exp)
{
    RSPXCFun *new_xc;  /* new XC functional */
    RSPXCFun *cur_xc;  /* current XC functional */
    QErrorCode ierr;   /* error information */
    /* creates the new XC functional */
    ierr = RSPXCFunCreate(&new_xc,
                          num_pert_lab,
                          pert_labels,
                          pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                          user_ctx,
#endif
                          get_xc_fun_mat,
                          get_xc_fun_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPXCFunCreate()");
    /* walks to the last XC functional */
    cur_xc = xc_fun;
    while (cur_xc->next_xc!=NULL) {
        cur_xc = cur_xc->next_xc;
    }
    /* inserts the new XC functional to the tail of the linked list */
    cur_xc->next_xc = new_xc;
    return QSUCCESS;
}

/* <function name='RSPXCFunAssemble'
             attr='private'
             author='Bin Gao'
             date='2015-06-23'>
     Assembles the linked list of XC functionals
     <param name='xc_fun' direction='inout'>
       The linked list of XC functionals
     </param>
     <param name='rsp_pert' direction='in'>
       The context of perturbations
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPXCFunAssemble(RSPXCFun *xc_fun, const RSPPert *rsp_pert)
{
    QInt ixc;          /* incremental recorder over XC functionals */
    RSPXCFun *cur_xc;  /* current XC functional */
    QErrorCode ierr;   /* error information */
    ixc = 0;
    cur_xc = xc_fun;
    do {
        if (cur_xc->num_pert_lab>0 &&
            (cur_xc->pert_labels==NULL || cur_xc->pert_max_orders==NULL)) {
            printf("RSPXCFunAssemble>> %"QINT_FMT"-th XC functional\n",
                   ixc);
            QErrorExit(FILE_AND_LINE, "perturbations of XC functional not set");
        }
        if (cur_xc->get_xc_fun_mat==NULL || cur_xc->get_xc_fun_exp==NULL) {
            printf("RSPXCFunAssemble>> %"QINT_FMT"-th XC functional\n",
                   ixc);
            QErrorExit(FILE_AND_LINE, "callback functions of XC functional not set");
        }
        /* checks perturbation labels and allowed maximal orders against
           all known perturbations */
        ierr = RSPPertValidateLabelOrder(rsp_pert,
                                         cur_xc->num_pert_lab,
                                         cur_xc->pert_labels,
                                         cur_xc->pert_max_orders);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertValidateLabelOrder()");
        /* moves to the next XC functional */
        ixc++;
        cur_xc = cur_xc->next_xc;
    } while (cur_xc!=NULL);
    return QSUCCESS;
}

/* <function name='RSPXCFunWrite'
             attr='private'
             author='Bin Gao'
             date='2015-06-23'>
     Writes the linked list of XC functionals
     <param name='xc_fun' direction='in'>
       The linked list of XC functionals
     </param>
     <param name='fp_xc' direction='inout'>File pointer</param>
     <return>Error information</return>
   </function> */
QErrorCode RSPXCFunWrite(RSPXCFun *xc_fun, FILE *fp_xc)
{
    QInt ixc;          /* incremental recorder over XC functionals */
    RSPXCFun *cur_xc;  /* current XC functional */
    QInt ilab;         /* incremental recorder over perturbation labels */
    ixc = 0;
    cur_xc = xc_fun;
    do {
        fprintf(fp_xc, "RSPXCFunWrite>> XC functional %"QINT_FMT"\n", ixc);
        fprintf(fp_xc,
                "RSPXCFunWrite>> number of pert. labels that XC functional depends on %"QINT_FMT"\n",
                cur_xc->num_pert_lab);
        fprintf(fp_xc, "RSPXCFunWrite>> label           maximum-order\n");
        for (ilab=0; ilab<cur_xc->num_pert_lab; ilab++) {
            fprintf(fp_xc,
                    "RSPXCFunWrite>>       %"QCPERTINT_FMT"                  %"QINT_FMT"\n",
                    cur_xc->pert_labels[ilab],
                    cur_xc->pert_max_orders[ilab]);
        }
#if defined(OPENRSP_C_USER_CONTEXT)
        if (cur_xc->user_ctx!=NULL) {
            fprintf(fp_xc, "RSPXCFunWrite>> user-defined function context given\n");
        }
#endif
        /* moves to the next XC functional */
        ixc++;
        cur_xc = cur_xc->next_xc;
    } while (cur_xc!=NULL);
    return QSUCCESS;
}

/* <function name='RSPXCFunGetMat'
             attr='private'
             author='Bin Gao'
             date='2015-10-15'>
     Calculates integral matrices of the linked list of XC functionals
     <param name='xc_fun' direction='inout'>
       The linked list of XC functionals
     </param>
     <param name='xc_len_tuple' direction='in'>
       Length of the perturbation tuple on the linked list of XC functionals
     </param>
     <param name='xc_pert_tuple' direction='in'>
       Perturbation tuple on the linked list of XC functionals
     </param>
     <param name='num_freq_configs' direction='in'>
       The number of different frequency configurations to be considered for
       the perturbation tuple
     </param>
     <param name='dmat_num_tuple' direction='in'>
       The number of different perturbation tuples of the atomic orbital (AO)
       based density matrices passed
     </param>
     <param name='dmat_idx_tuple' direction='in'>
       Indices of the density matrix perturbation tuples passed (canonically
       ordered)
     </param>
     <param name='num_dmat' direction='in'>
       Number of collected AO based density matrices for the passed density
       matrix perturbation tuples and all frequency configurations
     </param>
     <param name='dens_mat' direction='in'>
       The collected AO based density matrices
     </param>
     <param name='num_int' direction='in'>
       Number of the integral matrices
     </param>
     <param name='val_int' direction='inout'>
       The integral matrices
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPXCFunGetMat(RSPXCFun *xc_fun,
                          const QInt xc_len_tuple,
                          const QcPertInt *xc_pert_tuple,
                          const QInt num_freq_configs,
                          const QInt dmat_num_tuple,
                          const QInt *dmat_idx_tuple,
                          const QInt num_dmat,
                          QcMat *dens_mat[],
                          const QInt num_int,
                          QcMat *val_int[])
{
    RSPXCFun *cur_xc;  /* current XC functional */
    QErrorCode ierr;   /* error information */
    cur_xc = xc_fun;
    do {
        /* gets the host program's perturbation tuple on the XC functional */
        ierr = RSPPertInternTupleToHostTuple(xc_len_tuple,
                                             xc_pert_tuple,
                                             cur_xc->num_pert_lab,
                                             cur_xc->pert_labels,
                                             cur_xc->pert_max_orders,
                                             &cur_xc->xc_len_tuple,
                                             cur_xc->xc_pert_tuple);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertInternTupleToHostTuple()");
        /* checks if the perturbations on the XC functional result in
           zero values */
        if (cur_xc->xc_len_tuple<0) continue;
        /* calculates integral matrices using the callback function */
        cur_xc->get_xc_fun_mat(cur_xc->xc_len_tuple,
                               cur_xc->xc_pert_tuple,
                               num_freq_configs,
                               dmat_num_tuple,
                               dmat_idx_tuple,
                               num_dmat,
                               dens_mat,
#if defined(OPENRSP_C_USER_CONTEXT)
                               cur_xc->user_ctx,
#endif
                               num_int,
                               val_int);
        /* moves to the next XC functional */
        cur_xc = cur_xc->next_xc;
    } while (cur_xc!=NULL);
    return QSUCCESS;
}

/* <function name='RSPXCFunGetExp'
             attr='private'
             author='Bin Gao'
             date='2015-10-15'>
     Calculates expectation values of the linked list of XC functionals
     <param name='xc_fun' direction='inout'>
       The linked list of XC functionals
     </param>
     <param name='xc_len_tuple' direction='in'>
       Length of the perturbation tuple on the linked list of XC functionals
     </param>
     <param name='xc_pert_tuple' direction='in'>
       Perturbation tuple on the linked list of XC functionals
     </param>
     <param name='num_freq_configs' direction='in'>
       The number of different frequency configurations to be considered for
       the perturbation tuple
     </param>
     <param name='dmat_num_tuple' direction='in'>
       The number of different perturbation tuples of the atomic orbital (AO)
       based density matrices passed
     </param>
     <param name='dmat_idx_tuple' direction='in'>
       Indices of the density matrix perturbation tuples passed (canonically
       ordered)
     </param>
     <param name='num_dmat' direction='in'>
       Number of collected AO based density matrices for the passed density
       matrix perturbation tuples and all frequency configurations
     </param>
     <param name='dens_mat' direction='in'>
       The collected AO based density matrices
     </param>
     <param name='num_exp' direction='in'>
       Number of the expectation values
     </param>
     <param name='val_exp' direction='inout'>
       The expectation values
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPXCFunGetExp(RSPXCFun *xc_fun,
                          const QInt xc_len_tuple,
                          const QcPertInt *xc_pert_tuple,
                          const QInt num_freq_configs,
                          const QInt dmat_num_tuple,
                          const QInt *dmat_idx_tuple,
                          const QInt num_dmat,
                          QcMat *dens_mat[],
                          const QInt num_exp,
                          QReal *val_exp)
{
    RSPXCFun *cur_xc;  /* current XC functional */
    QErrorCode ierr;   /* error information */
    cur_xc = xc_fun;
    do {
        /* gets the host program's perturbation tuple on the XC functional */
        ierr = RSPPertInternTupleToHostTuple(xc_len_tuple,
                                             xc_pert_tuple,
                                             cur_xc->num_pert_lab,
                                             cur_xc->pert_labels,
                                             cur_xc->pert_max_orders,
                                             &cur_xc->xc_len_tuple,
                                             cur_xc->xc_pert_tuple);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertInternTupleToHostTuple()");
        /* checks if the perturbations on the XC functional result in
           zero values */
        if (cur_xc->xc_len_tuple<0) continue;
        /* calculates expectation values using the callback function */
        cur_xc->get_xc_fun_exp(cur_xc->xc_len_tuple,
                               cur_xc->xc_pert_tuple,
                               num_freq_configs,
                               dmat_num_tuple,
                               dmat_idx_tuple,
                               num_dmat,
                               dens_mat,
#if defined(OPENRSP_C_USER_CONTEXT)
                               cur_xc->user_ctx,
#endif
                               num_exp,
                               val_exp);
        /* moves to the next XC functional */
        cur_xc = cur_xc->next_xc;
    } while (cur_xc!=NULL);
    return QSUCCESS;
}

/* <function name='RSPXCFunDestroy'
             attr='private'
             author='Bin Gao'
             date='2015-06-23'>
     Destroys the linked list of XC functionals, should be called at the end
     <param name='xc_fun' direction='inout'>
       The linked list of XC functionals
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPXCFunDestroy(RSPXCFun **xc_fun)
{
    RSPXCFun *cur_xc;   /* current XC functional */
    RSPXCFun *next_xc;  /* next XC functional */
    /* walks to the last XC functional */
    cur_xc = *xc_fun;
    while (cur_xc!=NULL) {
        if (cur_xc->pert_max_orders!=NULL) {
            free(cur_xc->pert_max_orders);
            cur_xc->pert_max_orders = NULL;
        }
        if (cur_xc->pert_labels!=NULL) {
            free(cur_xc->pert_labels);
            cur_xc->pert_labels = NULL;
        }
        if (cur_xc->xc_pert_tuple!=NULL) {
            free(cur_xc->xc_pert_tuple);
            cur_xc->xc_pert_tuple = NULL;
        }
#if defined(OPENRSP_C_USER_CONTEXT)
        cur_xc->user_ctx = NULL;
#endif
        cur_xc->get_xc_fun_mat = NULL;
        cur_xc->get_xc_fun_exp = NULL;
        next_xc = cur_xc->next_xc;
        free(cur_xc);
        cur_xc = NULL;
        cur_xc = next_xc;
    }
    return QSUCCESS;
}

