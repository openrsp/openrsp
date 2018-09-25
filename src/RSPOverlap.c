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

#include "RSPOverlap.h"

/* <function name='RSPOverlapCreate'
             attr='private'
             author='Bin Gao'
             date='2014-08-05'>
     Create the context of overlap operator, should be called at first
     <param name='overlap' direction='inout'>
       The context of overlap operator
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels that can act as
       perturbations on the overlap operator
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
     <param name='get_overlap_mat' direction='in'>
       User-specified function for calculating integral matrices of the
       overlap operator and its derivatives
     </param>
     <param name='get_overlap_exp' direction='in'>
       User-specified function for calculating expectation values of the
       overlap operator and its derivatives
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPOverlapCreate(RSPOverlap *overlap,
                            const QInt num_pert_lab,
                            const QcPertInt *pert_labels,
                            const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                            void *user_ctx,
#endif
                            const GetOverlapMat get_overlap_mat,
                            const GetOverlapExp get_overlap_exp)
{
    QInt ilab;  /* incremental recorders over perturbation labels */
    QInt jlab;
    if (num_pert_lab<0) {
        printf("RSPOverlapCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        QErrorExit(FILE_AND_LINE, "invalid number of perturbation labels");
    }
    else if (num_pert_lab>OPENRSP_PERT_LABEL_MAX) {
        printf("RSPOverlapCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        printf("RSPOverlapCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
               OPENRSP_PERT_LABEL_MAX);
        QErrorExit(FILE_AND_LINE, "too many perturbation labels");
    }
    overlap->num_pert_lab = num_pert_lab;
    if (overlap->num_pert_lab>0) {
        overlap->pert_max_orders = (QInt *)malloc(num_pert_lab*sizeof(QInt));
        if (overlap->pert_max_orders==NULL) {
            printf("RSPOverlapCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for allowed maximal orders");
        }
        overlap->bra_pert_orders = (QInt *)malloc(num_pert_lab*sizeof(QInt));
        if (overlap->bra_pert_orders==NULL) {
            printf("RSPOverlapCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for pert. orders on bra center");
        }
        overlap->ket_pert_orders = (QInt *)malloc(num_pert_lab*sizeof(QInt));
        if (overlap->ket_pert_orders==NULL) {
            printf("RSPOverlapCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for pert. orders on ket center");
        }
        overlap->oper_pert_orders = (QInt *)malloc(num_pert_lab*sizeof(QInt));
        if (overlap->oper_pert_orders==NULL) {
            printf("RSPOverlapCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for pert. orders on overlap operator");
        }
        overlap->pert_labels = (QcPertInt *)malloc(num_pert_lab*sizeof(QcPertInt));
        if (overlap->pert_labels==NULL) {
            printf("RSPOverlapCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for perturbation labels");
        }
        overlap->bra_pert_labels = (QcPertInt *)malloc(num_pert_lab*sizeof(QcPertInt));
        if (overlap->bra_pert_labels==NULL) {
            printf("RSPOverlapCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for pert. labels on bra center");
        }
        overlap->ket_pert_labels = (QcPertInt *)malloc(num_pert_lab*sizeof(QcPertInt));
        if (overlap->ket_pert_labels==NULL) {
            printf("RSPOverlapCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for pert. labels on ket center");
        }
        overlap->oper_pert_labels = (QcPertInt *)malloc(num_pert_lab*sizeof(QcPertInt));
        if (overlap->oper_pert_labels==NULL) {
            printf("RSPOverlapCreate>> number of perturbation labels %"QINT_FMT"\n",
                   num_pert_lab);
            QErrorExit(FILE_AND_LINE, "allocates memory for pert. labels on overlap operator");
        }
        for (ilab=0; ilab<num_pert_lab; ilab++) {
            if (pert_labels[ilab]>OPENRSP_PERT_LABEL_MAX) {
                printf("RSPOverlapCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ilab,
                       pert_labels[ilab]);
                printf("RSPOverlapCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
                       OPENRSP_PERT_LABEL_MAX);
                QErrorExit(FILE_AND_LINE, "invalid perturbation label");
            }
            /* each element of <pert_labels> should be unique */
            for (jlab=0; jlab<ilab; jlab++) {
                if (pert_labels[jlab]==pert_labels[ilab]) {
                    printf("RSPOverlapCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                           jlab,
                           pert_labels[jlab]);
                    printf("RSPOverlapCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                           ilab,
                           pert_labels[ilab]);
                    QErrorExit(FILE_AND_LINE, "repeated perturbation labels not allowed");
                }
            }
            overlap->pert_labels[ilab] = pert_labels[ilab];
            if (pert_max_orders[ilab]<1) {
                printf("RSPOverlapCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ilab,
                       pert_labels[ilab]);
                printf("RSPOverlapCreate>> allowed maximal order is %"QINT_FMT"\n",
                       pert_max_orders[ilab]);
                QErrorExit(FILE_AND_LINE, "only positive order allowed");
            }
            overlap->pert_max_orders[ilab] = pert_max_orders[ilab];
        }
    }
    else {
        overlap->pert_max_orders = NULL;
        overlap->bra_pert_orders = NULL;
        overlap->ket_pert_orders = NULL;
        overlap->oper_pert_orders = NULL;
        overlap->pert_labels = NULL;
        overlap->bra_pert_labels = NULL;
        overlap->ket_pert_labels = NULL;
        overlap->oper_pert_labels = NULL;
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    overlap->user_ctx = user_ctx;
#endif
    overlap->get_overlap_mat = get_overlap_mat;
    overlap->get_overlap_exp = get_overlap_exp;
    return QSUCCESS;
}

/* <function name='RSPOverlapAssemble'
             attr='private'
             author='Bin Gao'
             date='2014-08-05'>
     Assembles the context of overlap operator
     <param name='overlap' direction='inout'>
       The context of overlap operator
     </param>
     <param name='rsp_pert' direction='in'>
       The context of perturbations
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPOverlapAssemble(RSPOverlap *overlap, const RSPPert *rsp_pert)
{
    QErrorCode ierr;  /* error information */
    if (overlap->num_pert_lab>0 &&
        (overlap->pert_labels==NULL || overlap->pert_max_orders==NULL)) {
        QErrorExit(FILE_AND_LINE, "perturbations of overlap operator not set");
    }
    if (overlap->get_overlap_mat==NULL || overlap->get_overlap_exp==NULL) {
        QErrorExit(FILE_AND_LINE, "callback functions of overlap operator not set");
    }
    /* checks perturbation labels and allowed maximal orders against
       all known perturbations */
    ierr = RSPPertValidateLabelOrder(rsp_pert,
                                     overlap->num_pert_lab,
                                     overlap->pert_labels,
                                     overlap->pert_max_orders);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertValidateLabelOrder()");
    return QSUCCESS;
}

/* <function name='RSPOverlapWrite'
             attr='private'
             author='Bin Gao'
             date='2014-08-05'>
     Writes the context of overlap operator
     <param name='overlap' direction='in'>
       The context of overlap operator
     </param>
     <param name='fp_overlap' direction='inout'>File pointer</param>
     <return>Error information</return>
   </function> */
QErrorCode RSPOverlapWrite(const RSPOverlap *overlap, FILE *fp_overlap)
{
    QInt ilab;  /* incremental recorder over perturbation labels */
    fprintf(fp_overlap,
            "RSPOverlapWrite>> number of pert. labels that overlap operator depends on %"QINT_FMT"\n",
            overlap->num_pert_lab);
    fprintf(fp_overlap, "RSPOverlapWrite>> label           maximum-order\n");
    for (ilab=0; ilab<overlap->num_pert_lab; ilab++) {
        fprintf(fp_overlap,
                "RSPOverlapWrite>>       %"QCPERTINT_FMT"                  %"QINT_FMT"\n",
                overlap->pert_labels[ilab],
                overlap->pert_max_orders[ilab]);
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    if (overlap->user_ctx!=NULL) {
        fprintf(fp_overlap, "RSPOverlapWrite>> user-defined function context given\n");
    }
#endif
    return QSUCCESS;
}

/* <function name='RSPOverlapGetMat'
             attr='private'
             author='Bin Gao'
             date='2015-10-15'>
     Calculates integral matrices of the overlap operator
     <param name='overlap' direction='inout'>
       The context of overlap operator
     </param>
     <param name='bra_len_tuple' direction='in'>
       Length of the perturbation tuple on the bra center
     </param>
     <param name='bra_pert_tuple' direction='in'>
       Perturbation tuple on the bra center
     </param>
     <param name='ket_len_tuple' direction='in'>
       Length of the perturbation tuple on the ket center
     </param>
     <param name='ket_pert_tuple' direction='in'>
       Perturbation tuple on the ket center
     </param>
     <param name='oper_len_tuple' direction='in'>
       Length of the perturbation tuple on the overlap operator
     </param>
     <param name='oper_pert_tuple' direction='in'>
       Perturbation tuple on the overlap operator
     </param>
     <param name='num_int' direction='in'>
       Number of the integral matrices
     </param>
     <param name='val_int' direction='inout'>
       The integral matrices
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPOverlapGetMat(RSPOverlap *overlap,
                            const QInt bra_len_tuple,
                            const QcPertInt *bra_pert_tuple,
                            const QInt ket_len_tuple,
                            const QcPertInt *ket_pert_tuple,
                            const QInt oper_len_tuple,
                            const QcPertInt *oper_pert_tuple,
                            const QInt num_int,
                            QcMat *val_int[])
{
    QErrorCode ierr;  /* error information */
    /* gets perturbation labels and corresponding orders out of the internal
       perturbation tuple on the bra center */
    ierr = RSPPertInternTupleToHostLabelOrder(bra_len_tuple,
                                              bra_pert_tuple,
                                              overlap->num_pert_lab,
                                              overlap->pert_labels,
                                              overlap->pert_max_orders,
                                              &overlap->bra_num_pert,
                                              overlap->bra_pert_labels,
                                              overlap->bra_pert_orders);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertInternTupleToHostLabelOrder()");
    /* checks if the perturbations on the bra center result in zero values */
    if (overlap->bra_num_pert<0) return QSUCCESS;
    /* performs the same procedure for the perturbations on the ket center */
    ierr = RSPPertInternTupleToHostLabelOrder(ket_len_tuple,
                                              ket_pert_tuple,
                                              overlap->num_pert_lab,
                                              overlap->pert_labels,
                                              overlap->pert_max_orders,
                                              &overlap->ket_num_pert,
                                              overlap->ket_pert_labels,
                                              overlap->ket_pert_orders);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertInternTupleToHostLabelOrder()");
    if (overlap->ket_num_pert<0) return QSUCCESS;
    /* performs the same procedure for the perturbations on the overlap operator */
    ierr = RSPPertInternTupleToHostLabelOrder(oper_len_tuple,
                                              oper_pert_tuple,
                                              overlap->num_pert_lab,
                                              overlap->pert_labels,
                                              overlap->pert_max_orders,
                                              &overlap->oper_num_pert,
                                              overlap->oper_pert_labels,
                                              overlap->oper_pert_orders);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertInternTupleToHostLabelOrder()");
    if (overlap->oper_num_pert<0) return QSUCCESS;
    /* calculates integral matrices using the callback function */
    overlap->get_overlap_mat(overlap->bra_num_pert,
                             overlap->bra_pert_labels,
                             overlap->bra_pert_orders,
                             overlap->ket_num_pert,
                             overlap->ket_pert_labels,
                             overlap->ket_pert_orders,
                             overlap->oper_num_pert,
                             overlap->oper_pert_labels,
                             overlap->oper_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             overlap->user_ctx,
#endif
                             num_int,
                             val_int);
    return QSUCCESS;
}

/* <function name='RSPOverlapGetExp'
             attr='private'
             author='Bin Gao'
             date='2015-10-15'>
     Calculates expectation values of the overlap operator
     <param name='overlap' direction='inout'>
       The context of overlap operator
     </param>
     <param name='bra_len_tuple' direction='in'>
       Length of the perturbation tuple on the bra center
     </param>
     <param name='bra_pert_tuple' direction='in'>
       Perturbation tuple on the bra center
     </param>
     <param name='ket_len_tuple' direction='in'>
       Length of the perturbation tuple on the ket center
     </param>
     <param name='ket_pert_tuple' direction='in'>
       Perturbation tuple on the ket center
     </param>
     <param name='oper_len_tuple' direction='in'>
       Length of the perturbation tuple on the overlap operator
     </param>
     <param name='oper_pert_tuple' direction='in'>
       Perturbation tuple on the overlap operator
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
QErrorCode RSPOverlapGetExp(RSPOverlap *overlap,
                            const QInt bra_len_tuple,
                            const QcPertInt *bra_pert_tuple,
                            const QInt ket_len_tuple,
                            const QcPertInt *ket_pert_tuple,
                            const QInt oper_len_tuple,
                            const QcPertInt *oper_pert_tuple,
                            const QInt num_dmat,
                            QcMat *dens_mat[],
                            const QInt num_exp,
                            QReal *val_exp)
{
    QErrorCode ierr;  /* error information */
    /* gets perturbation labels and corresponding orders out of the internal
       perturbation tuple on the bra center */
    ierr = RSPPertInternTupleToHostLabelOrder(bra_len_tuple,
                                              bra_pert_tuple,
                                              overlap->num_pert_lab,
                                              overlap->pert_labels,
                                              overlap->pert_max_orders,
                                              &overlap->bra_num_pert,
                                              overlap->bra_pert_labels,
                                              overlap->bra_pert_orders);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertInternTupleToHostLabelOrder()");
    /* checks if the perturbations on the bra center result in zero values */
    if (overlap->bra_num_pert<0) return QSUCCESS;
    /* performs the same procedure for the perturbations on the ket center */
    ierr = RSPPertInternTupleToHostLabelOrder(ket_len_tuple,
                                              ket_pert_tuple,
                                              overlap->num_pert_lab,
                                              overlap->pert_labels,
                                              overlap->pert_max_orders,
                                              &overlap->ket_num_pert,
                                              overlap->ket_pert_labels,
                                              overlap->ket_pert_orders);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertInternTupleToHostLabelOrder()");
    if (overlap->ket_num_pert<0) return QSUCCESS;
    /* performs the same procedure for the perturbations on the overlap operator */
    ierr = RSPPertInternTupleToHostLabelOrder(oper_len_tuple,
                                              oper_pert_tuple,
                                              overlap->num_pert_lab,
                                              overlap->pert_labels,
                                              overlap->pert_max_orders,
                                              &overlap->oper_num_pert,
                                              overlap->oper_pert_labels,
                                              overlap->oper_pert_orders);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertInternTupleToHostLabelOrder()");
    if (overlap->oper_num_pert<0) return QSUCCESS;
    /* calculates expectation values using the callback function */
    overlap->get_overlap_exp(overlap->bra_num_pert,
                             overlap->bra_pert_labels,
                             overlap->bra_pert_orders,
                             overlap->ket_num_pert,
                             overlap->ket_pert_labels,
                             overlap->ket_pert_orders,
                             overlap->oper_num_pert,
                             overlap->oper_pert_labels,
                             overlap->oper_pert_orders,
                             num_dmat,
                             dens_mat,
#if defined(OPENRSP_C_USER_CONTEXT)
                             overlap->user_ctx,
#endif
                             num_exp,
                             val_exp);
    return QSUCCESS;
}

/* <function name='RSPOverlapDestroy'
             attr='private'
             author='Bin Gao'
             date='2014-08-05'>
     Destroys the context of overlap operator, should be called at the end
     <param name='overlap' direction='inout'>
       The context of overlap operator
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPOverlapDestroy(RSPOverlap *overlap)
{
    if (overlap->pert_max_orders!=NULL) {
        free(overlap->pert_max_orders);
        overlap->pert_max_orders = NULL;
    }
    if (overlap->bra_pert_orders!=NULL) {
        free(overlap->bra_pert_orders);
        overlap->bra_pert_orders = NULL;
    }
    if (overlap->ket_pert_orders!=NULL) {
        free(overlap->ket_pert_orders);
        overlap->ket_pert_orders = NULL;
    }
    if (overlap->oper_pert_orders!=NULL) {
        free(overlap->oper_pert_orders);
        overlap->oper_pert_orders = NULL;
    }
    if (overlap->pert_labels!=NULL) {
        free(overlap->pert_labels);
        overlap->pert_labels = NULL;
    }
    if (overlap->bra_pert_labels!=NULL) {
        free(overlap->bra_pert_labels);
        overlap->bra_pert_labels = NULL;
    }
    if (overlap->ket_pert_labels!=NULL) {
        free(overlap->ket_pert_labels);
        overlap->ket_pert_labels = NULL;
    }
    if (overlap->oper_pert_labels!=NULL) {
        free(overlap->oper_pert_labels);
        overlap->oper_pert_labels = NULL;
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    overlap->user_ctx = NULL;
#endif
    overlap->get_overlap_mat = NULL;
    overlap->get_overlap_exp = NULL;
    return QSUCCESS;
}

