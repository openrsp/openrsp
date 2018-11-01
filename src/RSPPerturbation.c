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

#include "RSPPerturbation.h"

const QcPertInt OPENRSP_PERT_LABEL_MAX = (1<<OPENRSP_PERT_LABEL_BIT)-1;
const QcPertInt OPENRSP_PERT_ID_MAX = (QCPERTINT_MAX-(1<<OPENRSP_PERT_LABEL_BIT)+1)
                                    >> OPENRSP_PERT_LABEL_BIT;

/* see https://scaryreasoner.wordpress.com/2009/02/28/checking-sizeof-at-compile-time
   accessing date Oct. 6, 2015 */
#define QC_BUILD_BUG_ON(condition) ((void)sizeof(char[1 - 2*!!(condition)]))
void RSPPertCheckLabelBit()
{
    QC_BUILD_BUG_ON(sizeof(QcPertInt)*CHAR_BIT<=OPENRSP_PERT_LABEL_BIT);
    QC_BUILD_BUG_ON(QINT_MAX<OPENRSP_PERT_LABEL_MAX);
}
/* <function name='RSPPertCreate'
             attr='private'
             author='Bin Gao'
             date='2015-06-28'>
     Sets all perturbations involved in response theory calculations, should be
     called at first
     <param name='rsp_pert' direction='inout'>
       The context of perturbations
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels involved in calculations
     </param>
     <param name='pert_labels' direction='in'>
       All the different perturbation labels involved
     </param>
     <param name='pert_max_orders' direction='in'>
       Allowed maximal order of a perturbation described by exactly one of
       the above different labels
     </param>
     <param name='pert_num_comps' direction='in'>
       Number of components of a perturbation described by exactly one of
       the above different labels, up to the allowed maximal order, size
       is therefore the sum of <pert_max_orders>
     </param>
     <param name='user_ctx' direction='in'>
       User-defined callback function context
     </param>
     <param name='get_pert_concatenation' direction='in'>
       User-specified function for getting the ranks of components of
       sub-perturbation tuples (with the same perturbation label) for given
       components of the corresponding concatenated perturbation tuple
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPPertCreate(RSPPert *rsp_pert,
                         const QInt num_pert_lab,
                         const QcPertInt *pert_labels,
                         const QInt *pert_max_orders,
                         const QInt *pert_num_comps,
#if defined(OPENRSP_C_USER_CONTEXT)
                         void *user_ctx,
#endif
                         const GetPertCat get_pert_concatenation)
{
    QInt ilab;    /* incremental recorders over perturbation labels */
    QInt jlab;
    QInt iorder;  /* incremental recorder over orders */
    QInt icomp;   /* incremental recorder over components */
    if (num_pert_lab<1) {
        printf("RSPPertCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        QErrorExit(FILE_AND_LINE, "invalid number of perturbation labels");
    }
    else if (num_pert_lab>OPENRSP_PERT_LABEL_MAX) {
        printf("RSPPertCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        printf("RSPPertCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
               OPENRSP_PERT_LABEL_MAX);
        QErrorExit(FILE_AND_LINE, "too many perturbation labels");
    }
    rsp_pert->num_pert_lab = num_pert_lab;
    rsp_pert->pert_labels = (QcPertInt *)malloc(num_pert_lab*sizeof(QcPertInt));
    if (rsp_pert->pert_labels==NULL) {
        printf("RSPPertCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        QErrorExit(FILE_AND_LINE, "allocates memory for perturbation labels");
    }
    rsp_pert->pert_max_orders = (QInt *)malloc(num_pert_lab*sizeof(QInt));
    if (rsp_pert->pert_max_orders==NULL) {
        printf("RSPPertCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        QErrorExit(FILE_AND_LINE, "allocates memory for allowed maximal orders");
    }
    rsp_pert->ptr_ncomp = (QInt *)malloc((num_pert_lab+1)*sizeof(QInt));
    if (rsp_pert->ptr_ncomp==NULL) {
        printf("RSPPertCreate>> number of perturbation labels %"QINT_FMT"\n",
               num_pert_lab);
        QErrorExit(FILE_AND_LINE, "allocates memory for pointers to components");
    }
    rsp_pert->ptr_ncomp[0] = 0;
    for (ilab=0; ilab<num_pert_lab; ilab++) {
        if (pert_labels[ilab]>OPENRSP_PERT_LABEL_MAX) {
            printf("RSPPertCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                   ilab,
                   pert_labels[ilab]);
            printf("RSPPertCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
                   OPENRSP_PERT_LABEL_MAX);
            QErrorExit(FILE_AND_LINE, "invalid perturbation label");
        }
        /* each element of <pert_labels> should be unique */
        for (jlab=0; jlab<ilab; jlab++) {
            if (pert_labels[jlab]==pert_labels[ilab]) {
                printf("RSPPertCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       jlab,
                       pert_labels[jlab]);
                printf("RSPPertCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ilab,
                       pert_labels[ilab]);
                QErrorExit(FILE_AND_LINE, "repeated perturbation labels not allowed");
            }
        }
        rsp_pert->pert_labels[ilab] = pert_labels[ilab];
        if (pert_max_orders[ilab]<1) {
            printf("RSPPertCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                   ilab,
                   pert_labels[ilab]);
            printf("RSPPertCreate>> allowed maximal order is %"QINT_FMT"\n",
                   pert_max_orders[ilab]);
            QErrorExit(FILE_AND_LINE, "only positive order allowed");
        }
        rsp_pert->pert_max_orders[ilab] = pert_max_orders[ilab];
        /* <c>rsp_pert->ptr_ncomp[ilab]</c> points to the number of components
           of <c>rsp_pert->pert_labels[ilab]</c> */
        rsp_pert->ptr_ncomp[ilab+1] = rsp_pert->ptr_ncomp[ilab]+pert_max_orders[ilab];
    }
    /* <c>rsp_pert->ptr_ncomp[num_pert_lab]</c> equals to the size of
       <c>rsp_pert->pert_num_comps</c> */
    rsp_pert->pert_num_comps = (QInt *)malloc(rsp_pert->ptr_ncomp[num_pert_lab]
                                              *sizeof(QInt));
    if (rsp_pert->pert_num_comps==NULL) {
        printf("RSPPertCreate>> size of numbers of components %"QINT_FMT"\n",
               rsp_pert->ptr_ncomp[num_pert_lab]);
        QErrorExit(FILE_AND_LINE, "allocates memory for numbers of components");
    }
    for (ilab=0,icomp=0; ilab<num_pert_lab; ilab++) {
        for (iorder=1; iorder<=rsp_pert->pert_max_orders[ilab]; iorder++,icomp++) {
            if (pert_num_comps[icomp]<1) {
                printf("RSPPertCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ilab,
                       pert_labels[ilab]);
                printf("RSPPertCreate>> allowed maximal order is %"QINT_FMT"\n",
                       pert_max_orders[ilab]);
                printf("RSPPertCreate>> %"QINT_FMT"-th No. of comps. is %"QINT_FMT"\n",
                       iorder,
                       pert_num_comps[icomp]);
                QErrorExit(FILE_AND_LINE, "incorrect number of components");
            }
            rsp_pert->pert_num_comps[icomp] = pert_num_comps[icomp];
        }
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    rsp_pert->user_ctx = user_ctx;
#endif
    rsp_pert->get_pert_concatenation = get_pert_concatenation;
    return QSUCCESS;
}

/* <function name='RSPPertAssemble'
             attr='private'
             author='Bin Gao'
             date='2015-06-28'>
     Assembles the context of perturbations involved in calculations
     <param name='rsp_pert' direction='inout'>
       The context of perturbations
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPPertAssemble(RSPPert *rsp_pert)
{
    if (rsp_pert->pert_labels==NULL ||
        rsp_pert->pert_max_orders==NULL ||
        rsp_pert->ptr_ncomp==NULL ||
        rsp_pert->pert_num_comps==NULL ||
        rsp_pert->get_pert_concatenation==NULL) {
        QErrorExit(FILE_AND_LINE, "perturbations are not correctly set");
    }
    return QSUCCESS;
}

/* <function name='RSPPertWrite'
             attr='private'
             author='Bin Gao'
             date='2015-06-28'>
     Writes the context of perturbations involved in calculations
     <param name='rsp_pert' direction='in'>
       The context of perturbations
     </param>
     <param name='fp_pert' direction='inout'>File pointer</param>
     <return>Error information</return>
   </function> */
QErrorCode RSPPertWrite(const RSPPert *rsp_pert, FILE *fp_pert)
{
    QInt ilab;   /* incremental recorder over perturbation labels */
    QInt icomp;  /* incremental recorder over components */
    fprintf(fp_pert,
            "RSPPertWrite>> number of all perturbation labels %"QINT_FMT"\n",
            rsp_pert->num_pert_lab);
    fprintf(fp_pert,
            "RSPPertWrite>> label           maximum-order    numbers-of-components\n");
    for (ilab=0; ilab<rsp_pert->num_pert_lab; ilab++) {
        fprintf(fp_pert,
                "RSPPertWrite>>  %"QCPERTINT_FMT"               %"QINT_FMT"               ",
                rsp_pert->pert_labels[ilab],
                rsp_pert->pert_max_orders[ilab]);
        for (icomp=rsp_pert->ptr_ncomp[ilab]; icomp<rsp_pert->ptr_ncomp[ilab+1]; icomp++) {
            fprintf(fp_pert, " %"QINT_FMT"", rsp_pert->pert_num_comps[icomp]);
        }
        fprintf(fp_pert, "\n");
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    if (rsp_pert->user_ctx!=NULL) {
        fprintf(fp_pert, "RSPPertWrite>> user-defined function context given\n");
    }
#endif
    return QSUCCESS;
}

/* <function name='RSPPertDestroy'
             attr='private'
             author='Bin Gao'
             date='2015-06-28'>
     Destroys the context of perturbations involved in calculations, should be
     called at the end
     <param name='rsp_pert' direction='inout'>
       The context of perturbations
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPPertDestroy(RSPPert *rsp_pert)
{
    rsp_pert->num_pert_lab = 0;
    free(rsp_pert->pert_labels);
    rsp_pert->pert_labels = NULL;
    free(rsp_pert->pert_max_orders);
    rsp_pert->pert_max_orders = NULL;
    free(rsp_pert->ptr_ncomp);
    rsp_pert->ptr_ncomp = NULL;
    free(rsp_pert->pert_num_comps);
    rsp_pert->pert_num_comps = NULL;
#if defined(OPENRSP_C_USER_CONTEXT)
    rsp_pert->user_ctx = NULL;
#endif
    rsp_pert->get_pert_concatenation = NULL;
    return QSUCCESS;
}

/* <function name='RSPPertValidateLabelOrder'
             attr='private'
             author='Bin Gao'
             date='2015-10-15'>
     Check the validity of given perturbation labels and corresponding allowed
     maximal orders
     <param name='rsp_pert' direction='in'>
       The context of perturbations
     </param>
     <param name='num_pert_lab' direction='in'>
       Number of all different perturbation labels
     </param>
     <param name='pert_labels' direction='in'>
       Different perturbation labels that will be checked
     </param>
     <param name='pert_max_orders' direction='in'>
       Allowed maximal order of a perturbation described by exactly one of
       the above different labels
     </param>
     <return>
       <QSUCCESS> if perturbation labels and orders are valid,
       <QFAILURE> otherwise
     </return>
   </function> */
QErrorCode RSPPertValidateLabelOrder(const RSPPert *rsp_pert,
                                     const QInt num_pert_lab,
                                     const QcPertInt *pert_labels,
                                     const QInt *pert_max_orders)
{
    QInt ilab;         /* incremental recorders over perturbation labels */
    QInt jlab;
    QBool pert_valid;  /* validity of the perturbations */
    for (ilab=0; ilab<num_pert_lab; ilab++) {
        pert_valid = QFALSE;
        for (jlab=0; jlab<rsp_pert->num_pert_lab; jlab++) {
            /* valid perturbation label, checks its allowed maximal order */
            if (pert_labels[ilab]==rsp_pert->pert_labels[jlab]) {
                if (pert_max_orders[ilab]<=rsp_pert->pert_max_orders[jlab]) {
                    pert_valid = QTRUE;
                }
                break;
            }
        }
        if (pert_valid==QFALSE) return QFAILURE;
    }
    return QSUCCESS;
}

/* <function name='RSPPertHostToInternal'
             attr='private'
             author='Bin Gao'
             date='2015-10-08'>
     Check, convert and sort a host program's perturbation tuple and
     corresponding frequencies
     <param name='rsp_pert' direction='in'>
       The context of perturbations
     </param>
     <param name='len_tuple' direction='in'>
       Length of the host program's and the internal perturbation tuples
     </param>
     <param name='pert_tuple' direction='in'>
       The host program's perturbation tuple, in which the first label
       is the perturbation $a$
     </param>
     <param name='num_freq_configs' direction='in'>
       Number of different frequency configurations
     </param>
     <param name='pert_freqs' direction='in'>
       Complex frequencies of each perturbation label (except for the
       perturbation $a$) over all frequency configurations, size is therefore
       $2\times[(<len_tuple>-1)\times<num_freq_configs>]$, and arranged as
       <c>[num_freq_configs][len_tuple-1][2]</c> in memory (that is, the real
       and imaginary parts of each frequency are consecutive in memory)
     </param>
     <param name='intern_pert_tuple' direction='out'>
       The internal perturbation tuple, in which identical perturbation labels
       are consecutive, and the first one is the perturbation $a$
     </param>
     <param name='intern_pert_freqs' direction='out'>
       Internal complex frequencies (in ascending order among identical
       perturbation labels) of each perturbation label (including the
       perturbation $a$) over all frequency configurations, size is therefore
       $2\times<len_tuple>\times<num_freq_configs>$, and arranged as
       <c>[num_freq_configs][len_tuple][2]</c> in memory (that is, the real and
       imaginary parts of each frequency are consecutive in memory)
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPPertHostToInternal(const RSPPert *rsp_pert,
                                 const QInt len_tuple,
                                 QcPertInt *pert_tuple,
                                 const QInt num_freq_configs,
                                 QReal *pert_freqs)
{
//    QcPertInt id_pert;  /* identifiers of different perturbation labels */
//    QBool lab_valid;    /* validity of the perturbation labels */
//
//    QInt ipert,jpert;  /* incremental recorders */
//    QInt first_id;     /* first identical pertubation label in the tuple */
//    QInt last_id;      /* last identical pertubation label in the tuple */
//    QBool non_id;      /* indicates if non-identical label found */
//
//    id_pert = 0;
//
//    for (ipert=0,jpert=0; ipert<len_tuple; ) {
//
//      /* checks the current perturbation label against all known
//           perturbation labels */
//        lab_valid = QFALSE;
//        for (ilab=0; ilab<rsp_pert->num_pert_lab; ilab++) {
//            if (pert_tuple[ipert]==rsp_pert->pert_labels[ilab]) {
//                lab_valid = QTRUE;
//                break;
//            }
//        }
//        if (lab_valid==QTRUE) {
//            /* converts the current perturbation label to internal one */
//            intern_pert_tuple[jpert] = (id_pert<<OPENRSP_PERT_LABEL_BIT)
//                                     + pert_tuple[ipert];
//            /* finds the other same perturbation labels */
//
//
//            /* updates the identifier */
//            id_pert++;
//        else {
//            printf("RSPPertCreate>> %"QINT_FMT"-th perturbation %"QCPERTINT_FMT"\n",
//                   ipert,
//                   pert_tuple[ipert]);
//            QErrorExit(FILE_AND_LINE, "invalid perturbation label");
//        }
//
//
//    }
//
//    /* we first try to find consecutive identical pertubation labels */
//    first_id = 0;
//    non_id = QFALSE;
//    for (ipert=first_id; ipert<len_tuple-1; ipert++) {
//        if (pert_tuple[ipert]!=pert_tuple[ipert+1]) {
//            last_id = ipert;
//            non_id = QTRUE;
//            break;
//        }
//    }
//    if (non_id=QTRUE) {
//    }
//    else {
//    }

    return QSUCCESS;
}

/* <function name='RSPPertInternTupleToHostLabelOrder'
             attr='private'
             author='Bin Gao'
             date='2015-10-08'>
     Convert an internal perturbation tuple to host program's pertubation
     labels and the corresponding orders
     <param name='len_intern_tuple' direction='in'>
       Length of the internal perturbation tuple
     </param>
     <param name='intern_pert_tuple' direction='in'>
       The internal perturbation tuple
     </param>
     <param name='num_allowed_labels' direction='in'>
       Number of allowed different host program's perturbation labels
     </param>
     <param name='allowed_pert_labels' direction='in'>
       All the allowed different host program's perturbation labels
     </param>
     <param name='allowed_max_orders' direction='in'>
       Allowed maximal order of a perturbation described by exactly one of
       the above different host program's labels
     </param>
     <param name='num_pert' direction='out'>
       Number of different perturbations from the internal perturbation tuple,
       $-1$ indicates there are perturbation labels/orders not allowed
     </param>
     <param name='pert_labels' direction='out'>
       Host program's pertubation labels of the resulted perturbations
     </param>
     <param name='pert_orders' direction='out'>
       Orders of the resulted perturbations
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPPertInternTupleToHostLabelOrder(const QInt len_intern_tuple,
                                              const QcPertInt *intern_pert_tuple,
                                              const QInt num_allowed_labels,
                                              const QcPertInt *allowed_pert_labels,
                                              const QInt *allowed_max_orders,
                                              QInt *num_pert,
                                              QcPertInt *pert_labels,
                                              QInt *pert_orders)
{
    QcPertInt host_pert_label;  /* host program's perturbation label */
    QBool pert_allowed;         /* resulted perturbations allowed or not */
    QInt ipert;                 /* incremental recorder for different perturbations */
    QInt ilab;                  /* incremental recorder over perturbation labels */
    QInt idx_allowed;           /* index of the allowed perturbation label */
    ipert = 0;
    for (ilab=0; ilab<len_intern_tuple; ) {
        /* converts to the host program's label */
        host_pert_label = intern_pert_tuple[ilab] & OPENRSP_PERT_LABEL_MAX;
        /* checks if the label is allowed */
        pert_allowed = QFALSE;
        for (idx_allowed=0; idx_allowed<num_allowed_labels; idx_allowed++) {
            if (host_pert_label==allowed_pert_labels[idx_allowed]) {
                pert_allowed = QTRUE;
                break;
            }
        }
        /* returns a negative number if the label is not allowed */
        if (pert_allowed==QFALSE) {
            *num_pert = -1;
            return QSUCCESS;
        }
        /* finds consecutive identical internal perturbation labels */
        pert_labels[ipert] = intern_pert_tuple[ilab];
        pert_orders[ipert] = 1;
        ilab++;
        for (; ilab<len_intern_tuple; ) {
            if (pert_labels[ipert]==intern_pert_tuple[ilab]) {
                pert_orders[ipert]++;
            }
            else {
                break;
            }
            ilab++;
        }
        /* checks if the order is allowed */
        if (pert_orders[ipert]>allowed_max_orders[idx_allowed]) {
            *num_pert = -1;
            return QSUCCESS;
        }
        /* saves the host program's label */
        pert_labels[ipert] = host_pert_label;
        ipert++;
    }
    *num_pert = ipert;
    return QSUCCESS;
}

/* <function name='RSPPertInternTupleToHostTuple'
             attr='private'
             author='Bin Gao'
             date='2015-10-12'>
     Convert an internal perturbation tuple to the corresponding host
     program's one
     <param name='len_intern_tuple' direction='in'>
       Length of the internal perturbation tuple
     </param>
     <param name='intern_pert_tuple' direction='in'>
       The internal perturbation tuple
     </param>
     <param name='num_allowed_labels' direction='in'>
       Number of allowed different host program's perturbation labels
     </param>
     <param name='allowed_pert_labels' direction='in'>
       All the allowed different host program's perturbation labels
     </param>
     <param name='allowed_max_orders' direction='in'>
       Allowed maximal order of a perturbation described by exactly one of
       the above different host program's labels
     </param>
     <param name='len_tuple' direction='out'>
       Length of the resulted host program's perturbation tuple,
       $-1$ indicates there are perturbation labels/orders not allowed
     </param>
     <param name='pert_tuple' direction='out'>
       The resulted host program's perturbation tuple
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPPertInternTupleToHostTuple(const QInt len_intern_tuple,
                                         const QcPertInt *intern_pert_tuple,
                                         const QInt num_allowed_labels,
                                         const QcPertInt *allowed_pert_labels,
                                         const QInt *allowed_max_orders,
                                         QInt *len_tuple,
                                         QcPertInt *pert_tuple)
{
    QcPertInt host_pert_label;  /* host program's perturbation label */
    QInt host_pert_order;       /* order of host program's perturbation */
    QBool pert_allowed;         /* resulted perturbations allowed or not */
    QInt ipert;                 /* incremental recorder for different perturbations */
    QInt ilab;                  /* incremental recorder over perturbation labels */
    QInt idx_allowed;           /* index of the allowed perturbation label */
    ipert = 0;
    for (ilab=0; ilab<len_intern_tuple; ) {
        /* converts to the host program's label */
        host_pert_label = intern_pert_tuple[ilab] & OPENRSP_PERT_LABEL_MAX;
        /* checks if the label is allowed */
        pert_allowed = QFALSE;
        for (idx_allowed=0; idx_allowed<num_allowed_labels; idx_allowed++) {
            if (host_pert_label==allowed_pert_labels[idx_allowed]) {
                pert_allowed = QTRUE;
                break;
            }
        }
        /* returns a negative number if the label is not allowed */
        if (pert_allowed==QFALSE) {
            *len_tuple = -1;
            return QSUCCESS;
        }
        /* finds consecutive identical internal perturbation labels */
        pert_tuple[ipert] = intern_pert_tuple[ilab];
        host_pert_order = 1;
        ilab++;
        for (; ilab<len_intern_tuple; ) {
            if (pert_tuple[ipert]==intern_pert_tuple[ilab]) {
                host_pert_order++;
            }
            else {
                break;
            }
            ilab++;
        }
        /* checks if the order is allowed */
        if (host_pert_order>allowed_max_orders[idx_allowed]) {
            *len_tuple = -1;
            return QSUCCESS;
        }
        /* saves the host program's labels */
        for (; host_pert_order>0; host_pert_order--) {
            pert_tuple[ipert] = host_pert_label;
            ipert++;
        }
    }
    *len_tuple = ipert;
    return QSUCCESS;
}

/* <function name='RSPPertGetConcatenation'
             attr='private'
             author='Bin Gao'
             date='2015-06-28'>
     Gets the ranks of components of sub-perturbation tuples
     <param name='rsp_pert' direction='inout'>
       The context of perturbations
     </param>
     <param name='intern_pert_label' direction='in'>
       The internal perturbation label
     </param>
     <param name='first_cat_comp' direction='in'>
       Rank of the first component of the concatenated perturbation tuple
     </param>
     <param name='num_cat_comps' direction='in'>
       Number of components of the concatenated perturbation tuple
     </param>
     <param name='num_sub_tuples' direction='in'>
       Number of sub-perturbation tuples to construct the concatenated
       perturbation tuple
     </param>
     <param name='len_sub_tuples' direction='in'>
       Length of each sub-perturbation tuple, size is <num_sub_tuples> so
       that the length of the concatenated perturbation tuple is the sum
       of <len_sub_tuples>
     </param>
     <param name='rank_sub_comps' direction='out'>
       Ranks of components of sub-perturbation tuples for the corresponding
       component of the concatenated perturbation tuple, i.e. <num_cat_comps>
       components starting from the one with the rank <first_cat_comp>; size
       is therefore the product of <num_sub_tuples> and <num_cat_comps>, and
       is arranged as <c>[num_cat_comps][num_sub_tuples]</c> in memory
     </param>
     <return>Error information</return>
   </function> */
QErrorCode RSPPertGetConcatenation(const RSPPert *rsp_pert,
                                   const QcPertInt intern_pert_label,
                                   const QInt first_cat_comp,
                                   const QInt num_cat_comps,
                                   const QInt num_sub_tuples,
                                   const QInt *len_sub_tuples,
                                   QInt *rank_sub_comps)
{
    QcPertInt pert_label;
    /* converts to host program's perturbation label */
    pert_label = intern_pert_label & OPENRSP_PERT_LABEL_MAX;
#if defined(OPENRSP_ZERO_BASED)
    rsp_pert->get_pert_concatenation(pert_label,
                                     first_cat_comp,
                                     num_cat_comps,
                                     num_sub_tuples,
                                     len_sub_tuples,
#if defined(OPENRSP_C_USER_CONTEXT)
                                     rsp_pert->user_ctx,
#endif
                                     rank_sub_comps);
#else
    QInt icomp;  /* incremental recorder over ranks of components */
    rsp_pert->get_pert_concatenation(pert_label,
                                     first_cat_comp+1,
                                     num_cat_comps,
                                     num_sub_tuples,
                                     len_sub_tuples,
#if defined(OPENRSP_C_USER_CONTEXT)
                                     rsp_pert->user_ctx,
#endif
                                     rank_sub_comps);
    for (icomp=0; icomp<num_cat_comps*num_sub_tuples; icomp++) {
        rank_sub_comps[icomp]--;
    }
#endif
    return QSUCCESS;
}

