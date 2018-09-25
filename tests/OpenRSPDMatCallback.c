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

#include "OpenRSPDMatCallback.h"

void get_overlap_mat(const QInt bra_num_pert,
                     const QcPertInt *bra_pert_labels,
                     const QInt *bra_pert_orders,
                     const QInt ket_num_pert,
                     const QcPertInt *ket_pert_labels,
                     const QInt *ket_pert_orders,
                     const QInt oper_num_pert,
                     const QcPertInt *oper_pert_labels,
                     const QInt *oper_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                     void *user_ctx,
#endif
                     const QInt num_int,
                     QcMat *val_int[])
{
#include "OpenRSPTestPerturbations.h"
#if defined(OPENRSP_C_USER_CONTEXT)
    char *overlap_context;
#endif
    QcMat **tmp_int;
    const QReal NUMBER_ONE[]={1.0,0.0};
    QBool zero_int;
    QErrorCode ierr;
    QInt ilab;
#if defined(OPENRSP_C_USER_CONTEXT)
    overlap_context = (char *)user_ctx;
    if (strcmp(overlap_context, OVERLAP_CONTEXT)!=0) {
        printf("get_overlap_mat>> incorrect operator at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#endif
    switch (oper_num_pert) {
    case 0:
        switch (bra_num_pert) {
        case 0:
            switch (ket_num_pert) {
            case 0:
                /* unperturbed AO-based overlap integral matrix */
                tmp_int = (QcMat **)malloc(sizeof(QcMat *));
                if (tmp_int==NULL) {
                    printf("get_overlap_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "allocates memory for temporary matrix array");
                    exit(QFAILURE);
                }
                tmp_int[0] = (QcMat *)malloc(sizeof(QcMat));
                if (tmp_int[0]==NULL) {
                    printf("get_overlap_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "allocates memory for temporary matrices");
                    exit(QFAILURE);
                }
                ierr = OpenRSPTestReadMat(DMAT_OVERLAP, 0, 1, tmp_int);
                if (ierr!=QSUCCESS) {
                    printf("get_overlap_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling OpenRSPTestReadMat()");
                    exit(ierr);
                }
                ierr = QcMatAXPY(NUMBER_ONE, tmp_int[0], val_int[0]);
                if (ierr!=QSUCCESS) {
                    printf("get_overlap_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatAXPY()");
                    exit(ierr);
                }
                ierr = QcMatDestroy(tmp_int[0]);
                if (ierr!=QSUCCESS) {
                    printf("get_overlap_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatDestroy()");
                    exit(ierr);
                }
                free(tmp_int[0]);
                tmp_int[0] = NULL;
                free(tmp_int);
                tmp_int = NULL;
                break;
            default:
                zero_int = QFALSE;
                for (ilab=0; ilab<ket_num_pert; ilab++) {
                    if (ket_pert_labels[ilab]!=PERT_GEOMETRIC &&
                        ket_pert_labels[ilab]!=PERT_MAGNETIC) {
                        zero_int = QTRUE;
                        break;
                    }
                }
                if (zero_int==QFALSE) {
                    printf("get_overlap_mat>> not implemented at %s\n",
                           FILE_AND_LINE);
                    exit(QFAILURE);
                }
            }
            break;
        default:
            zero_int = QFALSE;
            for (ilab=0; ilab<bra_num_pert; ilab++) {
                if (bra_pert_labels[ilab]!=PERT_GEOMETRIC &&
                    bra_pert_labels[ilab]!=PERT_MAGNETIC) {
                    zero_int = QTRUE;
                    break;
                }
            }
            if (zero_int==QFALSE) {
                printf("get_overlap_mat>> not implemented at %s\n",
                       FILE_AND_LINE);
                exit(QFAILURE);
            }
        }
        break;
    default:
        zero_int = QFALSE;
        for (ilab=0; ilab<oper_num_pert; ilab++) {
            if (oper_pert_labels[ilab]!=PERT_GEOMETRIC &&
                oper_pert_labels[ilab]!=PERT_MAGNETIC) {
                zero_int = QTRUE;
                break;
            }
        }
        if (zero_int==QFALSE) {
            printf("get_overlap_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
}

void get_overlap_exp(const QInt bra_num_pert,
                     const QcPertInt *bra_pert_labels,
                     const QInt *bra_pert_orders,
                     const QInt ket_num_pert,
                     const QcPertInt *ket_pert_labels,
                     const QInt *ket_pert_orders,
                     const QInt oper_num_pert,
                     const QcPertInt *oper_pert_labels,
                     const QInt *oper_pert_orders,
                     const QInt num_dens,
                     QcMat *ao_dens[],
#if defined(OPENRSP_C_USER_CONTEXT)
                     void *user_ctx,
#endif
                     const QInt num_exp,
                     QReal *val_exp)
{
#include "OpenRSPTestPerturbations.h"
#if defined(OPENRSP_C_USER_CONTEXT)
    char *overlap_context;
#endif
    QcMat *val_int[1];
    QInt idens;
    QBool zero_int;
    QErrorCode ierr;
    QInt ilab;
#if defined(OPENRSP_C_USER_CONTEXT)
    overlap_context = (char *)user_ctx;
    if (strcmp(overlap_context, OVERLAP_CONTEXT)!=0) {
        printf("get_overlap_mat>> incorrect operator at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#endif
    switch (oper_num_pert) {
    case 0:
        switch (bra_num_pert) {
        case 0:
            switch (ket_num_pert) {
            case 0:
                /* unperturbed AO-based overlap integral matrix */
                val_int[0] = (QcMat *)malloc(sizeof(QcMat));
                if (val_int[0]==NULL) {
                    printf("get_overlap_exp>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "allocats memory for matrix");
                    exit(QFAILURE);
                }
                ierr = OpenRSPTestReadMat(DMAT_OVERLAP, 0, 1, val_int);
                if (ierr!=QSUCCESS) {
                    printf("get_overlap_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling OpenRSPTestReadMat()");
                    exit(ierr);
                }
                for (idens=0; idens<num_dens; idens++) {
                    ierr = QcMatGetMatProdTrace(val_int[0],
                                                ao_dens[idens],
                                                MAT_NO_OPERATION,
                                                1,
                                                &val_exp[2*idens]);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_exp>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatGetMatProdTrace");
                        exit(ierr);
                    }
                }
                ierr = QcMatDestroy(val_int[0]);
                if (ierr!=QSUCCESS) {
                    printf("get_overlap_exp>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatDestroy()");
                    exit(ierr);
                }
                free(val_int[0]);
                val_int[0] = NULL;
                break;
            default:
                zero_int = QFALSE;
                for (ilab=0; ilab<ket_num_pert; ilab++) {
                    if (ket_pert_labels[ilab]!=PERT_GEOMETRIC &&
                        ket_pert_labels[ilab]!=PERT_MAGNETIC) {
                        zero_int = QTRUE;
                        break;
                    }
                }
                if (zero_int==QFALSE) {
                    printf("get_overlap_exp>> not implemented at %s\n",
                           FILE_AND_LINE);
                    exit(QFAILURE);
                }
            }
            break;
        default:
            zero_int = QFALSE;
            for (ilab=0; ilab<bra_num_pert; ilab++) {
                if (bra_pert_labels[ilab]!=PERT_GEOMETRIC &&
                    bra_pert_labels[ilab]!=PERT_MAGNETIC) {
                    zero_int = QTRUE;
                    break;
                }
            }
            if (zero_int==QFALSE) {
                printf("get_overlap_exp>> not implemented at %s\n",
                       FILE_AND_LINE);
                exit(QFAILURE);
            }
        }
        break;
    default:
        zero_int = QFALSE;
        for (ilab=0; ilab<oper_num_pert; ilab++) {
            if (oper_pert_labels[ilab]!=PERT_GEOMETRIC &&
                oper_pert_labels[ilab]!=PERT_MAGNETIC) {
                zero_int = QTRUE;
                break;
            }
        }
        if (zero_int==QFALSE) {
            printf("get_overlap_exp>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
}

void get_one_oper_mat(const QInt oper_num_pert,
                      const QcPertInt *oper_pert_labels,
                      const QInt *oper_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                      void *user_ctx,
#endif
                      const QInt num_int,
                      QcMat *val_int[])
{
#include "OpenRSPTestPerturbations.h"
#if defined(OPENRSP_C_USER_CONTEXT)
    char *one_oper_context;
#endif
    QcMat **tmp_int;
    const QReal NUMBER_ONE[]={1.0,0.0};
    QBool zero_int;
    QErrorCode ierr;
    QInt ilab;
#if defined(OPENRSP_C_USER_CONTEXT)
    one_oper_context = (char *)user_ctx;
    if (strcmp(one_oper_context, ONEHAM_CONTEXT)==0) {
        zero_int = QFALSE;
        for (ilab=0; ilab<oper_num_pert; ilab++) {
            if (oper_pert_labels[ilab]!=PERT_GEOMETRIC &&
                oper_pert_labels[ilab]!=PERT_MAGNETIC) {
                zero_int = QTRUE;
                break;
            }
        }
        if (zero_int==QFALSE) {
            printf("get_one_oper_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
    else if (strcmp(one_oper_context, EXT_FIELD_CONTEXT)==0) {
        zero_int = QFALSE;
        for (ilab=0; ilab<oper_num_pert; ilab++) {
            if ((oper_pert_labels[ilab]!=PERT_GEOMETRIC &&
                 oper_pert_labels[ilab]!=PERT_MAGNETIC &&
                 oper_pert_labels[ilab]!=PERT_DIPOLE) ||
                (oper_pert_labels[ilab]==PERT_DIPOLE &&
                 oper_pert_orders[0]>1)) {
                zero_int = QTRUE;
                break;
            }
        }
        if (zero_int==QFALSE) {
            switch (oper_num_pert) {
            case 1:
                /* unperturbed AO-based dipole length integral matrices */
                if (oper_pert_labels[0]==PERT_DIPOLE && oper_pert_orders[0]==1) {
                    tmp_int = (QcMat **)malloc(3*sizeof(QcMat *));
                    if (tmp_int==NULL) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "allocates memory for temporary matrix array");
                        exit(QFAILURE);
                    }
                    for (ilab=0; ilab<3; ilab++) {
                        tmp_int[ilab] = (QcMat *)malloc(sizeof(QcMat));
                        if (tmp_int[ilab]==NULL) {
                            printf("get_one_oper_mat>> error happened at %s: %s\n",
                                   FILE_AND_LINE,
                                   "allocates memory for temporary matrices");
                            exit(QFAILURE);
                        }
                    }
                    ierr = OpenRSPTestReadMat(DMAT_DIPLEN, 0, 3, tmp_int);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling OpenRSPTestReadMat()");
                        exit(ierr);
                    }
                    for (ilab=0; ilab<3; ilab++) {
                        ierr = QcMatAXPY(NUMBER_ONE, tmp_int[ilab], val_int[ilab]);
                        if (ierr!=QSUCCESS) {
                            printf("get_one_oper_mat>> error happened at %s: %s\n",
                                   FILE_AND_LINE,
                                   "calling QcMatAXPY()");
                            exit(ierr);
                        }
                        ierr = QcMatDestroy(tmp_int[ilab]);
                        if (ierr!=QSUCCESS) {
                            printf("get_overlap_mat>> error happened at %s: %s\n",
                                   FILE_AND_LINE,
                                   "calling QcMatDestroy()");
                            exit(ierr);
                        }
                        free(tmp_int[ilab]);
                        tmp_int[ilab] = NULL;
                    }
                    free(tmp_int);
                    tmp_int = NULL;
                }
                else {
                    printf("get_one_oper_mat>> not implemented at %s\n",
                           FILE_AND_LINE);
                    exit(QFAILURE);
                }
                break;
            default:
                printf("get_one_oper_mat>> not implemented at %s\n",
                       FILE_AND_LINE);
                exit(QFAILURE);
            }
        }
    }
    else {
        printf("get_one_oper_mat>> unknown one-electron operator at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#else
    zero_int = QFALSE;
    for (ilab=0; ilab<oper_num_pert; ilab++) {
        if ((oper_pert_labels[ilab]!=PERT_GEOMETRIC &&
             oper_pert_labels[ilab]!=PERT_MAGNETIC &&
             oper_pert_labels[ilab]!=PERT_DIPOLE) ||
            (oper_pert_labels[ilab]==PERT_DIPOLE &&
             oper_pert_orders[0]>1)) {
            zero_int = QTRUE;
            break;
        }
    }
    if (zero_int==QFALSE) {
        switch (oper_num_pert) {
        case 1:
            /* unperturbed AO-based dipole length integral matrices */
            if (oper_pert_labels[0]==PERT_DIPOLE && oper_pert_orders[0]==1) {
                tmp_int = (QcMat **)malloc(3*sizeof(QcMat *));
                if (tmp_int==NULL) {
                    printf("get_one_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "allocates memory for temporary matrix array");
                    exit(QFAILURE);
                }
                for (ilab=0; ilab<3; ilab++) {
                    tmp_int[ilab] = (QcMat *)malloc(sizeof(QcMat));
                    if (tmp_int[ilab]==NULL) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "allocates memory for temporary matrices");
                        exit(QFAILURE);
                    }
                }
                ierr = OpenRSPTestReadMat(DMAT_DIPLEN, 0, 3, tmp_int);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling OpenRSPTestReadMat()");
                    exit(ierr);
                }
                for (ilab=0; ilab<3; ilab++) {
                    ierr = QcMatAXPY(NUMBER_ONE, tmp_int[ilab], val_int[ilab]);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatAXPY()");
                        exit(ierr);
                    }
                    ierr = QcMatDestroy(tmp_int[ilab]);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatDestroy()");
                        exit(ierr);
                    }
                    free(tmp_int[ilab]);
                    tmp_int[ilab] = NULL;
                }
                free(tmp_int);
                tmp_int = NULL;
            }
            else {
                printf("get_one_oper_mat>> not implemented at %s\n",
                       FILE_AND_LINE);
                exit(QFAILURE);
            }
            break;
        default:
            printf("get_one_oper_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
#endif
}

void get_one_oper_exp(const QInt oper_num_pert,
                      const QcPertInt *oper_pert_labels,
                      const QInt *oper_pert_orders,
                      const QInt num_dmat,
                      QcMat *dens_mat[],
#if defined(OPENRSP_C_USER_CONTEXT)
                      void *user_ctx,
#endif
                      const QInt num_exp,
                      QReal *val_exp)
{
#include "OpenRSPTestPerturbations.h"
#if defined(OPENRSP_C_USER_CONTEXT)
    char *one_oper_context;
#endif
    QcMat *val_int[1];
    QBool zero_int;
    QInt offset_exp;
    QInt ilab;
    QInt idens;
    QErrorCode ierr;
#if defined(OPENRSP_C_USER_CONTEXT)
    one_oper_context = (char *)user_ctx;
    if (strcmp(one_oper_context, ONEHAM_CONTEXT)==0) {
        zero_int = QFALSE;
        for (ilab=0; ilab<oper_num_pert; ilab++) {
            if (oper_pert_labels[ilab]!=PERT_GEOMETRIC &&
                oper_pert_labels[ilab]!=PERT_MAGNETIC) {
                zero_int = QTRUE;
                break;
            }
        }
        if (zero_int==QFALSE) {
            printf("get_one_oper_exp>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
    else if (strcmp(one_oper_context, EXT_FIELD_CONTEXT)==0) {
        zero_int = QFALSE;
        for (ilab=0; ilab<oper_num_pert; ilab++) {
            if ((oper_pert_labels[ilab]!=PERT_GEOMETRIC &&
                 oper_pert_labels[ilab]!=PERT_MAGNETIC &&
                 oper_pert_labels[ilab]!=PERT_DIPOLE) ||
                (oper_pert_labels[ilab]==PERT_DIPOLE &&
                 oper_pert_orders[0]>1)) {
                zero_int = QTRUE;
                break;
            }
        }
        if (zero_int==QFALSE) {
            switch (oper_num_pert) {
            case 1:
                /* unperturbed AO-based dipole length integral matrices */
                if (oper_pert_labels[0]==PERT_DIPOLE && oper_pert_orders[0]==1) {
                    val_int[0] = (QcMat *)malloc(sizeof(QcMat));
                    if (val_int[0]==NULL) {
                        printf("get_one_oper_exp>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "allocats memory for matrix");
                        exit(QFAILURE);
                    }
                    offset_exp = 0;
                    for (ilab=0; ilab<3; ilab++) {
                        ierr = OpenRSPTestReadMat(DMAT_DIPLEN, ilab, 1, val_int);
                        if (ierr!=QSUCCESS) {
                            printf("get_one_oper_exp>> error happened at %s: %s\n",
                                   FILE_AND_LINE,
                                   "calling OpenRSPTestReadMat()");
                            exit(ierr);
                        }
                        for (idens=0; idens<num_dmat; idens++) {
                            ierr = QcMatGetMatProdTrace(val_int[0],
                                                        dens_mat[idens],
                                                        MAT_NO_OPERATION,
                                                        1,
                                                        &val_exp[offset_exp+2*idens]);
                            if (ierr!=QSUCCESS) {
                                printf("get_one_oper_exp>> error happened at %s: %s\n",
                                       FILE_AND_LINE,
                                       "calling QcMatGetMatProdTrace()");
                                exit(ierr);
                            }
                        }
                        ierr = QcMatDestroy(val_int[0]);
                        if (ierr!=QSUCCESS) {
                            printf("get_one_oper_exp>> error happened at %s: %s\n",
                                   FILE_AND_LINE,
                                   "calling QcMatDestroy()");
                            exit(ierr);
                        }
                        offset_exp += 2*num_dmat;
                    }
                    free(val_int[0]);
                    val_int[0] = NULL;
                }
                else {
                    printf("get_one_oper_exp>> not implemented at %s\n",
                           FILE_AND_LINE);
                    exit(QFAILURE);
                }
                break;
            default:
                printf("get_one_oper_exp>> not implemented at %s\n",
                       FILE_AND_LINE);
                exit(QFAILURE);
            }
        }
    }
    else {
        printf("get_one_oper_exp>> unknown one-electron operator at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#else
    zero_int = QFALSE;
    for (ilab=0; ilab<oper_num_pert; ilab++) {
        if ((oper_pert_labels[ilab]!=PERT_GEOMETRIC &&
             oper_pert_labels[ilab]!=PERT_MAGNETIC &&
             oper_pert_labels[ilab]!=PERT_DIPOLE) ||
            (oper_pert_labels[ilab]==PERT_DIPOLE &&
             oper_pert_orders[0]>1)) {
            zero_int = QTRUE;
            break;
        }
    }
    if (zero_int==QFALSE) {
        switch (oper_num_pert) {
        case 1:
            /* unperturbed AO-based dipole length integral matrices */
            if (oper_pert_labels[0]==PERT_DIPOLE && oper_pert_orders[0]==1) {
                val_int[0] = (QcMat *)malloc(sizeof(QcMat));
                if (val_int[0]==NULL) {
                    printf("get_one_oper_exp>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "allocats memory for matrix");
                    exit(QFAILURE);
                }
                offset_exp = 0;
                for (ilab=0; ilab<3; ilab++) {
                    ierr = OpenRSPTestReadMat(DMAT_DIPLEN, ilab, 1, val_int);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_exp>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling OpenRSPTestReadMat()");
                        exit(ierr);
                    }
                    for (idens=0; idens<num_dmat; idens++) {
                        ierr = QcMatGetMatProdTrace(val_int[0],
                                                    dens_mat[idens],
                                                    MAT_NO_OPERATION,
                                                    1,
                                                    &val_exp[offset_exp+2*idens]);
                        if (ierr!=QSUCCESS) {
                            printf("get_one_oper_exp>> error happened at %s: %s\n",
                                   FILE_AND_LINE,
                                   "calling QcMatGetMatProdTrace()");
                            exit(ierr);
                        }
                    }
                    ierr = QcMatDestroy(val_int[0]);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_exp>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatDestroy()");
                        exit(ierr);
                    }
                    offset_exp += 2*num_dmat;
                }
                free(val_int[0]);
                val_int[0] = NULL;
            }
            else {
                printf("get_one_oper_exp>> not implemented at %s\n",
                       FILE_AND_LINE);
                exit(QFAILURE);
            }
            break;
        default:
            printf("get_one_oper_exp>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
#endif
}

void get_two_oper_mat(const QInt oper_num_pert,
                      const QcPertInt *oper_pert_labels,
                      const QInt *oper_pert_orders,
                      const QInt num_dmat,
                      QcMat *dens_mat[],
#if defined(OPENRSP_C_USER_CONTEXT)
                      void *user_ctx,
#endif
                      const QInt num_int,
                      QcMat *val_int[])
{
#include "OpenRSPTestPerturbations.h"
#if defined(OPENRSP_C_USER_CONTEXT)
    char *two_oper_context;
#endif
    QcMat **tmp_int;
    const QReal NUMBER_ONE[]={1.0,0.0};
    QBool zero_int;
    QErrorCode ierr;
    QInt ilab;
    static QInt id_gmat = -1;
#if defined(OPENRSP_C_USER_CONTEXT)
    two_oper_context = (char *)user_ctx;
    if (strcmp(two_oper_context, TWO_OPER_CONTEXT)!=0) {
        printf("get_two_oper_mat>> not implemented at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#endif
    zero_int = QFALSE;
    for (ilab=0; ilab<oper_num_pert; ilab++) {
        if (oper_pert_labels[ilab]!=PERT_GEOMETRIC &&
            oper_pert_labels[ilab]!=PERT_MAGNETIC) {
            zero_int = QTRUE;
            break;
        }
    }
    if (zero_int==QFALSE) {
        /* two-electron integrals contracted with perturbed AO-based density
           matrices for solving polarizability at frequency as 0.072 */
        if (oper_num_pert==0) {
            id_gmat++;
            if (id_gmat>5) {
                printf("get_two_oper_mat>> not implemented at %s\n",
                       FILE_AND_LINE);
                exit(QFAILURE);
            }
            tmp_int = (QcMat **)malloc(sizeof(QcMat *));
            if (tmp_int==NULL) {
                printf("get_two_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "allocates memory for temporary matrix array");
                exit(QFAILURE);
            }
            tmp_int[0] = (QcMat *)malloc(sizeof(QcMat));
            if (tmp_int[0]==NULL) {
                printf("get_two_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "allocates memory for temporary matrices");
                exit(QFAILURE);
            }
            ierr = OpenRSPTestReadMat(ALPHA_GMAT_AO_HF, id_gmat, 1, tmp_int);
            if (ierr!=QSUCCESS) {
                printf("get_two_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling OpenRSPTestReadMat()");
                exit(ierr);
            }
            ierr = QcMatAXPY(NUMBER_ONE, tmp_int[0], val_int[0]);
            if (ierr!=QSUCCESS) {
                printf("get_two_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatAXPY()");
                exit(ierr);
            }
            ierr = QcMatDestroy(tmp_int[0]);
            if (ierr!=QSUCCESS) {
                printf("get_two_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatDestroy()");
                exit(ierr);
            }
            free(tmp_int[0]);
            tmp_int[0] = NULL;
            free(tmp_int);
            tmp_int = NULL;
        }
        else {
            printf("get_two_oper_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
}

void get_two_oper_exp(const QInt oper_num_pert,
                      const QcPertInt *oper_pert_labels,
                      const QInt *oper_pert_orders,
                      const QInt dmat_len_tuple,
                      const QInt *num_LHS_dmat,
                      QcMat *LHS_dens_mat[],
                      const QInt *num_RHS_dmat,
                      QcMat *RHS_dens_mat[],
#if defined(OPENRSP_C_USER_CONTEXT)
                      void *user_ctx,
#endif
                      const QInt num_exp,
                      QReal *val_exp)
{
#include "OpenRSPTestPerturbations.h"
#if defined(OPENRSP_C_USER_CONTEXT)
    char *two_oper_context;
#endif
    QBool zero_int;
    QErrorCode ierr;
    QInt ilab;
#if defined(OPENRSP_C_USER_CONTEXT)
    two_oper_context = (char *)user_ctx;
    if (strcmp(two_oper_context, TWO_OPER_CONTEXT)!=0) {
        printf("get_two_oper_exp>> not implemented at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#endif
    zero_int = QFALSE;
    for (ilab=0; ilab<oper_num_pert; ilab++) {
        if (oper_pert_labels[ilab]!=PERT_GEOMETRIC &&
            oper_pert_labels[ilab]!=PERT_MAGNETIC) {
            zero_int = QTRUE;
            break;
        }
    }
    if (zero_int==QFALSE) {
        printf("get_two_oper_exp>> not implemented at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
}

void get_linear_rsp_solution(const QInt num_pert,
                             const QInt *num_comps,
                             const QInt *num_freq_sums,
                             const QReal *freq_sums,
                             QcMat *RHS_mat[],
#if defined(OPENRSP_C_USER_CONTEXT)
                             void *user_ctx,
#endif
                             QcMat *rsp_param[])
{
#if defined(OPENRSP_C_USER_CONTEXT)
    char *solver_context;
#endif
    static QInt id_rsparam = -1;
    QErrorCode ierr;
#if defined(OPENRSP_C_USER_CONTEXT)
    solver_context = (char *)user_ctx;
    if (strcmp(solver_context, SOLVER_CONTEXT)!=0) {
        printf("get_linear_rsp_solution>> invalid context at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#endif
    id_rsparam++;
    if (id_rsparam>2) {
        printf("get_linear_rsp_solution>> not implemented at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
    ierr = OpenRSPTestReadMat(ALPHA_RSPARAM_AO_HF, id_rsparam, 1, rsp_param);
    if (ierr!=QSUCCESS) {
        printf("get_linear_rsp_solution>> error happened at %s: %s\n",
               FILE_AND_LINE,
               "calling OpenRSPTestReadMat()");
        exit(ierr);
    }
}

