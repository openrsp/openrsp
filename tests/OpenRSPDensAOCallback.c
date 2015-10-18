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

#include "OpenRSPDensAOCallback.h"

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
    QErrorCode ierr;
#if defined(OPENRSP_C_USER_CONTEXT)
    overlap_context = (char *)user_ctx;
    if (strcmp(overlap_context, OVERLAP_CONTEXT)!=0) {
        printf("get_overlap_mat>> incorrect operator at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#endif
    /* unperturbed AO-based overlap integral matrix */
    if (oper_num_pert==0 && bra_num_pert==0 && ket_num_pert==0) {
        ierr = OpenRSPTestReadMat(OVERLAP_AO, 0, 1, val_int);
        if (ierr!=QSUCCESS) {
            printf("get_overlap_mat>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling OpenRSPTestReadMat()");
            exit(ierr);
        }
    }
    else if (oper_num_pert==1 && bra_num_pert==0 && ket_num_pert==0) {
        if (oper_pert_labels[0]==PERT_GEOMETRIC) {
            printf("get_overlap_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        /* zero matrices */
        else if (oper_pert_labels[0]==PERT_DIPOLE) {
            ierr = OpenRSPTestZeroMat(num_int, val_int);
            if (ierr!=QSUCCESS) {
                printf("get_overlap_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling OpenRSPTestZeroMat()");
                exit(ierr);
            }
        }
        else if (oper_pert_labels[0]==PERT_MAGNETIC) {
            printf("get_overlap_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        else {
            printf("get_overlap_mat>> unknown perturbations at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
    else if (oper_num_pert==0 && bra_num_pert==1 && ket_num_pert==0) {
        if (bra_pert_labels[0]==PERT_GEOMETRIC) {
            printf("get_overlap_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        /* zero matrices */
        else if (bra_pert_labels[0]==PERT_DIPOLE) {
            ierr = OpenRSPTestZeroMat(num_int, val_int);
            if (ierr!=QSUCCESS) {
                printf("get_overlap_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling OpenRSPTestZeroMat()");
                exit(ierr);
            }
        }
        else if (bra_pert_labels[0]==PERT_MAGNETIC) {
            printf("get_overlap_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        else {
            printf("get_overlap_mat>> unknown perturbations at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
    else if (oper_num_pert==0 && bra_num_pert==0 && ket_num_pert==1) {
        if (ket_pert_labels[0]==PERT_GEOMETRIC) {
            printf("get_overlap_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        /* zero matrices */
        else if (ket_pert_labels[0]==PERT_DIPOLE) {
            ierr = OpenRSPTestZeroMat(num_int, val_int);
            if (ierr!=QSUCCESS) {
                printf("get_overlap_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling OpenRSPTestZeroMat()");
                exit(ierr);
            }
        }
        else if (ket_pert_labels[0]==PERT_MAGNETIC) {
            printf("get_overlap_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        else {
            printf("get_overlap_mat>> unknown perturbations at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
    else {
        printf("get_overlap_mat>> not implemented at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
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
    QErrorCode ierr;
#if defined(OPENRSP_C_USER_CONTEXT)
    overlap_context = (char *)user_ctx;
    if (strcmp(overlap_context, OVERLAP_CONTEXT)!=0) {
        printf("get_overlap_mat>> incorrect operator at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#endif
    /* unperturbed AO-based overlap integral matrix */
    if (oper_num_pert==0 && bra_num_pert==0 && ket_num_pert==0) {
        val_int[0] = (QcMat *)malloc(sizeof(QcMat));
        if (val_int[0]==NULL) {
            printf("get_overlap_exp>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "allocats memory for matrix");
            exit(QFAILURE);
        }
        ierr = QcMatCreate(val_int[0]);
        if (ierr!=QSUCCESS) {
            printf("get_overlap_exp>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatCreate()");
            exit(ierr);
        }
        ierr = OpenRSPTestReadMat(OVERLAP_AO, 0, 1, val_int);
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
    }
    else if (oper_num_pert==1 && bra_num_pert==0 && ket_num_pert==0) {
        if (oper_pert_labels[0]==PERT_GEOMETRIC) {
            printf("get_overlap_exp>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        /* zero matrices */
        else if (oper_pert_labels[0]==PERT_DIPOLE) {
            for (idens=0; idens<2*num_exp; idens++) {
                val_exp[idens] = 0;
            }
        }
        else if (oper_pert_labels[0]==PERT_MAGNETIC) {
            printf("get_overlap_exp>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        else {
            printf("get_overlap_exp>> unknown perturbations at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
    else if (oper_num_pert==0 && bra_num_pert==1 && ket_num_pert==0) {
        if (bra_pert_labels[0]==PERT_GEOMETRIC) {
            printf("get_overlap_exp>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        /* zero matrices */
        else if (bra_pert_labels[0]==PERT_DIPOLE) {
            for (idens=0; idens<2*num_exp; idens++) {
                val_exp[idens] = 0;
            }
        }
        else if (bra_pert_labels[0]==PERT_MAGNETIC) {
            printf("get_overlap_exp>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        else {
            printf("get_overlap_exp>> unknown perturbations at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
    else if (oper_num_pert==0 && bra_num_pert==0 && ket_num_pert==1) {
        if (ket_pert_labels[0]==PERT_GEOMETRIC) {
            printf("get_overlap_exp>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        /* zero matrices */
        else if (ket_pert_labels[0]==PERT_DIPOLE) {
            for (idens=0; idens<2*num_exp; idens++) {
                val_exp[idens] = 0;
            }
        }
        else if (ket_pert_labels[0]==PERT_MAGNETIC) {
            printf("get_overlap_exp>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        else {
            printf("get_overlap_exp>> unknown perturbations at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
    else if (oper_num_pert==0 && bra_num_pert==1 && ket_num_pert==1) {
        /* zero matrices */
        if (bra_pert_labels[0]==PERT_DIPOLE || ket_pert_labels[0]==PERT_DIPOLE) {
            for (idens=0; idens<2*num_exp; idens++) {
                val_exp[idens] = 0;
            }
        }
        else {
            printf("get_overlap_exp>> unknown perturbations at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
    else {
         printf("get_overlap_exp>> not implemented at %s\n",
                FILE_AND_LINE);
         exit(QFAILURE);
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
    QErrorCode ierr;
#if defined(OPENRSP_C_USER_CONTEXT)
    one_oper_context = (char *)user_ctx;
    if (strcmp(one_oper_context, ONEHAM_CONTEXT)==0) {
        /* zero matrices */
        if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
            ierr = OpenRSPTestZeroMat(num_int, val_int);
            if (ierr!=QSUCCESS) {
                printf("get_one_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling OpenRSPTestZeroMat()");
                exit(ierr);
            }
        }
        else {
            printf("get_one_oper_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
    else if (strcmp(one_oper_context, EXT_FIELD_CONTEXT)==0) {
        if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
            /* unperturbed AO-based dipole length integral matrices */
            if (oper_pert_orders[0]==1) {
                ierr = OpenRSPTestReadMat(DIPLEN_AO, 0, 3, val_int);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling OpenRSPTestReadMat()");
                    exit(ierr);
                }
            }
            /* zero matrices */
            else {
                ierr = OpenRSPTestZeroMat(num_int, val_int);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling OpenRSPTestZeroMat()");
                    exit(ierr);
                }
            }
        }
        else {
            printf("get_one_oper_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
    else {
        printf("get_one_oper_mat>> unknown one-electron operator at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#else
    if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
        /* unperturbed AO-based dipole length integral matrices */
        if (oper_pert_orders[0]==1) {
            ierr = OpenRSPTestReadMat(DIPLEN_AO, 0, 3, val_int);
            if (ierr!=QSUCCESS) {
                printf("get_one_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling OpenRSPTestReadMat()");
                exit(ierr);
            }
        }
        /* zero matrices */
        else {
            ierr = OpenRSPTestZeroMat(num_int, val_int);
            if (ierr!=QSUCCESS) {
                printf("get_one_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling OpenRSPTestZeroMat()");
                exit(ierr);
            }
        }
    }
    else {
        printf("get_one_oper_mat>> not implemented at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
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
    QInt imat;
    QInt offset_exp;
    QInt idens;
    QErrorCode ierr;
#if defined(OPENRSP_C_USER_CONTEXT)
    one_oper_context = (char *)user_ctx;
    if (strcmp(one_oper_context, ONEHAM_CONTEXT)==0) {
        /* zero matrices */
        if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
            for (idens=0; idens<2*num_exp; idens++) {
                val_exp[idens] = 0;
            }
        }
        else {
            printf("get_one_oper_exp>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
    else if (strcmp(one_oper_context, EXT_FIELD_CONTEXT)==0) {
        if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
            /* unperturbed AO-based dipole length integral matrices */
            if (oper_pert_orders[0]==1) {
                val_int[0] = (QcMat *)malloc(sizeof(QcMat));
                if (val_int[0]==NULL) {
                    printf("get_one_oper_exp>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "allocats memory for matrix");
                    exit(QFAILURE);
                }
                offset_exp = 0;
                for (imat=0; imat<3; imat++) {
                    ierr = QcMatCreate(val_int[0]);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_exp>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatCreate()");
                        exit(ierr);
                    }
                    ierr = OpenRSPTestReadMat(DIPLEN_AO, imat, 1, val_int);
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
            /* zero matrices */
            else {
                for (idens=0; idens<2*num_exp; idens++) {
                    val_exp[idens] = 0;
                }
            }
        }
        else {
            printf("get_one_oper_exp>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
    else {
        printf("get_one_oper_exp>> unknown one-electron operator at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#else
    if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
        /* unperturbed AO-based dipole length integral matrices */
        if (oper_pert_orders[0]==1) {
            offset_exp = 0;
            for (imat=0; imat<3; imat++) {
                ierr = QcMatCreate(&val_int[0]);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_exp>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatCreate()");
                    exit(ierr);
                }
                ierr = OpenRSPTestReadMat(DIPLEN_AO, imat, 1, val_int);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_exp>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling OpenRSPTestReadMat()");
                    exit(ierr);
                }
                for (idens=0; idens<num_dmat; idens++) {
                    ierr = QcMatGetMatProdTrace(&val_int[0],
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
        }
        /* zero matrices */
        else {
            for (idens=0; idens<2*num_exp; idens++) {
                val_exp[idens] = 0;
            }
        }
    else {
        printf("get_one_oper_exp>> not implemented at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
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
    QErrorCode ierr;
    static QInt id_gmat = -1;
#if defined(OPENRSP_C_USER_CONTEXT)
    two_oper_context = (char *)user_ctx;
    if (strcmp(two_oper_context, TWO_OPER_CONTEXT)!=0) {
        printf("get_two_oper_mat>> not implemented at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
    else {
        /* two-electron integrals contracted with perturbed AO-based density
           matrices for solving polarizability at frequency as 0.072 */
        if (oper_num_pert==0) {
            id_gmat++;
            if (id_gmat>5) {
                printf("get_two_oper_mat>> not implemented at %s\n",
                       FILE_AND_LINE);
                exit(QFAILURE);
            }
            ierr = OpenRSPTestReadMat(ALPHA_GMAT_AO_HF, id_gmat, 1, val_int);
            if (ierr!=QSUCCESS) {
                printf("get_two_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling OpenRSPTestReadMat()");
                exit(ierr);
            }
        }
        /* zero matrices */
        else if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
            ierr = OpenRSPTestZeroMat(num_int, val_int);
            if (ierr!=QSUCCESS) {
                printf("get_two_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling OpenRSPTestZeroMat()");
                exit(ierr);
            }
        }
        else {
            printf("get_two_oper_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
#else
    /* two-electron integrals contracted with perturbed AO-based density
       matrices for solving polarizability at frequency as 0.072 */
    if (oper_num_pert==0) {
        id_gmat++;
        if (id_gmat>5) {
            printf("get_two_oper_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        ierr = OpenRSPTestReadMat(ALPHA_GMAT_AO_HF, id_gmat, 1, val_int);
        if (ierr!=QSUCCESS) {
            printf("get_two_oper_mat>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling OpenRSPTestReadMat()");
            exit(ierr);
        }
    }
    /* zero matrices */
    else if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
        ierr = OpenRSPTestZeroMat(num_int, val_int);
        if (ierr!=QSUCCESS) {
            printf("get_two_oper_mat>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling OpenRSPTestZeroMat()");
            exit(ierr);
        }
    }
    else {
        printf("get_two_oper_mat>> not implemented at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#endif
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
    QInt ival;
    QErrorCode ierr;
#if defined(OPENRSP_C_USER_CONTEXT)
    two_oper_context = (char *)user_ctx;
    if (strcmp(two_oper_context, TWO_OPER_CONTEXT)!=0) {
        printf("get_two_oper_exp>> not implemented at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
    else {
        /* zero matrices */
        if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
            for (ival=0; ival<2*num_exp; ival++) {
                val_exp[ival] = 0;
            }
        }
        else {
            printf("get_two_oper_exp>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
#else
    /* zero matrices */
    if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
        for (ival=0; ival<2*num_exp; ival++) {
            val_exp[ival] = 0;
        }
    }
    else {
        printf("get_two_oper_exp>> not implemented at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#endif
}

void get_linear_rsp_solution(const QInt size_pert,
                             const QReal *freq_sums,
                             const QInt num_freq_sums,
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

