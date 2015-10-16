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

#include "OpenRSPAODensCallback.h"

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
#include "ao_dens_ground_state_hf/num_atomic_orbitals.h"
#include "ao_dens_ground_state_hf/overlap_integrals.h"
#if defined(OPENRSP_C_USER_CONTEXT)
    char *overlap_context;
#endif
    QInt idx_block_row[1] = {IDX_BLOCK_ROW};
    QInt idx_block_col[1] = {IDX_BLOCK_COL};
    QcDataType data_type[1] = {QREALMAT};
    QBool assembled;
    QInt imat;
    QErrorCode ierr;
#if defined(OPENRSP_C_USER_CONTEXT)
    overlap_context = (char *)user_ctx;
    if (strcmp(overlap_context, "OVERLAP")!=0) {
        printf("get_overlap_mat>> incorrect operator at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#endif
    /* overlap integrals */
    if (oper_num_pert==0 && bra_num_pert==0 && ket_num_pert==0) {
        /* checks if the matrix is assembled or not */
        ierr = QcMatIsAssembled(val_int[0], &assembled);
        if (ierr!=QSUCCESS) {
            printf("get_overlap_mat>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatIsAssembled");
            exit(ierr);
        }
        if (assembled==QFALSE) {
            ierr = QcMatBlockCreate(val_int[0], 1);
            if (ierr!=QSUCCESS) {
                printf("get_overlap_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatBlockCreate");
                exit(ierr);
            }
            ierr = QcMatSetSymType(val_int[0], QSYMMAT);
            if (ierr!=QSUCCESS) {
                printf("get_overlap_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatSetSymType");
                exit(ierr);
            }
            ierr = QcMatSetDataType(val_int[0],
                                    1,
                                    idx_block_row,
                                    idx_block_col,
                                    data_type);
            if (ierr!=QSUCCESS) {
                printf("get_overlap_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatSetDataType");
                exit(ierr);
            }
            ierr = QcMatSetDimMat(val_int[0], NUM_AO, NUM_AO);
            if (ierr!=QSUCCESS) {
                printf("get_overlap_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatSetDimMat");
                exit(ierr);
            }
            ierr = QcMatAssemble(val_int[0]);
            if (ierr!=QSUCCESS) {
                printf("get_overlap_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatAssemble");
                exit(ierr);
            }
        }
        ierr = QcMatSetValues(val_int[0],
                              IDX_BLOCK_ROW,
                              IDX_BLOCK_COL,
                              IDX_FIRST_ROW,
                              NUM_AO,
                              IDX_FIRST_COL,
                              NUM_AO,
                              values_overlap,
                              NULL);
        if (ierr!=QSUCCESS) {
            printf("get_overlap_mat>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatSetValues");
            exit(ierr);
        }
    }
    else if (oper_num_pert==1 && bra_num_pert==0 && ket_num_pert==0) {
        if (oper_pert_labels[0]==PERT_GEOMETRIC) {
            printf("get_overlap_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        /* zero integrals */
        else if (oper_pert_labels[0]==PERT_DIPOLE) {
            for (imat=0; imat<num_int; imat++) {
                /* checks if the matrix is assembled or not */
                ierr = QcMatIsAssembled(val_int[imat], &assembled);
                if (ierr!=QSUCCESS) {
                    printf("get_overlap_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatIsAssembled");
                    exit(ierr);
                }
                if (assembled==QFALSE) {
                    ierr = QcMatBlockCreate(val_int[imat], 1);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatBlockCreate");
                        exit(ierr);
                    }
                    ierr = QcMatSetSymType(val_int[imat], QSYMMAT);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetSymType");
                        exit(ierr);
                    }
                    ierr = QcMatSetDataType(val_int[imat],
                                            1,
                                            idx_block_row,
                                            idx_block_col,
                                            data_type);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetDataType");
                        exit(ierr);
                    }
                    ierr = QcMatSetDimMat(val_int[imat], NUM_AO, NUM_AO);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetDimMat");
                        exit(ierr);
                    }
                    ierr = QcMatAssemble(val_int[imat]);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatAssemble");
                        exit(ierr);
                    }
                }
                ierr = QcMatZeroEntries(val_int[imat]);
                if (ierr!=QSUCCESS) {
                    printf("get_overlap_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatZeroEntries");
                    exit(ierr);
                }
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
        /* zero integrals */
        else if (bra_pert_labels[0]==PERT_DIPOLE) {
            for (imat=0; imat<num_int; imat++) {
                /* checks if the matrix is assembled or not */
                ierr = QcMatIsAssembled(val_int[imat], &assembled);
                if (ierr!=QSUCCESS) {
                    printf("get_overlap_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatIsAssembled");
                    exit(ierr);
                }
                if (assembled==QFALSE) {
                    ierr = QcMatBlockCreate(val_int[imat], 1);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatBlockCreate");
                        exit(ierr);
                    }
                    ierr = QcMatSetSymType(val_int[imat], QSYMMAT);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetSymType");
                        exit(ierr);
                    }
                    ierr = QcMatSetDataType(val_int[imat],
                                            1,
                                            idx_block_row,
                                            idx_block_col,
                                            data_type);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetDataType");
                        exit(ierr);
                    }
                    ierr = QcMatSetDimMat(val_int[imat], NUM_AO, NUM_AO);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetDimMat");
                        exit(ierr);
                    }
                    ierr = QcMatAssemble(val_int[imat]);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatAssemble");
                        exit(ierr);
                    }
                }
                ierr = QcMatZeroEntries(val_int[imat]);
                if (ierr!=QSUCCESS) {
                    printf("get_overlap_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatZeroEntries");
                    exit(ierr);
                }
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
        /* zero integrals */
        else if (ket_pert_labels[0]==PERT_DIPOLE) {
            for (imat=0; imat<num_int; imat++) {
                /* checks if the matrix is assembled or not */
                ierr = QcMatIsAssembled(val_int[imat], &assembled);
                if (ierr!=QSUCCESS) {
                    printf("get_overlap_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatIsAssembled");
                    exit(ierr);
                }
                if (assembled==QFALSE) {
                    ierr = QcMatBlockCreate(val_int[imat], 1);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatBlockCreate");
                        exit(ierr);
                    }
                    ierr = QcMatSetSymType(val_int[imat], QSYMMAT);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetSymType");
                        exit(ierr);
                    }
                    ierr = QcMatSetDataType(val_int[imat],
                                            1,
                                            idx_block_row,
                                            idx_block_col,
                                            data_type);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetDataType");
                        exit(ierr);
                    }
                    ierr = QcMatSetDimMat(val_int[imat], NUM_AO, NUM_AO);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetDimMat");
                        exit(ierr);
                    }
                    ierr = QcMatAssemble(val_int[imat]);
                    if (ierr!=QSUCCESS) {
                        printf("get_overlap_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatAssemble");
                        exit(ierr);
                    }
                }
                ierr = QcMatZeroEntries(val_int[imat]);
                if (ierr!=QSUCCESS) {
                    printf("get_overlap_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatZeroEntries");
                    exit(ierr);
                }
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
#include "ao_dens_ground_state_hf/num_atomic_orbitals.h"
#include "ao_dens_ground_state_hf/overlap_integrals.h"
#if defined(OPENRSP_C_USER_CONTEXT)
    char *overlap_context;
#endif
    QInt idx_block_row[1] = {IDX_BLOCK_ROW};
    QInt idx_block_col[1] = {IDX_BLOCK_COL};
    QcDataType data_type[1] = {QREALMAT};
    QcMat val_int[1];
    QInt idens;
    QErrorCode ierr;
#if defined(OPENRSP_C_USER_CONTEXT)
    overlap_context = (char *)user_ctx;
    if (strcmp(overlap_context, "OVERLAP")!=0) {
        printf("get_overlap_mat>> incorrect operator at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#endif
    /* overlap integrals */
    if (oper_num_pert==0 && bra_num_pert==0 && ket_num_pert==0) {
        ierr = QcMatCreate(&val_int[0]);
        if (ierr!=QSUCCESS) {
            printf("get_overlap_exp>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatCreate");
            exit(ierr);
        }
        ierr = QcMatBlockCreate(&val_int[0], 1);
        if (ierr!=QSUCCESS) {
            printf("get_overlap_exp>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatBlockCreate");
            exit(ierr);
        }
        ierr = QcMatSetSymType(&val_int[0], QSYMMAT);
        if (ierr!=QSUCCESS) {
            printf("get_overlap_exp>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatSetSymType");
            exit(ierr);
        }
        ierr = QcMatSetDataType(&val_int[0],
                                1,
                                idx_block_row,
                                idx_block_col,
                                data_type);
        if (ierr!=QSUCCESS) {
            printf("get_overlap_exp>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatSetDataType");
            exit(ierr);
        }
        ierr = QcMatSetDimMat(&val_int[0], NUM_AO, NUM_AO);
        if (ierr!=QSUCCESS) {
            printf("get_overlap_exp>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatSetDimMat");
            exit(ierr);
        }
        ierr = QcMatAssemble(&val_int[0]);
        if (ierr!=QSUCCESS) {
            printf("get_overlap_exp>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatAssemble");
            exit(ierr);
        }
        ierr = QcMatSetValues(&val_int[0],
                              IDX_BLOCK_ROW,
                              IDX_BLOCK_COL,
                              IDX_FIRST_ROW,
                              NUM_AO,
                              IDX_FIRST_COL,
                              NUM_AO,
                              values_overlap,
                              NULL);
        if (ierr!=QSUCCESS) {
            printf("get_overlap_exp>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatSetValues");
            exit(ierr);
        }
        for (idens=0; idens<num_dens; idens++) {
            ierr = QcMatGetMatProdTrace(&val_int[0],
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
        ierr = QcMatDestroy(&val_int[0]);
        if (ierr!=QSUCCESS) {
            printf("get_overlap_exp>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatDestroy");
            exit(ierr);
        }
    }
    else if (oper_num_pert==1 && bra_num_pert==0 && ket_num_pert==0) {
        if (oper_pert_labels[0]==PERT_GEOMETRIC) {
            printf("get_overlap_exp>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        /* zero integrals */
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
        /* zero integrals */
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
        /* zero integrals */
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
        /* zero integrals */
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
#include "ao_dens_ground_state_hf/num_atomic_orbitals.h"
#include "ao_dens_alpha_hf/dipole_length_integrals.h"
#if defined(OPENRSP_C_USER_CONTEXT)
    char *oneham_context = "ONEHAM";
    char *ext_field_context = "EXT_FIELD";
    char *one_oper_context;
#endif
    QInt idx_block_row[1] = {IDX_BLOCK_ROW};
    QInt idx_block_col[1] = {IDX_BLOCK_COL};
    QcDataType data_type[1] = {QREALMAT};
    QBool assembled;
    QInt imat;
    QErrorCode ierr;
#if defined(OPENRSP_C_USER_CONTEXT)
    one_oper_context = (char *)user_ctx;
    if (strcmp(one_oper_context, oneham_context)==0) {
        /* electric fields (zero integrals) */
        if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
            for (imat=0; imat<num_int; imat++) {
                /* checks if the matrix is assembled or not */
                ierr = QcMatIsAssembled(val_int[imat], &assembled);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatIsAssembled");
                    exit(ierr);
                }
                if (assembled==QFALSE) {
                    ierr = QcMatBlockCreate(val_int[imat], 1);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatBlockCreate");
                        exit(ierr);
                    }
                    ierr = QcMatSetSymType(val_int[imat], QSYMMAT);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetSymType");
                        exit(ierr);
                    }
                    ierr = QcMatSetDataType(val_int[imat],
                                            1,
                                            idx_block_row,
                                            idx_block_col,
                                            data_type);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetDataType");
                        exit(ierr);
                    }
                    ierr = QcMatSetDimMat(val_int[imat], NUM_AO, NUM_AO);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetDimMat");
                        exit(ierr);
                    }
                    ierr = QcMatAssemble(val_int[imat]);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatAssemble");
                        exit(ierr);
                    }
                }
                ierr = QcMatZeroEntries(val_int[imat]);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatZeroEntries");
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
    else if (strcmp(one_oper_context, ext_field_context)==0) {
        /* electric fields */
        if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
            /* checks if the matrix is assembled or not */
            for (imat=0; imat<num_int; imat++) {
                ierr = QcMatIsAssembled(val_int[imat], &assembled);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatIsAssembled");
                    exit(ierr);
                }
                if (assembled==QFALSE) {
                    ierr = QcMatBlockCreate(val_int[imat], 1);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatBlockCreate");
                        exit(ierr);
                    }
                    ierr = QcMatSetSymType(val_int[imat], QSYMMAT);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetSymType");
                        exit(ierr);
                    }
                    ierr = QcMatSetDataType(val_int[imat],
                                            1,
                                            idx_block_row,
                                            idx_block_col,
                                            data_type);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetDataType");
                        exit(ierr);
                    }
                    ierr = QcMatSetDimMat(val_int[imat], NUM_AO, NUM_AO);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetDimMat");
                        exit(ierr);
                    }
                    ierr = QcMatAssemble(val_int[imat]);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatAssemble");
                        exit(ierr);
                    }
                }
            }
            /* dipole length integrals */
            if (oper_pert_orders[0]==1) {
                for (imat=0; imat<3; imat++) {
                    ierr = QcMatSetValues(val_int[imat],
                                          IDX_BLOCK_ROW,
                                          IDX_BLOCK_COL,
                                          IDX_FIRST_ROW,
                                          NUM_AO,
                                          IDX_FIRST_COL,
                                          NUM_AO,
                                          &values_diplen[imat*SIZE_AO_MAT],
                                          NULL);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetValues");
                        exit(ierr);
                    }
                }
            }
            /* zero integrals */
            else {
                for (imat=0; imat<3; imat++) {
                    ierr = QcMatZeroEntries(val_int[imat]);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_mat>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatZeroEntries");
                        exit(ierr);
                    }
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
    /* electric fields */
    if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
        /* checks if the matrix is assembled or not */
        for (imat=0; imat<num_int; imat++) {
            ierr = QcMatIsAssembled(val_int[imat], &assembled);
            if (ierr!=QSUCCESS) {
                printf("get_one_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatIsAssembled");
                exit(ierr);
            }
            if (assembled==QFALSE) {
                ierr = QcMatBlockCreate(val_int[imat], 1);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatBlockCreate");
                    exit(ierr);
                }
                ierr = QcMatSetSymType(val_int[imat], QSYMMAT);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatSetSymType");
                    exit(ierr);
                }
                ierr = QcMatSetDataType(val_int[imat],
                                        num_blocks,
                                        idx_block_row,
                                        idx_block_col,
                                        data_type);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatSetDataType");
                    exit(ierr);
                }
                ierr = QcMatSetDimMat(val_int[imat], NUM_AO, NUM_AO);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatSetDimMat");
                    exit(ierr);
                }
                ierr = QcMatAssemble(val_int[imat]);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatAssemble");
                    exit(ierr);
                }
            }
        }
        /* dipole length integrals */
        if (oper_pert_orders[0]==1) {
            for (imat=0; imat<3; imat++) {
                ierr = QcMatSetValues(val_int[imat],
                                      IDX_BLOCK_ROW,
                                      IDX_BLOCK_COL,
                                      IDX_FIRST_ROW,
                                      NUM_AO,
                                      IDX_FIRST_COL,
                                      NUM_AO,
                                      values_diplen[imat*SIZE_AO_MAT],
                                      NULL);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatSetValues");
                    exit(ierr);
                }
            }
        }
        /* zero integrals */
        else {
            for (imat=0; imat<3; imat++) {
                ierr = QcMatZeroEntries(val_int[imat]);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatZeroEntries");
                    exit(ierr);
                }
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
#include "ao_dens_ground_state_hf/num_atomic_orbitals.h"
#include "ao_dens_alpha_hf/dipole_length_integrals.h"
#if defined(OPENRSP_C_USER_CONTEXT)
    char *oneham_context = "ONEHAM";
    char *ext_field_context = "EXT_FIELD";
    char *one_oper_context;
#endif
    QInt idx_block_row[1] = {IDX_BLOCK_ROW};
    QInt idx_block_col[1] = {IDX_BLOCK_COL};
    QcDataType data_type[1] = {QREALMAT};
    QcMat val_int[1];
    QInt offset_exp;
    QInt imat;
    QInt idens;
    QErrorCode ierr;
#if defined(OPENRSP_C_USER_CONTEXT)
    one_oper_context = (char *)user_ctx;
    if (strcmp(one_oper_context, oneham_context)==0) {
        /* electric fields (zero integrals) */
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
    else if (strcmp(one_oper_context, ext_field_context)==0) {
        /* electric fields */
        if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
            /* dipole length integrals */
            if (oper_pert_orders[0]==1) {
                offset_exp = 0;
                for (imat=0; imat<3; imat++) {
                    ierr = QcMatCreate(&val_int[0]);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_exp>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatCreate");
                        exit(ierr);
                    }
                    ierr = QcMatBlockCreate(&val_int[0], 1);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_exp>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatBlockCreate");
                        exit(ierr);
                    }
                    ierr = QcMatSetSymType(&val_int[0], QSYMMAT);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_exp>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetSymType");
                        exit(ierr);
                    }
                    ierr = QcMatSetDataType(&val_int[0],
                                            1,
                                            idx_block_row,
                                            idx_block_col,
                                            data_type);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_exp>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetDataType");
                        exit(ierr);
                    }
                    ierr = QcMatSetDimMat(&val_int[0], NUM_AO, NUM_AO);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_exp>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetDimMat");
                        exit(ierr);
                    }
                    ierr = QcMatAssemble(&val_int[0]);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_exp>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatAssemble");
                        exit(ierr);
                    }
                    ierr = QcMatSetValues(&val_int[0],
                                          IDX_BLOCK_ROW,
                                          IDX_BLOCK_COL,
                                          IDX_FIRST_ROW,
                                          NUM_AO,
                                          IDX_FIRST_COL,
                                          NUM_AO,
                                          &values_diplen[imat*SIZE_AO_MAT],
                                          NULL);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_exp>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatSetValues");
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
                                   "calling QcMatGetMatProdTrace");
                            exit(ierr);
                        }
                    }
                    ierr = QcMatDestroy(&val_int[0]);
                    if (ierr!=QSUCCESS) {
                        printf("get_one_oper_exp>> error happened at %s: %s\n",
                               FILE_AND_LINE,
                               "calling QcMatDestroy");
                        exit(ierr);
                    }
                    offset_exp += 2*num_dmat;
                }
            }
            /* zero integrals */
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
    /* electric fields */
    if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
        /* dipole length integrals */
        if (oper_pert_orders[0]==1) {
            offset_exp = 0;
            for (imat=0; imat<3; imat++) {
                ierr = QcMatCreate(&val_int[0]);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_exp>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatCreate");
                    exit(ierr);
                }
                ierr = QcMatBlockCreate(&val_int[0], 1);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_exp>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatBlockCreate");
                    exit(ierr);
                }
                ierr = QcMatSetSymType(&val_int[0], QSYMMAT);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_exp>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatSetSymType");
                    exit(ierr);
                }
                ierr = QcMatSetDataType(&val_int[0],
                                        1,
                                        idx_block_row,
                                        idx_block_col,
                                        data_type);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_exp>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatSetDataType");
                    exit(ierr);
                }
                ierr = QcMatSetDimMat(&val_int[0], NUM_AO, NUM_AO);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_exp>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatSetDimMat");
                    exit(ierr);
                }
                ierr = QcMatAssemble(&val_int[0]);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_exp>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatAssemble");
                    exit(ierr);
                }
                ierr = QcMatSetValues(&val_int[0],
                                      IDX_BLOCK_ROW,
                                      IDX_BLOCK_COL,
                                      IDX_FIRST_ROW,
                                      NUM_AO,
                                      IDX_FIRST_COL,
                                      NUM_AO,
                                      values_diplen[imat*SIZE_AO_MAT],
                                      NULL);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_exp>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatSetValues");
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
                               "calling QcMatGetMatProdTrace");
                        exit(ierr);
                    }
                }
                ierr = QcMatDestroy(val_int[0]);
                if (ierr!=QSUCCESS) {
                    printf("get_one_oper_exp>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatDestroy");
                    exit(ierr);
                }
                offset_exp += 2*num_dmat;
            }
        }
        /* zero integrals */
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
#include "ao_dens_ground_state_hf/num_atomic_orbitals.h"
#include "ao_dens_alpha_hf/two_electron_integrals.h"
#if defined(OPENRSP_C_USER_CONTEXT)
    char *two_oper_context;
#endif
    QInt idx_block_row[1] = {IDX_BLOCK_ROW};
    QInt idx_block_col[1] = {IDX_BLOCK_COL};
    QcDataType data_type[1] = {QREALMAT};
    QBool assembled;
    QInt imat;
    QErrorCode ierr;
    static QInt id_gmat = -1;
#if defined(OPENRSP_C_USER_CONTEXT)
    two_oper_context = (char *)user_ctx;
    if (strcmp(two_oper_context, "NONLAO")!=0) {
        printf("get_two_oper_mat>> not implemented at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
    else {
        /* checks if the matrix is assembled or not */
        for (imat=0; imat<num_int; imat++) {
            ierr = QcMatIsAssembled(val_int[imat], &assembled);
            if (ierr!=QSUCCESS) {
                printf("get_two_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatIsAssembled");
                exit(ierr);
            }
            if (assembled==QFALSE) {
                ierr = QcMatBlockCreate(val_int[imat], 1);
                if (ierr!=QSUCCESS) {
                    printf("get_two_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatBlockCreate");
                    exit(ierr);
                }
                ierr = QcMatSetSymType(val_int[imat], QSYMMAT);
                if (ierr!=QSUCCESS) {
                    printf("get_two_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatSetSymType");
                    exit(ierr);
                }
                ierr = QcMatSetDataType(val_int[imat],
                                        1,
                                        idx_block_row,
                                        idx_block_col,
                                        data_type);
                if (ierr!=QSUCCESS) {
                    printf("get_two_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatSetDataType");
                    exit(ierr);
                }
                ierr = QcMatSetDimMat(val_int[imat], NUM_AO, NUM_AO);
                if (ierr!=QSUCCESS) {
                    printf("get_two_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatSetDimMat");
                    exit(ierr);
                }
                ierr = QcMatAssemble(val_int[imat]);
                if (ierr!=QSUCCESS) {
                    printf("get_two_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatAssemble");
                    exit(ierr);
                }
            }
        }
        /* unperturbed two-electron integrals */
        if (oper_num_pert==0) {
            id_gmat++;
            if (id_gmat>5) {
                printf("get_two_oper_mat>> not implemented at %s\n",
                       FILE_AND_LINE);
                exit(QFAILURE);
            }
            ierr = QcMatSetValues(val_int[0],
                                  IDX_BLOCK_ROW,
                                  IDX_BLOCK_COL,
                                  IDX_FIRST_ROW,
                                  NUM_AO,
                                  IDX_FIRST_COL,
                                  NUM_AO,
                                  &values_gmat[id_gmat*SIZE_AO_MAT],
                                  NULL);
            if (ierr!=QSUCCESS) {
                printf("get_two_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatSetValues");
                exit(ierr);
            }
        }
        /* electric fields (zero integrals) */
        else if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
            for (imat=0; imat<num_int; imat++) {
                ierr = QcMatZeroEntries(val_int[imat]);
                if (ierr!=QSUCCESS) {
                    printf("get_two_oper_mat>> error happened at %s: %s\n",
                           FILE_AND_LINE,
                           "calling QcMatZeroEntries");
                    exit(ierr);
                }
            }
        }
        else {
            printf("get_two_oper_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
    }
#else
    /* checks if the matrix is assembled or not */
    for (imat=0; imat<num_int; imat++) {
        ierr = QcMatIsAssembled(val_int[imat], &assembled);
        if (ierr!=QSUCCESS) {
            printf("get_two_oper_mat>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatIsAssembled");
            exit(ierr);
        }
        if (assembled==QFALSE) {
            ierr = QcMatBlockCreate(val_int[imat], 1);
            if (ierr!=QSUCCESS) {
                printf("get_two_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatBlockCreate");
                exit(ierr);
            }
            ierr = QcMatSetSymType(val_int[imat], QSYMMAT);
            if (ierr!=QSUCCESS) {
                printf("get_two_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatSetSymType");
                exit(ierr);
            }
            ierr = QcMatSetDataType(val_int[imat],
                                    num_blocks,
                                    idx_block_row,
                                    idx_block_col,
                                    data_type);
            if (ierr!=QSUCCESS) {
                printf("get_two_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatSetDataType");
                exit(ierr);
            }
            ierr = QcMatSetDimMat(val_int[imat], NUM_AO, NUM_AO);
            if (ierr!=QSUCCESS) {
                printf("get_two_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatSetDimMat");
                exit(ierr);
            }
            ierr = QcMatAssemble(val_int[imat]);
            if (ierr!=QSUCCESS) {
                printf("get_two_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatAssemble");
                exit(ierr);
            }
        }
    }
    /* unperturbed two-electron integrals */
    if (oper_num_pert==0) {
        id_gmat++;
        if (id_gmat>5) {
            printf("get_two_oper_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        ierr = QcMatSetValues(val_int[0],
                              IDX_BLOCK_ROW,
                              IDX_BLOCK_COL,
                              IDX_FIRST_ROW,
                              NUM_AO,
                              IDX_FIRST_COL,
                              NUM_AO,
                              values_gmat[id_gmat*SIZE_AO_MAT],
                              NULL);
        if (ierr!=QSUCCESS) {
            printf("get_two_oper_mat>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatSetValues");
            exit(ierr);
        }
    }
    /* electric fields (zero integrals) */
    else if (oper_num_pert==1 && oper_pert_labels[0]==PERT_DIPOLE) {
        for (imat=0; imat<num_int; imat++) {
            ierr = QcMatZeroEntries(val_int[imat]);
            if (ierr!=QSUCCESS) {
                printf("get_two_oper_mat>> error happened at %s: %s\n",
                       FILE_AND_LINE,
                       "calling QcMatZeroEntries");
                exit(ierr);
            }
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
    if (strcmp(two_oper_context, "NONLAO")!=0) {
        printf("get_two_oper_exp>> not implemented at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
    else {
        /* electric fields */
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
    /* electric fields */
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

void get_nuc_contrib(const QInt nuc_num_pert,
                     const QcPertInt *nuc_pert_labels,
                     const QInt *nuc_pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT) 
                     void *user_ctx,
#endif
                     const QInt size_pert,
                     QReal *val_nuc)
{
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
#include "ao_dens_ground_state_hf/num_atomic_orbitals.h"
#include "ao_dens_alpha_hf/response_parameters.h"
    QInt idx_block_row[1] = {IDX_BLOCK_ROW};
    QInt idx_block_col[1] = {IDX_BLOCK_COL};
    QcDataType data_type[1] = {QREALMAT};
    QBool assembled;
    QErrorCode ierr;
    static QInt id_rsp_param = -1;
    id_rsp_param++;
    if (id_rsp_param>2) {
        printf("get_linear_rsp_solution>> not implemented at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
    /* checks if the matrix is assembled or not */
    ierr = QcMatIsAssembled(rsp_param[0], &assembled);
    if (ierr!=QSUCCESS) {
        printf("get_linear_rsp_solution>> error happened at %s: %s\n",
               FILE_AND_LINE,
               "calling QcMatIsAssembled");
        exit(ierr);
    }
    if (assembled==QFALSE) {
        ierr = QcMatBlockCreate(rsp_param[0], 1);
        if (ierr!=QSUCCESS) {
            printf("get_linear_rsp_solution>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatBlockCreate");
            exit(ierr);
        }
        ierr = QcMatSetDataType(rsp_param[0],
                                1,
                                idx_block_row,
                                idx_block_col,
                                data_type);
        if (ierr!=QSUCCESS) {
            printf("get_linear_rsp_solution>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatSetDataType");
            exit(ierr);
        }
        ierr = QcMatSetDimMat(rsp_param[0], NUM_AO, NUM_AO);
        if (ierr!=QSUCCESS) {
            printf("get_linear_rsp_solution>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatSetDimMat");
            exit(ierr);
        }
        ierr = QcMatAssemble(rsp_param[0]);
        if (ierr!=QSUCCESS) {
            printf("get_linear_rsp_solution>> error happened at %s: %s\n",
                   FILE_AND_LINE,
                   "calling QcMatAssemble");
            exit(ierr);
        }
    }
    ierr = QcMatSetValues(rsp_param[0],
                          IDX_BLOCK_ROW,
                          IDX_BLOCK_COL,
                          IDX_FIRST_ROW,
                          NUM_AO,
                          IDX_FIRST_COL,
                          NUM_AO,
                          &alpha_rsp_param[id_rsp_param*SIZE_AO_MAT],
                          NULL);
    if (ierr!=QSUCCESS) {
        printf("get_linear_rsp_solution>> error happened at %s: %s\n",
               FILE_AND_LINE,
               "calling QcMatSetValues");
        exit(ierr);
    }
}

