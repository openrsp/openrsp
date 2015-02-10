/* OpenRSP: open-ended library for response theory
   Copyright 2014

   OpenRSP is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OpenRSP is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with OpenRSP. If not, see <http://www.gnu.org/licenses/>.

   This file implements the callback function get_overlap_mat.

   2014-07-31, Bin Gao:
   * first version
*/

/* QcMatrix library */
#include "qcmatrix.h"

#include "tests/ao_dens/openrsp_c_callback.h"
#include "tests/openrsp_c_perturbations.h"

QVoid get_overlap_mat(const QInt bra_num_pert,
                      const QInt *bra_pert_labels,
                      const QInt *bra_pert_orders,
                      const QInt ket_num_pert,
                      const QInt *ket_pert_labels,
                      const QInt *ket_pert_orders,
                      const QInt num_pert,
                      const QInt *pert_labels,
                      const QInt *pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                      QVoid *user_ctx,
#endif
                      const QInt num_int,
                      QcMat *val_int[])
{
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *overlap_context;
#endif
#include "tests/ao_dens/openrsp_c_ao_dims.h"
#include "tests/ao_dens/openrsp_c_ao_overlap.h"
    QInt idx_block_row[1] = {IDX_BLOCK_ROW};
    QInt idx_block_col[1] = {IDX_BLOCK_COL};
    QcDataType data_type[1] = {QREALMAT};
    QBool assembled;
    QInt imat;
    QErrorCode ierr;
#if defined(OPENRSP_C_USER_CONTEXT)
    overlap_context = (QChar *)user_ctx;
    if (strcmp(overlap_context, "OVERLAP")!=0) {
        printf("get_overlap_mat>> incorrect operator at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
#endif
    /* overlap integrals */
    if (num_pert==0 && bra_num_pert==0 && ket_num_pert==0) {
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
    else if (num_pert==1 && bra_num_pert==0 && ket_num_pert==0) {
        if (pert_labels[0]==PERT_GEOMETRIC) {
            printf("get_overlap_mat>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        /* zero integrals */
        else if (pert_labels[0]==PERT_DIPOLE) {
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
        else if (pert_labels[0]==PERT_MAGNETIC) {
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
    else if (num_pert==0 && bra_num_pert==1 && ket_num_pert==0) {
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
    else if (num_pert==0 && bra_num_pert==0 && ket_num_pert==1) {
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
