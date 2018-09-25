/* OpenRSP: open-ended library for response theory
   Copyright 2014
  This source code form is subject to the terms of the
  GNU Lesser General Public License, version 2.1.
  If a copy of the GNU LGPL v2.1 was not distributed with this
  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

   This file implements the callback function get_two_oper_mat.

   2014-07-31, Bin Gao:
   * first version
*/

/* QcMatrix library */
#include "qcmatrix.h"

#include "tests/ao_dens/openrsp_c_callback.h"
#include "tests/openrsp_c_perturbations.h"

QVoid get_two_oper_mat(const QInt num_pert,
                       const QInt *pert_labels,
                       const QInt *pert_orders,
                       const QInt num_var_dens,
                       QcMat *var_ao_dens[],
#if defined(OPENRSP_C_USER_CONTEXT)
                       QVoid *user_ctx,
#endif
                       const QInt num_int,
                       QcMat *val_int[])
{
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *two_oper_context;
#endif
#include "tests/ao_dens/openrsp_c_ao_dims.h"
#include "tests/ao_dens/openrsp_c_ao_gmat.h"
    QInt idx_block_row[1] = {IDX_BLOCK_ROW};
    QInt idx_block_col[1] = {IDX_BLOCK_COL};
    QcDataType data_type[1] = {QREALMAT};
    QBool assembled;
    QInt imat;
    QErrorCode ierr;
    static QInt id_gmat = -1;
#if defined(OPENRSP_C_USER_CONTEXT)
    two_oper_context = (QChar *)user_ctx;
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
        if (num_pert==0) {
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
        else if (num_pert==1 && pert_labels[0]==PERT_DIPOLE) {
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
    if (num_pert==0) {
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
    else if (num_pert==1 && pert_labels[0]==PERT_DIPOLE) {
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
