/* OpenRSP: open-ended library for response theory
   Copyright 2014
  This source code form is subject to the terms of the
  GNU Lesser General Public License, version 2.1.
  If a copy of the GNU LGPL v2.1 was not distributed with this
  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

   This file implements the callback function get_one_oper_exp.

   2014-07-31, Bin Gao:
   * first version
*/

/* QcMatrix library */
#include "qcmatrix.h"

#include "tests/dens_mat/openrsp_c_callback.h"
#include "tests/openrsp_c_perturbations.h"

QVoid get_one_oper_exp(const QInt len_tuple,
                       const QInt *pert_tuple,
                       const QInt num_dmat,
                       QcMat *dens_mat[],
#if defined(OPENRSP_C_USER_CONTEXT)
                       QVoid *user_ctx,
#endif
                       const QInt num_exp,
                       QReal *val_exp)
{
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *oneham_context = "ONEHAM";
    QChar *ext_field_context = "EXT_FIELD";
    QChar *one_oper_context;
#endif
#include "tests/dens_mat/openrsp_c_ao_dims.h"
#include "tests/dens_mat/openrsp_c_ao_diplen.h"
    QInt idx_block_row[1] = {IDX_BLOCK_ROW};
    QInt idx_block_col[1] = {IDX_BLOCK_COL};
    QcDataType data_type[1] = {QREALMAT};
    QcMat val_int[1];
    QInt offset_exp;
    QInt imat;
    QInt idens;
    QErrorCode ierr;
#if defined(OPENRSP_C_USER_CONTEXT)
    one_oper_context = (QChar *)user_ctx;
    if (strcmp(one_oper_context, oneham_context)==0) {
        /* electric fields (zero integrals) */
        if (len_tuple==1 && pert_tuple[0]==PERT_DIPOLE) {
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
        if (len_tuple==1 && pert_tuple[0]==PERT_DIPOLE) {
            /* dipole length integrals */
            if (pert_orders[0]==1) {
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
    if (len_tuple==1 && pert_tuple[0]==PERT_DIPOLE) {
        /* dipole length integrals */
        if (pert_orders[0]==1) {
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
