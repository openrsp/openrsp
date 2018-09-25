/* OpenRSP: open-ended library for response theory
   Copyright 2014
  This source code form is subject to the terms of the
  GNU Lesser General Public License, version 2.1.
  If a copy of the GNU LGPL v2.1 was not distributed with this
  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

   This file implements the callback function get_overlap_exp.

   2014-07-31, Bin Gao:
   * first version
*/

/* QcMatrix library */
#include "qcmatrix.h"

#include "tests/ao_dens/openrsp_c_callback.h"
#include "tests/openrsp_c_perturbations.h"

QVoid get_overlap_exp(const QInt bra_num_pert,
                      const QInt *bra_pert_labels,
                      const QInt *bra_pert_orders,
                      const QInt ket_num_pert,
                      const QInt *ket_pert_labels,
                      const QInt *ket_pert_orders,
                      const QInt num_pert,
                      const QInt *pert_labels,
                      const QInt *pert_orders,
                      const QInt num_dens,
                      QcMat *ao_dens[],
#if defined(OPENRSP_C_USER_CONTEXT)
                      QVoid *user_ctx,
#endif
                      const QInt num_exp,
                      QReal *val_exp)
{
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *overlap_context;
#endif
#include "tests/ao_dens/openrsp_c_ao_dims.h"
#include "tests/ao_dens/openrsp_c_ao_overlap.h"
    QInt idx_block_row[1] = {IDX_BLOCK_ROW};
    QInt idx_block_col[1] = {IDX_BLOCK_COL};
    QcDataType data_type[1] = {QREALMAT};
    QcMat val_int[1];
    QInt idens;
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
    else if (num_pert==1 && bra_num_pert==0 && ket_num_pert==0) {
        if (pert_labels[0]==PERT_GEOMETRIC) {
            printf("get_overlap_exp>> not implemented at %s\n",
                   FILE_AND_LINE);
            exit(QFAILURE);
        }
        /* zero integrals */
        else if (pert_labels[0]==PERT_DIPOLE) {
            for (idens=0; idens<2*num_exp; idens++) {
                val_exp[idens] = 0;
            }
        }
        else if (pert_labels[0]==PERT_MAGNETIC) {
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
    else if (num_pert==0 && bra_num_pert==1 && ket_num_pert==0) {
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
    else if (num_pert==0 && bra_num_pert==0 && ket_num_pert==1) {
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
    else if (num_pert==0 && bra_num_pert==1 && ket_num_pert==1) {
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
