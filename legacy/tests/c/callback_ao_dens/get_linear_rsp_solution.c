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

   This file implements the callback function get_linear_rsp_solution.

   2014-07-31, Bin Gao:
   * first version
*/

/* QcMatrix library */
#include "qcmatrix.h"

#include "tests/ao_dens/openrsp_c_callback.h"

QVoid get_linear_rsp_solution(const QInt size_pert,
                              const QInt num_freq_sums,
                              const QReal *freq_sums,
                              QcMat *RHS_mat[],
#if defined(OPENRSP_C_USER_CONTEXT)
                              QVoid *user_ctx,
#endif
                              QcMat *rsp_param[])
{
#include "tests/ao_dens/openrsp_c_ao_dims.h"
#include "tests/ao_dens/openrsp_c_ao_alpha_rsp_param.h"
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
