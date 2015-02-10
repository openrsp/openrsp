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

   This file implements the callback function get_two_oper_exp.

   2014-07-31, Bin Gao:
   * first version
*/

/* QcMatrix library */
#include "qcmatrix.h"

#include "tests/ao_dens/openrsp_c_callback.h"
#include "tests/openrsp_c_perturbations.h"

QVoid get_two_oper_exp(const QInt num_pert,
                       const QInt *pert_labels,
                       const QInt *pert_orders,
                       const QInt num_var_dens,
                       QcMat *var_ao_dens[],
                       const QInt num_contr_dens,
                       QcMat *contr_ao_dens[],
#if defined(OPENRSP_C_USER_CONTEXT)
                       QVoid *user_ctx,
#endif
                       const QInt num_exp,
                       QReal *val_exp)
{
#if defined(OPENRSP_C_USER_CONTEXT)
    QChar *two_oper_context;
#endif
    QInt ival;
    QErrorCode ierr;
#if defined(OPENRSP_C_USER_CONTEXT)
    two_oper_context = (QChar *)user_ctx;
    if (strcmp(two_oper_context, "NONLAO")!=0) {
        printf("get_two_oper_exp>> not implemented at %s\n",
               FILE_AND_LINE);
        exit(QFAILURE);
    }
    else {
        /* electric fields */
        if (num_pert==1 && pert_labels[0]==PERT_DIPOLE) {
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
    if (num_pert==1 && pert_labels[0]==PERT_DIPOLE) {
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
