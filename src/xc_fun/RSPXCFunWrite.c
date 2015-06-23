/* OpenRSP: open-ended library for response theory
   Copyright 2015 Radovan Bast,
                  Daniel H. Friese,
                  Bin Gao,
                  Dan J. Jonsson,
                  Magnus Ringholm,
                  Kenneth Ruud,
                  Andreas Thorvaldsen

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

   This file implements the function RSPXCFunWrite().

   2015-06-23, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_xc_fun.h"

/*% \brief writes the linked list of XC functionals
    \author Bin Gao
    \date 2015-06-23
    \param[RSPXCFun:struct]{in} xc_fun the linked list of XC functionals
    \param[FILE]{inout} fp_xc file pointer
    \return[QErrorCode:int] error information
*/
QErrorCode RSPXCFunWrite(RSPXCFun *xc_fun, FILE *fp_xc)
{
    QInt ixc;          /* incremental recorder over XC functionals */
    RSPXCFun *cur_xc;  /* current XC functional */
    QInt ipert;        /* incremental recorder over perturbations */
    /* walks to the last XC functional */
    ixc = 0;
    cur_xc = xc_fun;
    do {
        fprintf(fp_xc, "RSPXCFunWrite>> XC functional %d\n", ixc);
        fprintf(fp_xc,
                "RSPXCFunWrite>> number of perturbations that the XC functional depends on %d\n",
                cur_xc->num_pert);
        fprintf(fp_xc, "RSPXCFunWrite>> label           maximum-order\n");
        for (ipert=0; ipert<cur_xc->num_pert; ipert++) {
            fprintf(fp_xc,
                    "RSPXCFunWrite>>       %d                  %d\n",
                    cur_xc->pert_labels[ipert],
                    cur_xc->pert_max_orders[ipert]);
        }
        if (cur_xc->user_ctx!=NULL) {
            fprintf(fp_xc, "RSPXCFunWrite>> user-defined function context given\n");
        }
        ixc++;
        cur_xc = cur_xc->next_xc;
    } while (cur_xc!=NULL);
    return QSUCCESS;
}
