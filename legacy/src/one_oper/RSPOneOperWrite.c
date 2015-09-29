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

   This file implements the function RSPOneOperWrite().

   2014-07-30, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_one_oper.h"

/*% \brief writes the linked list of one-electron operators
    \author Bin Gao
    \date 2014-07-30
    \param[RSPOneOper:struct]{in} one_oper the linked list of one-electron operators
    \param[FILE]{inout} fp_oper file pointer
    \return[QErrorCode:int] error information
*/
QErrorCode RSPOneOperWrite(RSPOneOper *one_oper, FILE *fp_oper)
{
    QInt ioper;            /* incremental recorder over opertors */
    RSPOneOper *cur_oper;  /* current operator */
    QInt ipert;            /* incremental recorder over perturbations */
    /* walks to the last operator */
    ioper = 0;
    cur_oper = one_oper;
    do {
        fprintf(fp_oper, "RSPOneOperWrite>> operator %"QINT_FMT"\n", ioper);
        fprintf(fp_oper,
                "RSPOneOperWrite>> number of perturbations that the operator depends on %"QINT_FMT"\n",
                cur_oper->num_pert);
        fprintf(fp_oper, "RSPOneOperWrite>> label           maximum-order\n");
        for (ipert=0; ipert<cur_oper->num_pert; ipert++) {
            fprintf(fp_oper,
                    "RSPOneOperWrite>>       %"QINT_FMT"                  %"QINT_FMT"\n",
                    cur_oper->pert_labels[ipert],
                    cur_oper->pert_max_orders[ipert]);
        }
#if defined(OPENRSP_C_USER_CONTEXT)
        if (cur_oper->user_ctx!=NULL) {
            fprintf(fp_oper, "RSPOneOperWrite>> user-defined function context given\n");
        }
#endif
        ioper++;
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}
