/* OpenRSP: open-ended library for response theory
   Copyright 2015 Radovan Bast,
                  Daniel H. Friese,
                  Bin Gao,
                  Dan J. Jonsson,
                  Magnus Ringholm,
                  Kenneth Ruud,
                  Andreas Thorvaldsen
  This source code form is subject to the terms of the
  GNU Lesser General Public License, version 2.1.
  If a copy of the GNU LGPL v2.1 was not distributed with this
  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

   This file implements the function RSPTwoOperWrite().

   2014-08-06, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_two_oper.h"

/*% \brief writes the linked list of two-electron operators
    \author Bin Gao
    \date 2014-08-06
    \param[RSPTwoOper:struct]{in} two_oper the linked list of two-electron operators
    \param[FILE]{inout} fp_oper file pointer
    \return[QErrorCode:int] error information
*/
QErrorCode RSPTwoOperWrite(RSPTwoOper *two_oper, FILE *fp_oper)
{
    QInt ioper;            /* incremental recorder over opertors */
    RSPTwoOper *cur_oper;  /* current operator */
    QInt ipert;            /* incremental recorder over perturbations */
    /* walks to the last operator */
    ioper = 0;
    cur_oper = two_oper;
    do {
        fprintf(fp_oper, "RSPTwoOperWrite>> operator %"QINT_FMT"\n", ioper);
        fprintf(fp_oper,
                "RSPTwoOperWrite>> number of perturbations that the operator depends on %"QINT_FMT"\n",
                cur_oper->num_pert);
        fprintf(fp_oper, "RSPTwoOperWrite>> label           maximum-order\n");
        for (ipert=0; ipert<cur_oper->num_pert; ipert++) {
            fprintf(fp_oper,
                    "RSPTwoOperWrite>>       %"QINT_FMT"                  %"QINT_FMT"\n",
                    cur_oper->pert_labels[ipert],
                    cur_oper->pert_max_orders[ipert]);
        }
#if defined(OPENRSP_C_USER_CONTEXT)
        if (cur_oper->user_ctx!=NULL) {
            fprintf(fp_oper, "RSPTwoOperWrite>> user-defined function context given\n");
        }
#endif
        ioper++;
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}
