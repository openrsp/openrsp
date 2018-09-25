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

   This file implements the function RSPOneOperDestroy().

   2014-07-30, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_one_oper.h"

/*% \brief destroys the linked list of one-electron operators, should be called at the end
    \author Bin Gao
    \date 2014-07-30
    \param[RSPOneOper:struct]{inout} one_oper the linked list of one-electron operators
    \return[QErrorCode:int] error information
*/
QErrorCode RSPOneOperDestroy(RSPOneOper **one_oper)
{
    RSPOneOper *cur_oper;   /* current operator */
    RSPOneOper *next_oper;  /* next operator */
    /* walks to the last operator */
    cur_oper = *one_oper;
    while (cur_oper!=NULL) {
        cur_oper->num_pert = 0;
        free(cur_oper->pert_labels);
        cur_oper->pert_labels = NULL;
        free(cur_oper->pert_max_orders);
        cur_oper->pert_max_orders = NULL;
#if defined(OPENRSP_C_USER_CONTEXT)
        cur_oper->user_ctx = NULL;
#endif
        cur_oper->get_one_oper_mat = NULL;
        cur_oper->get_one_oper_exp = NULL;
        next_oper = cur_oper->next_oper;
        free(cur_oper);
        cur_oper = NULL;
        cur_oper = next_oper;
    }
    return QSUCCESS;
}
