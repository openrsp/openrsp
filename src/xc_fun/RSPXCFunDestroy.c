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

   This file implements the function RSPTwoOperDestroy().

   2014-08-06, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_two_oper.h"

/*% \brief destroys the linked list of two-electron operators, should be called at the end
    \author Bin Gao
    \date 2014-08-06
    \param[RSPTwoOper:struct]{inout} two_oper the linked list of two-electron operators
    \return[QErrorCode:int] error information
*/
QErrorCode RSPTwoOperDestroy(RSPTwoOper **two_oper)
{
    RSPTwoOper *cur_oper;   /* current operator */
    RSPTwoOper *next_oper;  /* next operator */
    /* walks to the last operator */
    cur_oper = *two_oper;
    while (cur_oper!=NULL) {
        cur_oper->num_pert = 0;
        free(cur_oper->pert_labels);
        cur_oper->pert_labels = NULL;
        free(cur_oper->pert_max_orders);
        cur_oper->pert_max_orders = NULL;
        cur_oper->user_ctx = NULL;
        cur_oper->get_two_oper_mat = NULL;
        cur_oper->get_two_oper_exp = NULL;
        next_oper = cur_oper->next_oper;
        free(cur_oper);
        cur_oper = NULL;
        cur_oper = next_oper;
    }
    return QSUCCESS;
}
