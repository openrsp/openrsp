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

   This file implements the function RSPTwoOperAdd().

   2014-08-06, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_two_oper.h"

/*% \brief adds a two-electron operator to the linked list
    \author Bin Gao
    \date 2014-08-06
    \param[RSPTwoOper:struct]{inout} two_oper the linked list of two-electron operators
    \param[QInt:int]{in} num_pert number of perturbations that the two-electron
        operator depends on
    \param[QInt:int]{in} pert_labels labels of the perturbations
    \param[QInt:int]{in} pert_max_orders maximum allowed orders of the perturbations
    \param[QVoid:void]{in} user_ctx user-defined callback function context
    \param[GetTwoOperMat:void]{in} get_two_oper_mat user specified function for
        getting integral matrices
    \param[GetTwoOperExp:void]{in} get_two_oper_exp user specified function for
        getting expectation values
    \return[QErrorCode:int] error information
*/
QErrorCode RSPTwoOperAdd(RSPTwoOper *two_oper,
                         const QInt num_pert,
                         const QInt *pert_labels,
                         const QInt *pert_max_orders,
                         QVoid *user_ctx,
                         const GetTwoOperMat get_two_oper_mat,
                         const GetTwoOperExp get_two_oper_exp)
{
    RSPTwoOper *new_oper;  /* new operator */
    RSPTwoOper *cur_oper;  /* current operator */
    QErrorCode ierr;       /* error information */
    /* creates the new operator */
    ierr = RSPTwoOperCreate(&new_oper,
                            num_pert,
                            pert_labels,
                            pert_max_orders,
                            user_ctx,
                            get_two_oper_mat,
                            get_two_oper_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPTwoOperCreate");
    /* walks to the last operator */
    cur_oper = two_oper;
    while (cur_oper->next_oper!=NULL) {
        cur_oper = cur_oper->next_oper;
    }
    /* inserts the new operator to the tail of the linked list */
    cur_oper->next_oper = new_oper;
    return QSUCCESS;
}
