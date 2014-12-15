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

   This file implements the function RSPTwoOperGetMat().

   2014-08-06, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_two_oper.h"

/*% \brief gets integral matrices of the linked list of two-electron operators
    \author Bin Gao
    \date 2014-08-06
    \param[RSPTwoOper:struct]{in} two_oper the linked list of two-electron operators
    \param[QInt:int]{in} num_pert number of perturbations
    \param[QInt:int]{in} perturbations the perturbations
    \param[QInt:int]{in} pert_orders orders of the perturbations
    \param[QInt:int]{in} num_var_dens number of variable AO based density matrices
    \param[QMat:struct]{in} var_ao_dens the variable AO based density matrices (\math{D})
        for calculating \math{G(D)}
    \param[QInt:int]{in} num_int number of the integral matrices
    \param[QMat:struct]{inout} val_int the integral matrices
    \return[QErrorCode:int] error information
*/
QErrorCode RSPTwoOperGetMat(RSPTwoOper *two_oper,
                            const QInt num_pert,
                            const QInt *perturbations,
                            const QInt *pert_orders,
                            const QInt num_var_dens,
                            QMat *var_ao_dens[],
                            const QInt num_int,
                            QMat *val_int[])
{
    RSPTwoOper *cur_oper;  /* current operator */
    /* walks to the last operator */
    cur_oper = two_oper;
    do {
        cur_oper->get_two_oper_mat(num_pert,
                                   perturbations,
                                   pert_orders,
                                   num_var_dens,
                                   var_ao_dens,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   cur_oper->user_ctx,
#endif
                                   num_int,
                                   val_int);
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}
