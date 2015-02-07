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

   This file implements the function RSPTwoOperGetExp().

   2014-08-06, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_two_oper.h"

/*% \brief gets expectation values of the linked list of two-electron operators
    \author Bin Gao
    \date 2014-08-06
    \param[RSPTwoOper:struct]{in} two_oper the linked list of two-electron operators
    \param[QInt:int]{in} num_pert number of perturbations
    \param[QInt:int]{in} pert_labels labels of the perturbations
    \param[QInt:int]{in} pert_orders orders of the perturbations
    \param[QInt:int]{in} num_var_dens number of variable AO based density matrices
    \param[QMat:struct]{in} var_ao_dens the variable AO based density matrices (\math{D})
        for calculating \math{G(D)}
    \param[QInt:int]{in} num_contr_dens number of contracted AO based density matrices
    \param[QMat:struct]{in} contr_ao_dens the contracted AO based density matrices (\math{D})
        for calculating \math{\mathrm{tr}[GD]}
    \param[QInt:int]{in} num_exp number of expectation values
    \param[QReal:real]{out} val_exp the expectation values
    \return[QErrorCode:int] error information
*/
QErrorCode RSPTwoOperGetExp(RSPTwoOper *two_oper,
                            const QInt num_pert,
                            const QInt *pert_labels,
                            const QInt *pert_orders,
                            const QInt num_var_dens,
                            QMat *var_ao_dens[],
                            const QInt num_contr_dens,
                            QMat *contr_ao_dens[],
                            const QInt num_exp,
                            QReal *val_exp)
{
    RSPTwoOper *cur_oper;  /* current operator */
    /* walks to the last operator */
    cur_oper = two_oper;
    do {
        cur_oper->get_two_oper_exp(num_pert,
                                   pert_labels,
                                   pert_orders,
                                   num_var_dens,
                                   var_ao_dens,
                                   num_contr_dens,
                                   contr_ao_dens,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   cur_oper->user_ctx,
#endif
                                   num_exp,
                                   val_exp);
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}
