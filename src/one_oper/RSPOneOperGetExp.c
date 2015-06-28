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

   This file implements the function RSPOneOperGetExp().

   2014-08-03, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_one_oper.h"

/*% \brief gets expectation values of the linked list of one-electron operators
    \author Bin Gao
    \date 2014-08-03
    \param[RSPOneOper:struct]{in} one_oper the linked list of one-electron operators
    \param[QInt:int]{in} len_tuple length of perturbation tuple on the one-electron operator
    \param[QInt:int]{in} pert_tuple perturbation tuple on the one-electron operator
    \param[QInt:int]{in} num_dmat number of atomic orbital (AO) based density matrices
    \param[QcMat:struct]{in} dens_mat the AO based density matrices
    \param[QInt:int]{in} num_exp number of expectation values
    \param[QReal:real]{out} val_exp the expectation values
    \return[QErrorCode:int] error information
*/
QErrorCode RSPOneOperGetExp(RSPOneOper *one_oper,
                            const QInt len_tuple,
                            const QInt *pert_tuple,
                            const QInt num_dmat,
                            QcMat *dens_mat[],
                            const QInt num_exp,
                            QReal *val_exp)
{
    RSPOneOper *cur_oper;  /* current operator */
    /* walks to the last operator */
    cur_oper = one_oper;
    do {
/*FIXME: checks perturbations if resulting zero integrals*/
        cur_oper->get_one_oper_exp(len_tuple,
                                   pert_tuple,
                                   num_dmat,
                                   dens_mat,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   cur_oper->user_ctx,
#endif
                                   num_exp,
                                   val_exp);
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}
