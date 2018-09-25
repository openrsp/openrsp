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

   This file implements the function RSPTwoOperGetExp().

   2014-08-06, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_two_oper.h"

/*% \brief gets expectation values of the linked list of two-electron operators
    \author Bin Gao
    \date 2014-08-06
    \param[RSPTwoOper:struct]{in} two_oper the linked list of two-electron operators
    \param[QInt:int]{in} len_tuple length of perturbation tuple on the two-electron operator
    \param[QInt:int]{in} pert_tuple perturbation tuple on the two-electron operator
    \param[QInt:int]{in} len_dmat_tuple length of different perturbation tuples
        of the left-hand-side (LHS) and right-hand-side (RHS) AO based density
        matrices passed
    \param[QInt:int]{in} num_LHS_dmat number of LHS AO based density matrices
        passed for each LHS density matrix perturbation tuple
    \param[QcMat:struct]{in} LHS_dens_mat the LHS AO based density matrices
    \param[QInt:int]{in} num_RHS_dmat number of RHS AO based density matrices
        passed for each RHS density matrix perturbation tuple
    \param[QcMat:struct]{in} RHS_dens_mat the RHS AO based density matrices
    \param[QInt:int]{in} num_exp number of expectation values
    \param[QReal:real]{out} val_exp the expectation values
    \return[QErrorCode:int] error information
*/
QErrorCode RSPTwoOperGetExp(RSPTwoOper *two_oper,
                            const QInt len_tuple,
                            const QInt *pert_tuple,
                            const QInt len_dmat_tuple,
                            const QInt *num_LHS_dmat,
                            QcMat *LHS_dens_mat[],
                            const QInt *num_RHS_dmat,
                            QcMat *RHS_dens_mat[],
                            const QInt num_exp,
                            QReal *val_exp)
{
    RSPTwoOper *cur_oper;  /* current operator */
    /* walks to the last operator */
    cur_oper = two_oper;
    do {
/*FIXME: checks perturbations if resulting zero integrals*/
        cur_oper->get_two_oper_exp(len_tuple,
                                   pert_tuple,
                                   len_dmat_tuple,
                                   num_LHS_dmat,
                                   LHS_dens_mat,
                                   num_RHS_dmat,
                                   RHS_dens_mat,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   cur_oper->user_ctx,
#endif
                                   num_exp,
                                   val_exp);
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}
