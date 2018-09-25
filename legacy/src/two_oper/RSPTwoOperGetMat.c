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

   This file implements the function RSPTwoOperGetMat().

   2014-08-06, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_two_oper.h"

/*% \brief gets integral matrices of the linked list of two-electron operators
    \author Bin Gao
    \date 2014-08-06
    \param[RSPTwoOper:struct]{in} two_oper the linked list of two-electron operators
    \param[QInt:int]{in} len_tuple length of perturbation tuple on the two-electron operator
    \param[QInt:int]{in} pert_tuple perturbation tuple on the two-electron operator
    \param[QInt:int]{in} num_dmat number of AO based density matrices
    \param[QcMat:struct]{in} dens_mat the AO based density matrices
    \param[QInt:int]{in} num_int number of the integral matrices
    \param[QcMat:struct]{inout} val_int the integral matrices
    \return[QErrorCode:int] error information
*/
QErrorCode RSPTwoOperGetMat(RSPTwoOper *two_oper,
                            const QInt len_tuple,
                            const QInt *pert_tuple,
                            const QInt num_dmat,
                            QcMat *dens_mat[],
                            const QInt num_int,
                            QcMat *val_int[])
{
    RSPTwoOper *cur_oper;  /* current operator */
    /* walks to the last operator */
    cur_oper = two_oper;
    do {
/*FIXME: checks perturbations if resulting zero integrals*/
        cur_oper->get_two_oper_mat(len_tuple,
                                   pert_tuple,
                                   num_dmat,
                                   dens_mat,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   cur_oper->user_ctx,
#endif
                                   num_int,
                                   val_int);
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}
