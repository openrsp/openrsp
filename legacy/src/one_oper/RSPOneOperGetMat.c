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

   This file implements the function RSPOneOperGetMat().

   2014-07-31, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_one_oper.h"

/*% \brief gets integral matrices of the linked list of one-electron operators
    \author Bin Gao
    \date 2014-07-31
    \param[RSPOneOper:struct]{in} one_oper the linked list of one-electron operators
    \param[QInt:int]{in} len_tuple length of perturbation tuple on the one-electron operator
    \param[QInt:int]{in} pert_tuple perturbation tuple on the one-electron operator
    \param[QInt:int]{in} num_int number of the integral matrices
    \param[QcMat:struct]{inout} val_int the integral matrices
    \return[QErrorCode:int] error information
*/
QErrorCode RSPOneOperGetMat(RSPOneOper *one_oper,
                            const QInt len_tuple,
                            const QInt *pert_tuple,
                            const QInt num_int,
                            QcMat *val_int[])
{
    RSPOneOper *cur_oper;  /* current operator */
    /* walks to the last operator */
    cur_oper = one_oper;
    do {
/*FIXME: checks perturbations if resulting zero integrals*/
        cur_oper->get_one_oper_mat(len_tuple,
                                   pert_tuple,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   cur_oper->user_ctx,
#endif
                                   num_int,
                                   val_int);
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}
