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

   This file implements the function OpenRSPAddOneOper().

   2014-07-30, Bin Gao:
   * first version
*/

#include "openrsp.h"

/*@% \brief adds a one-electron operator to the Hamiltonian
     \author Bin Gao
     \date 2014-07-30
     \param[OneRSP:struct]{inout} open_rsp the context of response theory calculations
     \param[QInt:int]{in} num_pert number of perturbations that the one-electron
         operator depends on
     \param[QInt:int]{in} perturbations perturbations that the one-electron operator
         depends on
     \param[QInt:int]{in} pert_max_orders maximum allowed orders of the perturbations
     \param[QVoid:void]{in} user_ctx user-defined callback function context
     \param[GetOneOperMat:void]{in} get_one_oper_mat user specified function for
         getting integral matrices
     \param[GetOneOperExp:void]{in} get_one_oper_exp user specified function for
         getting expectation values
     \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPAddOneOper(OpenRSP *open_rsp,
                             const QInt num_pert,
                             const QInt *perturbations,
                             const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             QVoid *user_ctx,
#endif
                             const GetOneOperMat get_one_oper_mat,
                             const GetOneOperExp get_one_oper_exp)
{
    QErrorCode ierr;  /* error information */
    /* creates the linked list of one-electron operators */
    if (open_rsp->one_oper==NULL) {
        ierr = RSPOneOperCreate(&open_rsp->one_oper,
                                num_pert,
                                perturbations,
                                pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                user_ctx,
#endif
                                get_one_oper_mat,
                                get_one_oper_exp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOneOperCreate");
    }
    /* adds the one-electron operator to the linked list */
    else {
        ierr = RSPOneOperAdd(open_rsp->one_oper,
                             num_pert,
                             perturbations,
                             pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             user_ctx,
#endif
                             get_one_oper_mat,
                             get_one_oper_exp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOneOperAdd");
    }
    return QSUCCESS;
}
