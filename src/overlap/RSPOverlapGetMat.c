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

   This file implements the function RSPOverlapGetMat().

   2014-08-05, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_overlap.h"

/*% \brief gets integral matrices of the overlap integrals
    \author Bin Gao
    \date 2014-08-05
    \param[RSPOverlap:struct]{in} overlap the overlap integrals
    \param[QInt:int]{in} bra_num_pert number of perturbations on the bra
    \param[QInt:int]{in} bra_pert_labels labels of the perturbations on the bra
    \param[QInt:int]{in} bra_pert_orders orders of the perturbations on the bra
    \param[QInt:int]{in} ket_num_pert number of perturbations on the ket
    \param[QInt:int]{in} ket_pert_labels labels of the perturbations on the ket
    \param[QInt:int]{in} ket_pert_orders orders of the perturbations on the ket
    \param[QInt:int]{in} num_pert number of perturbations
    \param[QInt:int]{in} pert_labels labels of the perturbations
    \param[QInt:int]{in} pert_orders orders of the perturbations
    \param[QInt:int]{in} num_int number of the integral matrices
    \param[QMat:struct]{inout} val_int the integral matrices
    \return[QErrorCode:int] error information
*/
QErrorCode RSPOverlapGetMat(const RSPOverlap *overlap,
                            const QInt bra_num_pert,
                            const QInt *bra_pert_labels,
                            const QInt *bra_pert_orders,
                            const QInt ket_num_pert,
                            const QInt *ket_pert_labels,
                            const QInt *ket_pert_orders,
                            const QInt num_pert,
                            const QInt *pert_labels,
                            const QInt *pert_orders,
                            const QInt num_int,
                            QMat *val_int[])
{
    overlap->get_overlap_mat(bra_num_pert,
                             bra_pert_labels,
                             bra_pert_orders,
                             ket_num_pert,
                             ket_pert_labels,
                             ket_pert_orders,
                             num_pert,
                             pert_labels,
                             pert_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                             overlap->user_ctx,
#endif
                             num_int,
                             val_int);
    return QSUCCESS;
}
