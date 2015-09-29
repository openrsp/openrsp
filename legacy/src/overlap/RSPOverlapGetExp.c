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

   This file implements the function RSPOverlapGetExp().

   2014-08-05, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_overlap.h"

/*% \brief gets expectation values of the overlap integrals
    \author Bin Gao
    \date 2014-08-05
    \param[RSPOverlap:struct]{in} overlap the overlap integrals
    \param[QInt:int]{in} bra_len_tuple length of the perturbation tuple on the bra
    \param[QInt:int]{in} bra_pert_tuple perturbation tuple on the bra
    \param[QInt:int]{in} ket_len_tuple length of the perturbation tuple on the ket
    \param[QInt:int]{in} ket_pert_tuple perturbation tuple on the ket
    \param[QInt:int]{in} len_tuple length of perturbation tuple on the overlap integrals
    \param[QInt:int]{in} pert_tuple perturbation tuple on the overlap integrals
    \param[QInt:int]{in} num_dmat number of atomic orbital (AO) based density matrices
    \param[QcMat:struct]{in} dens_mat the AO based density matrices
    \param[QInt:int]{in} num_exp number of expectation values
    \param[QReal:real]{out} val_exp the expectation values
    \return[QErrorCode:int] error information
*/
QErrorCode RSPOverlapGetExp(const RSPOverlap *overlap,
                            const QInt bra_len_tuple,
                            const QInt *bra_pert_tuple,
                            const QInt ket_len_tuple,
                            const QInt *ket_pert_tuple,
                            const QInt len_tuple,
                            const QInt *pert_tuple,
                            const QInt num_dmat,
                            QcMat *dens_mat[],
                            const QInt num_exp,
                            QReal *val_exp)
{
/*FIXME: checks perturbations if resulting zero integrals*/
    overlap->get_overlap_exp(bra_len_tuple,
                             bra_pert_tuple,
                             ket_len_tuple,
                             ket_pert_tuple,
                             len_tuple,
                             pert_tuple,
                             num_dmat,
                             dens_mat,
#if defined(OPENRSP_C_USER_CONTEXT)
                             overlap->user_ctx,
#endif
                             num_exp,
                             val_exp);
    return QSUCCESS;
}
