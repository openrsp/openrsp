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

   This file implements the function RSPOverlapGetMat().

   2014-08-05, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_overlap.h"

/*% \brief gets integral matrices of the overlap integrals
    \author Bin Gao
    \date 2014-08-05
    \param[RSPOverlap:struct]{in} overlap the overlap integrals
    \param[QInt:int]{in} bra_len_tuple length of the perturbation tuple on the bra
    \param[QInt:int]{in} bra_pert_tuple perturbation tuple on the bra
    \param[QInt:int]{in} ket_len_tuple length of the perturbation tuple on the ket
    \param[QInt:int]{in} ket_pert_tuple perturbation tuple on the ket
    \param[QInt:int]{in} len_tuple length of perturbation tuple on the overlap integrals
    \param[QInt:int]{in} pert_tuple perturbation tuple on the overlap integrals
    \param[QInt:int]{in} num_int number of the integral matrices
    \param[QcMat:struct]{inout} val_int the integral matrices
    \return[QErrorCode:int] error information
*/
QErrorCode RSPOverlapGetMat(const RSPOverlap *overlap,
                            const QInt bra_len_tuple,
                            const QInt *bra_pert_tuple,
                            const QInt ket_len_tuple,
                            const QInt *ket_pert_tuple,
                            const QInt len_tuple,
                            const QInt *pert_tuple,
                            const QInt num_int,
                            QcMat *val_int[])
{
/*FIXME: checks perturbations if resulting zero integrals*/
    overlap->get_overlap_mat(bra_len_tuple,
                             bra_pert_tuple,
                             ket_len_tuple,
                             ket_pert_tuple,
                             len_tuple,
                             pert_tuple,
#if defined(OPENRSP_C_USER_CONTEXT)
                             overlap->user_ctx,
#endif
                             num_int,
                             val_int);
    return QSUCCESS;
}
