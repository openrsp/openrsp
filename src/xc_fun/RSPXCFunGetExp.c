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

   This file implements the function RSPXCFunGetExp().

   2015-06-23, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_xc_fun.h"

/*% \brief gets expectation values of the linked list of XC functionals
    \author Bin Gao
    \date 2015-06-23
    \param[RSPXCFun:struct]{in} xc_fun the linked list of XC functionals
    \param[QInt:int]{in} len_tuple length of perturbation tuple on the XC functional
    \param[QInt:int]{in} pert_tuple perturbation tuple on the XC functional
    \param[QInt:int]{in} num_freq_configs the number of different frequency
        configurations to be considered for the perturbation tuple
    \param[QInt:int]{in} len_dmat_tuple the number of different perturbation
        tuples of the AO based density matrices passed
    \param[QInt:int]{in} idx_dmat_tuple indices of the density matrix
        perturbation tuples passed (canonically ordered)
    \param[QInt:int]{in} num_dmat number of collected AO based density matrices for
        the passed density matrix perturbation tuples and all frequency configurations
    \param[QcMat:struct]{in} dens_mat the collected AO based density matrices
    \param[QInt:int]{in} num_exp number of expectation values
    \param[QReal:real]{out} val_exp the expectation values
    \return[QErrorCode:int] error information
*/
QErrorCode RSPXCFunGetExp(RSPXCFun *xc_fun,
                          const QInt len_tuple,
                          const QInt *pert_tuple,
                          const QInt num_freq_configs,
                          const QInt len_dmat_tuple,
                          const QInt *idx_dmat_tuple,
                          const QInt num_dmat,
                          QcMat *dens_mat[],
                          const QInt num_exp,
                          QReal *val_exp)
{
    RSPXCFun *cur_xc;  /* current XC functional */
    /* walks to the last XC functional */
    cur_xc = xc_fun;
    do {
/*FIXME: checks perturbations if resulting zero integrals*/
        cur_xc->get_xc_fun_exp(len_tuple,
                               pert_tuple,
                               num_freq_configs,
                               len_dmat_tuple,
                               idx_dmat_tuple,
                               num_dmat,
                               dens_mat,
#if defined(OPENRSP_C_USER_CONTEXT)
                               cur_xc->user_ctx,
#endif
                               num_exp,
                               val_exp);
        cur_xc = cur_xc->next_xc;
    } while (cur_xc!=NULL);
    return QSUCCESS;
}
