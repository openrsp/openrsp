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

   This file implements the function RSPNucHamiltonGetContributions().

   2015-02-12, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_nuc_contrib.h"

/*% \brief gets the nuclear contributions
    \author Bin Gao
    \date 2015-02-12
    \param[RSPNucHamilton:struct]{in} nuc_hamilton the context of nuclear Hamiltonian
    \param[QInt:int]{in} len_tuple length of perturbation tuple on the nuclear Hamiltonian
    \param[QInt:int]{in} pert_tuple perturbation tuple on the nuclear Hamiltonian
    \param[QInt:int]{in} size_pert size of the perturbations on the nuclear Hamiltonian
    \param[QReal:real]{inout} val_nuc the nuclear contributions
    \return[QErrorCode:int] error information
*/
QErrorCode RSPNucHamiltonGetContributions(const RSPNucHamilton *nuc_hamilton,
                                          const QInt len_tuple,
                                          const QInt *pert_tuple,
                                          const QInt size_pert,
                                          QReal *val_nuc)
{
    nuc_hamilton->get_nuc_contrib(len_tuple,
                                  pert_tuple,
#if defined(OPENRSP_C_USER_CONTEXT)
                                  nuc_hamilton->user_ctx,
#endif
                                  size_pert,
                                  val_nuc);
    return QSUCCESS;
}
