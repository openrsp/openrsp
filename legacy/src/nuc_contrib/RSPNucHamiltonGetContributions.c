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
/*FIXME: checks perturbations if resulting zero values*/
    nuc_hamilton->get_nuc_contrib(len_tuple,
                                  pert_tuple,
#if defined(OPENRSP_C_USER_CONTEXT)
                                  nuc_hamilton->user_ctx,
#endif
                                  size_pert,
                                  val_nuc);
    return QSUCCESS;
}
