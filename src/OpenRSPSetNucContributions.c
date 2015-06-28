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

   This file implements the function OpenRSPSetNucContributions().

   2015-02-12, Bin Gao:
   * first version
*/

#include "openrsp.h"

/*% \brief sets the context of nuclear Hamiltonian
    \author Bin Gao
    \date 2015-02-12
    \param[OpenRSP:struct]{inout} open_rsp the context of response theory calculations
     \param[QInt:int]{in} num_pert number of different perturbation labels that can
         act as perturbations on the nuclear Hamiltonian
     \param[QInt:int]{in} pert_labels all the different perturbation labels
     \param[QInt:int]{in} pert_max_orders maximum allowed order of each perturbation (label)
    \param[QVoid:void]{in} user_ctx user-defined callback function context
    \param[GetNucContrib:void]{in} get_nuc_contrib user specified function for
        getting nuclear contributions
    \param[QInt:int]{in} num_atoms number of atoms
    \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPSetNucContributions(OpenRSP *open_rsp,
                                      const QInt num_pert,
                                      const QInt *pert_labels,
                                      const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                      QVoid *user_ctx,
#endif

                                      const GetNucContrib get_nuc_contrib,
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
                                      const QInt num_atoms)
{
    QErrorCode ierr;  /* error information */
    /* creates the context of nuclear Hamiltonian */
    if (open_rsp->nuc_hamilton!=NULL) {
        ierr = RSPNucHamiltonDestroy(open_rsp->nuc_hamilton);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPNucHamiltonDestroy");
    }
    else {
        open_rsp->nuc_hamilton = (RSPNucHamilton *)malloc(sizeof(RSPNucHamilton));
        if (open_rsp->nuc_hamilton==NULL) {
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for nuc_hamiltonian");
        }
    }
    ierr = RSPNucHamiltonCreate(open_rsp->nuc_hamilton,
                                num_pert,
                                pert_labels,
                                pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                                user_ctx,
#endif
                                get_nuc_contrib,
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
                                num_atoms);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPNucHamiltonCreate");
    return QSUCCESS;
}
