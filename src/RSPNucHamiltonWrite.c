/*
  OpenRSP: open-ended library for response theory
  Copyright 2015 Radovan Bast,
                 Daniel H. Friese,
                 Bin Gao,
                 Dan J. Jonsson,
                 Magnus Ringholm,
                 Kenneth Ruud,
                 Andreas Thorvaldsen

  OpenRSP is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.

  OpenRSP is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with OpenRSP. If not, see <http://www.gnu.org/licenses/>.

*/

#include "RSPNucHamilton.h"

/*% \brief writes the context of nuclear Hamiltonian
    \author Bin Gao
    \date 2015-02-12
    \param[RSPNucHamilton:struct]{in} nuc_hamilton the context of nuclear Hamiltonian
    \param[FILE]{inout} fp_nuc file pointer
    \return[QErrorCode:int] error information
*/
QErrorCode RSPNucHamiltonWrite(const RSPNucHamilton *nuc_hamilton, FILE *fp_nuc)
{
    QInt ipert;  /* incremental recorder over perturbations */
    fprintf(fp_nuc,
            "RSPNucHamiltonWrite>> number of perturbations that nuclear Hamiltonian depend on %"QINT_FMT"\n",
            nuc_hamilton->num_pert);
    fprintf(fp_nuc, "RSPNucHamiltonWrite>> label           maximum-order\n");
    for (ipert=0; ipert<nuc_hamilton->num_pert; ipert++) {
        fprintf(fp_nuc,
                "RSPNucHamiltonWrite>>       %"QINT_FMT"                  %"QINT_FMT"\n",
                nuc_hamilton->pert_labels[ipert],
                nuc_hamilton->pert_max_orders[ipert]);
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    if (nuc_hamilton->user_ctx!=NULL) {
        fprintf(fp_nuc, "RSPNucHamiltonWrite>> user-defined function context given\n");
    }
#endif
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
    fprintf(fp_nuc,
            "RSPNucHamiltonWrite>> number of atoms %"QINT_FMT"\n",
            nuc_hamilton->num_atoms);
    return QSUCCESS;
}

