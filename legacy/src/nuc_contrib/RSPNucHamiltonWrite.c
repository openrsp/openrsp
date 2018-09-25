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

   This file implements the function RSPNucHamiltonWrite().

   2015-02-12, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_nuc_contrib.h"

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
