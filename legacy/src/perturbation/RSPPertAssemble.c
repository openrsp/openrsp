/* RSPPert: open-ended library for response theory
   Copyright 2015 Radovan Bast,
                  Daniel H. Friese,
                  Bin Gao,
                  Dan J. Jonsson,
                  Magnus Ringholm,
                  Kenneth Ruud,
                  Andreas Thorvaldsen

   RSPPert is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   RSPPert is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with RSPPert. If not, see <http://www.gnu.org/licenses/>.

   This file implements the function RSPPertAssemble().

   2015-06-28, Bin Gao:
   * first version
*/

#include "eom/openrsp_perturbation.h"

/*% \brief assembles the context of all perturbations involved in calculations
    \author Bin Gao
    \date 2015-06-28
    \param[RSPPert:struct]{inout} rsp_pert context of all perturbations involved in calculations
    \return[QErrorCode:int] error information
*/
QErrorCode RSPPertAssemble(RSPPert *rsp_pert)
{
    if (rsp_pert->num_pert<1 ||
        rsp_pert->pert_labels==NULL ||
        rsp_pert->pert_max_orders==NULL ||
        rsp_pert->ptr_ncomp==NULL ||
        rsp_pert->pert_num_comps==NULL ||
        rsp_pert->get_pert_concatenation==NULL) {
        QErrorExit(FILE_AND_LINE, "perturbations are not correctly set");
    }
    return QSUCCESS;
}
