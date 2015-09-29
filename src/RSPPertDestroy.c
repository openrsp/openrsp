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

#include "RSPPerturbation.h"

/*% \brief destroys the context of all perturbations involved in calculations
    \author Bin Gao
    \date 2015-06-28
    \param[RSPPert:struct]{inout} rsp_pert context of all perturbations involved in calculations
    \return[QErrorCode:int] error information
*/
QErrorCode RSPPertDestroy(RSPPert *rsp_pert)
{
    rsp_pert->num_pert = 0;
    free(rsp_pert->pert_labels);
    rsp_pert->pert_labels = NULL;
    free(rsp_pert->pert_max_orders);
    rsp_pert->pert_max_orders = NULL;
    free(rsp_pert->ptr_ncomp);
    rsp_pert->ptr_ncomp = NULL;
    free(rsp_pert->pert_num_comps);
    rsp_pert->pert_num_comps = NULL;
#if defined(OPENRSP_C_USER_CONTEXT)
    rsp_pert->user_ctx = NULL;
#endif
    rsp_pert->get_pert_concatenation = NULL;
    return QSUCCESS;
}

