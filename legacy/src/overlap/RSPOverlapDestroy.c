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

   This file implements the function RSPOverlapDestroy().

   2014-08-05, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_overlap.h"

/*% \brief destroys the overlap integrals, should be called at the end
    \author Bin Gao
    \date 2014-08-05
    \param[RSPOverlap:struct]{inout} overlap the overlap integrals
    \return[QErrorCode:int] error information
*/
QErrorCode RSPOverlapDestroy(RSPOverlap *overlap)
{
    overlap->num_pert = 0;
    free(overlap->pert_labels);
    overlap->pert_labels = NULL;
    free(overlap->pert_max_orders);
    overlap->pert_max_orders = NULL;
#if defined(OPENRSP_C_USER_CONTEXT)
    overlap->user_ctx = NULL;
#endif
    overlap->get_overlap_mat = NULL;
    overlap->get_overlap_exp = NULL;
    return QSUCCESS;
}
