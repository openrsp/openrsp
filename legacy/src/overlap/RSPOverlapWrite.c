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

   This file implements the function RSPOverlapWrite().

   2014-08-05, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_overlap.h"

/*% \brief writes the overlap integrals
    \author Bin Gao
    \date 2014-08-05
    \param[RSPOverlap:struct]{in} overlap the overlap integrals
    \param[FILE]{inout} fp_overlap file pointer
    \return[QErrorCode:int] error information
*/
QErrorCode RSPOverlapWrite(const RSPOverlap *overlap, FILE *fp_overlap)
{
    QInt ipert;  /* incremental recorder over perturbations */
    fprintf(fp_overlap,
            "RSPOverlapWrite>> number of perturbations that overlap integrals depend on %"QINT_FMT"\n",
            overlap->num_pert);
    fprintf(fp_overlap, "RSPOverlapWrite>> label           maximum-order\n");
    for (ipert=0; ipert<overlap->num_pert; ipert++) {
        fprintf(fp_overlap,
                "RSPOverlapWrite>>       %"QINT_FMT"                  %"QINT_FMT"\n",
                overlap->pert_labels[ipert],
                overlap->pert_max_orders[ipert]);
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    if (overlap->user_ctx!=NULL) {
        fprintf(fp_overlap, "RSPOverlapWrite>> user-defined function context given\n");
    }
#endif
    return QSUCCESS;
}
