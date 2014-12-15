/* OpenRSP: open-ended library for response theory
   Copyright 2014

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
    fprintf(fp_overlap, "RSPOverlapWrite>> perturbation    maximum-order\n");
    for (ipert=0; ipert<overlap->num_pert; ipert++) {
        fprintf(fp_overlap,
                "RSPOverlapWrite>>       %"QINT_FMT"                  %"QINT_FMT"\n",
                overlap->perturbations[ipert],
                overlap->pert_max_orders[ipert]);
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    if (overlap->user_ctx!=NULL) {
        fprintf(fp_overlap, "RSPOverlapWrite>> user-defined function context given\n");
    }
#endif
    return QSUCCESS;
}
