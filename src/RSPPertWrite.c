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

/*% \brief writes the context of all perturbations involved in calculations
    \author Bin Gao
    \date 2015-06-28
    \param[RSPPert:struct]{inout} rsp_pert context of all perturbations involved in calculations
    \param[FILE]{inout} fp_pert file pointer
    \return[QErrorCode:int] error information
*/
QErrorCode RSPPertWrite(const RSPPert *rsp_pert, FILE *fp_pert)
{
    QInt ipert,icomp;  /* incremental recorders */
    fprintf(fp_pert,
            "RSPPertWrite>> number of all perturbation lables %"QCPERTINT_FMT"\n",
            rsp_pert->num_pert);
    fprintf(fp_pert,
            "RSPPertWrite>> label           maximum-order    numbers-of-components\n");
    for (ipert=0; ipert<rsp_pert->num_pert; ipert++) {
        fprintf(fp_pert,
                "RSPPertWrite>>  %"QCPERTINT_FMT"               %"QINT_FMT"               ",
                rsp_pert->pert_labels[ipert],
                rsp_pert->pert_max_orders[ipert]);
        for (icomp=rsp_pert->ptr_ncomp[ipert]; icomp<rsp_pert->ptr_ncomp[ipert+1]; icomp++) {
            fprintf(fp_pert, " %"QINT_FMT"", rsp_pert->pert_num_comps[icomp]);
        }
        fprintf(fp_pert, "\n");
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    if (rsp_pert->user_ctx!=NULL) {
        fprintf(fp_pert, "RSPPertWrite>> user-defined function context given\n");
    }
#endif
    return QSUCCESS;
}

