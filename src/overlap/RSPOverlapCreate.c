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

   This file implements the function RSPOverlapCreate().

   2014-08-05, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_overlap.h"

/*% \brief creates the overlap integrals, should be called at first
    \author Bin Gao
    \date 2014-08-05
    \param[RSPOverlap:struct]{inout} overlap the overlap integrals
    \param[QInt:int]{in} num_pert number of perturbations that the overlap integrals depend on
    \param[QInt:int]{in} pert_labels labels of the perturbations
    \param[QInt:int]{in} pert_max_orders maximum allowed orders of the perturbations
    \param[QVoid:void]{in} user_ctx user-defined callback function context
    \param[GetOverlapMat:void]{in} get_overlap_mat user specified function for
        getting integral matrices
    \param[GetOverlapExp:void]{in} get_overlap_exp user specified function for
        getting expectation values
    \return[QErrorCode:int] error information
*/
QErrorCode RSPOverlapCreate(RSPOverlap *overlap,
                            const QInt num_pert,
                            const QInt *pert_labels,
                            const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                            QVoid *user_ctx,
#endif
                            const GetOverlapMat get_overlap_mat,
                            const GetOverlapExp get_overlap_exp)
{
    QInt ipert,jpert;      /* incremental recorder over perturbations */
    if (num_pert>0) {
        overlap->num_pert = num_pert;
    }
    else {
        printf("RSPOverlapCreate>> number of perturbations %"QINT_FMT"\n", num_pert);
        QErrorExit(FILE_AND_LINE, "invalid number of perturbations");
    }
    overlap->pert_labels = (QInt *)malloc(num_pert*sizeof(QInt));
    if (overlap->pert_labels==NULL) {
        printf("RSPOverlapCreate>> number of perturbations %"QINT_FMT"\n", num_pert);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for pert_labels");
    }
    overlap->pert_max_orders = (QInt *)malloc(num_pert*sizeof(QInt));
    if (overlap->pert_max_orders==NULL) {
        printf("RSPOverlapCreate>> number of perturbations %"QINT_FMT"\n", num_pert);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for pert_max_orders");
    }
    for (ipert=0; ipert<num_pert; ipert++) {
        /* each element of \var{pert_labels} should be unique */
        for (jpert=0; jpert<ipert; jpert++) {
            if (pert_labels[jpert]==pert_labels[ipert]) {
                printf("RSPOverlapCreate>> perturbation %"QINT_FMT" is %"QINT_FMT"\n",
                       jpert,
                       pert_labels[jpert]);
                printf("RSPOverlapCreate>> perturbation %"QINT_FMT" is %"QINT_FMT"\n",
                       ipert,
                       pert_labels[ipert]);
                QErrorExit(FILE_AND_LINE, "same perturbation not allowed");
            }
        }
        overlap->pert_labels[ipert] = pert_labels[ipert];
        if (pert_max_orders[ipert]<1) {
            printf("RSPOverlapCreate>> order of %"QINT_FMT"-th perturbation (%"QINT_FMT") is %"QINT_FMT"\n",
                   ipert,
                   pert_labels[ipert],
                   pert_max_orders[ipert]);
            QErrorExit(FILE_AND_LINE, "only positive order allowed");
        }
        overlap->pert_max_orders[ipert] = pert_max_orders[ipert];
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    overlap->user_ctx = user_ctx;
#endif
    overlap->get_overlap_mat = get_overlap_mat;
    overlap->get_overlap_exp = get_overlap_exp;
    return QSUCCESS;
}
