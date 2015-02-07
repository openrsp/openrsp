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

   This file implements the function OpenRSPSetPerturbations().

   2014-08-04, Bin Gao:
   * first version
*/

#include "openrsp.h"

/*@% \brief sets all perturbations involved in response theory calculations
     \author Bin Gao
     \date 2014-08-04
     \param[OneRSP:struct]{inout} open_rsp the context of response theory calculations
     \param[QInt:int]{in} num_pert number of all perturbations involved in calculations
     \param[QInt:int]{in} pert_labels labels of the perturbations
     \param[QInt:int]{in} pert_max_orders maximum allowed orders of all perturbations
     \param[QInt:int]{in} pert_sizes sizes of all perturbations up to their maximum orders,
         whose dimension is \sum{\var{pert_max_orders}}
     \param[QVoid:void]{in} user_ctx user-defined callback function context
     \param[GetPertComp:void]{in} get_pert_comp user specified function for
         getting components of a perturbation
     \param[GetPertRank:void]{in} get_pert_rank user specified function for
         getting rank of a perturbation
     \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPSetPerturbations(OpenRSP *open_rsp,
                                   const QInt num_pert,
                                   const QInt *pert_labels,
                                   const QInt *pert_max_orders,
                                   const QInt *pert_sizes,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   QVoid *user_ctx,
#endif
                                   const GetPertComp get_pert_comp,
                                   const GetPertRank get_pert_rank)
{
    QInt ipert,jpert,iorder;  /* incremental recorders */
    if (num_pert<1) {
        printf("OpenRSPSetPerturbations>> number of perturbations %"QINT_FMT"\n",
               num_pert);
        QErrorExit(FILE_AND_LINE, "invalid number of perturbations");
    }
    open_rsp->num_pert = num_pert;
    open_rsp->pert_labels = (QInt *)malloc(num_pert*sizeof(QInt));
    if (open_rsp->pert_labels==NULL) {
        printf("OpenRSPSetPerturbations>> number of perturbations %"QINT_FMT"\n",
               num_pert);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for pert_labels");
    }
    open_rsp->pert_max_orders = (QInt *)malloc(num_pert*sizeof(QInt));
    if (open_rsp->pert_max_orders==NULL) {
        printf("OpenRSPSetPerturbations>> number of perturbations %"QINT_FMT"\n",
               num_pert);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for pert_max_orders");
    }
    open_rsp->size_ptr = (QInt *)malloc((num_pert+1)*sizeof(QInt));
    if (open_rsp->size_ptr==NULL) {
        printf("OpenRSPSetPerturbations>> number of perturbations %"QINT_FMT"\n",
               num_pert);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for size_ptr");
    }
    open_rsp->size_ptr[0] = 0;
    for (ipert=0; ipert<num_pert; ipert++) {
        /* each element of \var{pert_labels} should be unique */
        for (jpert=0; jpert<ipert; jpert++) {
            if (pert_labels[jpert]==pert_labels[ipert]) {
                printf("OpenRSPSetPerturbations>> perturbation %"QINT_FMT" is %"QINT_FMT"\n",
                       jpert,
                       pert_labels[jpert]);
                printf("OpenRSPSetPerturbations>> perturbation %"QINT_FMT" is %"QINT_FMT"\n",
                       ipert,
                       pert_labels[ipert]);
                QErrorExit(FILE_AND_LINE, "repeated labels of perturbations not allowed");
            }
        }
        open_rsp->pert_labels[ipert] = pert_labels[ipert];
        if (pert_max_orders[ipert]<1) {
            printf("OpenRSPSetPerturbations>> order of %"QINT_FMT"-th perturbation (%"QINT_FMT") is %"QINT_FMT"\n",
                   ipert,
                   pert_labels[ipert],
                   pert_max_orders[ipert]);
            QErrorExit(FILE_AND_LINE, "only positive order allowed");
        }
        open_rsp->pert_max_orders[ipert] = pert_max_orders[ipert];
        /* \var{open_rsp->size_ptr[ipert]} indicates the start of sizes of
           \var{open_rsp->pert_labels[ipert]} */
        open_rsp->size_ptr[ipert+1] = open_rsp->size_ptr[ipert]+pert_max_orders[ipert];
    }
    /* the last element of \var{open_rsp->size_ptr[num_pert]} is the size of
       \var{open_rsp->pert_sizes} */
    open_rsp->pert_sizes = (QInt *)malloc(open_rsp->size_ptr[num_pert]*sizeof(QInt));
    if (open_rsp->pert_sizes==NULL) {
        printf("OpenRSPSetPerturbations>> size of pert_sizes %"QINT_FMT"\n",
               open_rsp->size_ptr[num_pert]);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for pert_sizes");
    }
    for (ipert=0,jpert=0; ipert<num_pert; ipert++) {
        for (iorder=1; iorder<=open_rsp->pert_max_orders[ipert]; iorder++,jpert++) {
            if (pert_sizes[jpert]<1) {
                printf("OpenRSPSetPerturbations>> size of %"QINT_FMT"-th perturbation (%"QINT_FMT", order %"QINT_FMT") is %"QINT_FMT"\n",
                       ipert,
                       pert_labels[ipert],
                       pert_max_orders[ipert],
                       pert_sizes[jpert]);
                QErrorExit(FILE_AND_LINE, "incorrect size");
            }
            open_rsp->pert_sizes[jpert] = pert_sizes[jpert];
        }
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    open_rsp->user_ctx = user_ctx;
#endif
    open_rsp->get_pert_comp = get_pert_comp;
    open_rsp->get_pert_rank = get_pert_rank;
    return QSUCCESS;
}
