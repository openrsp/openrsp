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

/*% \brief sets all perturbations involved in response theory calculations
    \author Bin Gao
    \date 2015-06-28
    \param[RSPPert:struct]{inout} rsp_pert context of all perturbations involved in calculations
    \param[QInt:int]{in} num_pert number of all different perturbation labels involved
        in calculations
    \param[QInt:int]{in} pert_labels all different perturbation labels involved
    \param[QInt:int]{in} pert_max_orders maximum allowed order of each perturbation (label)
    \param[QInt:int]{in} pert_num_comps number of components of each perturbation (label)
        up to its maximum order, size is \sum{\var{pert_max_orders}}
    \param[QVoid:void]{in} user_ctx user-defined callback function context
    \param[GetPertCat:void]{in} get_pert_concatenation user specified function for
        getting the ranks of components of sub-perturbation tuples (with same
        perturbation label) for given components of the corresponding concatenated
        perturbation tuple
    \return[QErrorCode:int] error information
*/
QErrorCode RSPPertCreate(RSPPert *rsp_pert,
                         const QcPertInt num_pert,
                         const QcPertInt *pert_labels,
                         const QInt *pert_max_orders,
                         const QInt *pert_num_comps,
#if defined(OPENRSP_C_USER_CONTEXT)
                         QVoid *user_ctx,
#endif
                         const GetPertCat get_pert_concatenation)
{
    QInt ipert;   /* incremental recorders over perturbations */
    QInt jpert;
    QInt iorder;  /* incremental recorder over orders */
    if (num_pert<1) {
        printf("RSPPertCreate>> number of perturbation labels %"QCPERTINT_FMT"\n",
               num_pert);
        QErrorExit(FILE_AND_LINE, "invalid number of perturbation labels");
    }
    else if (num_pert>OPENRSP_PERT_LABEL_MAX) {
        printf("RSPPertCreate>> number of perturbation labels %"QCPERTINT_FMT"\n",
               num_pert);
        printf("RSPPertCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
               OPENRSP_PERT_LABEL_MAX);
        QErrorExit(FILE_AND_LINE, "too many perturbation labels");
    }
    rsp_pert->num_pert = num_pert;
    rsp_pert->pert_labels = (QcPertInt *)malloc(num_pert*sizeof(QcPertInt));
    if (rsp_pert->pert_labels==NULL) {
        printf("RSPPertCreate>> number of perturbation labels %"QCPERTINT_FMT"\n",
               num_pert);
        QErrorExit(FILE_AND_LINE, "allocates memory for perturbation labels");
    }
    rsp_pert->pert_max_orders = (QInt *)malloc(num_pert*sizeof(QInt));
    if (rsp_pert->pert_max_orders==NULL) {
        printf("RSPPertCreate>> number of perturbation labels %"QCPERTINT_FMT"\n",
               num_pert);
        QErrorExit(FILE_AND_LINE, "allocates memory for allowed maximal orders");
    }
    rsp_pert->ptr_ncomp = (QInt *)malloc((num_pert+1)*sizeof(QInt));
    if (rsp_pert->ptr_ncomp==NULL) {
        printf("RSPPertCreate>> number of perturbation labels %"QCPERTINT_FMT"\n",
               num_pert);
        QErrorExit(FILE_AND_LINE, "allocates memory for pointers to components");
    }
    rsp_pert->ptr_ncomp[0] = 0;
    for (ipert=0; ipert<num_pert; ipert++) {
        if (pert_labels[ipert]>OPENRSP_PERT_LABEL_MAX) {
            printf("RSPPertCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                   ipert,
                   pert_labels[ipert]);
            printf("RSPPertCreate>> maximal value for pert. labels %"QCPERTINT_FMT"\n",
                   OPENRSP_PERT_LABEL_MAX);
            QErrorExit(FILE_AND_LINE, "invalid perturbation label");
        }
        /* each element of <pert_labels> should be unique */
        for (jpert=0; jpert<ipert; jpert++) {
            if (pert_labels[jpert]==pert_labels[ipert]) {
                printf("RSPPertCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       jpert,
                       pert_labels[jpert]);
                printf("RSPPertCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ipert,
                       pert_labels[ipert]);
                QErrorExit(FILE_AND_LINE, "repeated perturbation labels not allowed");
            }
        }
        rsp_pert->pert_labels[ipert] = pert_labels[ipert];
        if (pert_max_orders[ipert]<1) {
            printf("RSPPertCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                   ipert,
                   pert_labels[ipert]);
            printf("RSPPertCreate>> allowed maximal order is %"QINT_FMT"\n",
                   pert_max_orders[ipert]);
            QErrorExit(FILE_AND_LINE, "only positive order allowed");
        }
        rsp_pert->pert_max_orders[ipert] = pert_max_orders[ipert];
        /* <c>rsp_pert->ptr_ncomp[ipert]</c> points to the number of components
           of <c>rsp_pert->pert_labels[ipert]</c> */
        rsp_pert->ptr_ncomp[ipert+1] = rsp_pert->ptr_ncomp[ipert]
                                     + pert_max_orders[ipert];
    }
    /* <c>rsp_pert->ptr_ncomp[num_pert]</c> equals to the size of
       <c>rsp_pert->pert_num_comps</c> */
    rsp_pert->pert_num_comps = (QInt *)malloc(rsp_pert->ptr_ncomp[num_pert]
                                              *sizeof(QInt));
    if (rsp_pert->pert_num_comps==NULL) {
        printf("RSPPertCreate>> size of numbers of components %"QINT_FMT"\n",
               rsp_pert->ptr_ncomp[num_pert]);
        QErrorExit(FILE_AND_LINE, "allocates memory for numbers of components");
    }
    for (ipert=0,jpert=0; ipert<num_pert; ipert++) {
        for (iorder=1; iorder<=rsp_pert->pert_max_orders[ipert]; iorder++,jpert++) {
            if (pert_num_comps[jpert]<1) {
                printf("RSPPertCreate>> %"QINT_FMT"-th pert. label %"QCPERTINT_FMT"\n",
                       ipert,
                       pert_labels[ipert]);
                printf("RSPPertCreate>> allowed maximal order is %"QINT_FMT"\n",
                       pert_max_orders[ipert]);
                printf("RSPPertCreate>> number of components is %"QINT_FMT"\n",
                       pert_num_comps[jpert]);
                QErrorExit(FILE_AND_LINE, "incorrect number of components");
            }
            rsp_pert->pert_num_comps[jpert] = pert_num_comps[jpert];
        }
    }
#if defined(OPENRSP_C_USER_CONTEXT)
    rsp_pert->user_ctx = user_ctx;
#endif
    rsp_pert->get_pert_concatenation = get_pert_concatenation;
    return QSUCCESS;
}
