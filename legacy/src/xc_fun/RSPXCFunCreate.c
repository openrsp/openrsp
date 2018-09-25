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

   This file implements the function RSPXCFunCreate().

   2015-06-23, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_xc_fun.h"

/*% \brief creates the linked list of XC functionals, should be called at first
    \author Bin Gao
    \date 2015-06-23
    \param[RSPXCFun:struct]{inout} xc_fun the linked list of XC functionals
    \param[QInt:int]{in} num_pert number of different perturbation labels that can
        act as perturbations on the XC functional
    \param[QInt:int]{in} pert_labels all the different perturbation labels
    \param[QInt:int]{in} pert_max_orders maximum allowed order of each perturbation (label)
    \param[QVoid:void]{in} user_ctx user-defined callback function context
    \param[GetXCFunMat:void]{in} get_xc_fun_mat user specified function for
        getting integral matrices
    \param[GetXCFunExp:void]{in} get_xc_fun_exp user specified function for
        getting expectation values
    \return[QErrorCode:int] error information
*/
QErrorCode RSPXCFunCreate(RSPXCFun **xc_fun,
                          const QInt num_pert,
                          const QInt *pert_labels,
                          const QInt *pert_max_orders,
                          QVoid *user_ctx,
                          const GetXCFunMat get_xc_fun_mat,
                          const GetXCFunExp get_xc_fun_exp)
{
    RSPXCFun *new_xc;  /* new XC functional */
    QInt ipert,jpert;  /* incremental recorder over perturbations */
    new_xc = (RSPXCFun *)malloc(sizeof(RSPXCFun));
    if (new_xc==NULL) {
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for new_xc");
    }
    if (num_pert>0) {
        new_xc->num_pert = num_pert;
    }
    else {
        printf("RSPXCFunCreate>> number of perturbations %d\n", num_pert);
        QErrorExit(FILE_AND_LINE, "invalid number of perturbations");
    }
    new_xc->pert_labels = (QInt *)malloc(num_pert*sizeof(QInt));
    if (new_xc->pert_labels==NULL) {
        printf("RSPXCFunCreate>> number of perturbations %d\n", num_pert);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for pert_labels");
    }
    new_xc->pert_max_orders = (QInt *)malloc(num_pert*sizeof(QInt));
    if (new_xc->pert_max_orders==NULL) {
        printf("RSPXCFunCreate>> number of perturbations %d\n", num_pert);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for pert_max_orders");
    }
    for (ipert=0; ipert<num_pert; ipert++) {
        /* each element of \var{pert_labels} should be unique */
        for (jpert=0; jpert<ipert; jpert++) {
            if (pert_labels[jpert]==pert_labels[ipert]) {
                printf("RSPXCFunCreate>> perturbation %d is %d\n",
                       jpert,
                       pert_labels[jpert]);
                printf("RSPXCFunCreate>> perturbation %d is %d\n",
                       ipert,
                       pert_labels[ipert]);
                QErrorExit(FILE_AND_LINE, "same perturbation not allowed");
            }
        }
        new_xc->pert_labels[ipert] = pert_labels[ipert];
        if (pert_max_orders[ipert]<1) {
            printf("RSPXCFunCreate>> order of %d-th perturbation (%d) is %d\n",
                   ipert,
                   pert_labels[ipert],
                   pert_max_orders[ipert]);
            QErrorExit(FILE_AND_LINE, "only positive order allowed");
        }
        new_xc->pert_max_orders[ipert] = pert_max_orders[ipert];
    }
    new_xc->user_ctx = user_ctx;
    new_xc->get_xc_fun_mat = get_xc_fun_mat;
    new_xc->get_xc_fun_exp = get_xc_fun_exp;
    new_xc->next_xc = NULL;
    *xc_fun = new_xc;
    return QSUCCESS;
}
