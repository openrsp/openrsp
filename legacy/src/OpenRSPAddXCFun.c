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

   This file implements the function OpenRSPAddXCFun().

   2015-06-23, Bin Gao:
   * first version
*/

#include "openrsp.h"

/*@% \brief adds an XC functional  to the Hamiltonian
     \author Bin Gao
     \date 2015-06-23
     \param[OpenRSP:struct]{inout} open_rsp the context of response theory calculations
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
QErrorCode OpenRSPAddXCFun(OpenRSP *open_rsp,
                           const QInt num_pert,
                           const QInt *pert_labels,
                           const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                           QVoid *user_ctx,
#endif
                           const GetXCFunMat get_xc_fun_mat,
                           const GetXCFunExp get_xc_fun_exp)
{
    QErrorCode ierr;  /* error information */
    /* creates the linked list of XC functionals */
    if (open_rsp->xc_fun==NULL) {
        ierr = RSPXCFunCreate(&open_rsp->xc_fun,
                              num_pert,
                              pert_labels,
                              pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                              user_ctx,
#endif
                              get_xc_fun_mat,
                              get_xc_fun_exp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPXCFunCreate");
    }
    /* adds the XC functional to the linked list */
    else {
        ierr = RSPXCFunAdd(open_rsp->xc_fun,
                           num_pert,
                           pert_labels,
                           pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                           user_ctx,
#endif
                           get_xc_fun_mat,
                           get_xc_fun_exp);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPXCFunAdd");
    }
    return QSUCCESS;
}
