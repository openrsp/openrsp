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

   This file implements the function OpenRSPSetPDBS().

   2014-07-30, Bin Gao:
   * first version
*/

#include "openrsp.h"

/*@% \brief sets the context of perturbation dependent basis sets
     \author Bin Gao
     \date 2014-07-30
     \param[OpenRSP:struct]{inout} open_rsp the context of response theory calculations
     \param[QInt:int]{in} num_pert number of different perturbation labels that can
         act as perturbations on the basis sets
     \param[QInt:int]{in} pert_labels all the different perturbation labels
     \param[QInt:int]{in} pert_max_orders maximum allowed order of each perturbation (label)
     \param[QVoid:void]{in} user_ctx user-defined callback function context
     \param[GetOverlapMat:void]{in} get_overlap_mat user specified function for
         getting overlap integrals
     \param[GetOverlapExp:void]{in} get_overlap_exp user specified function for
         getting expectation values of overlap integrals
     \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPSetPDBS(OpenRSP *open_rsp,
                          const QInt num_pert,
                          const QInt *pert_labels,
                          const QInt *pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                          QVoid *user_ctx,
#endif
                          const GetOverlapMat get_overlap_mat,
                          const GetOverlapExp get_overlap_exp)
{
    QErrorCode ierr;  /* error information */
    /* creates the context of perturbation dependent basis sets */
    if (open_rsp->overlap!=NULL) {
        ierr = RSPOverlapDestroy(open_rsp->overlap);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOverlapDestroy");
    }
    else {
        open_rsp->overlap = (RSPOverlap *)malloc(sizeof(RSPOverlap));
        if (open_rsp->overlap==NULL) {
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for overlap");
        }
    }
    ierr = RSPOverlapCreate(open_rsp->overlap,
                            num_pert,
                            pert_labels,
                            pert_max_orders,
#if defined(OPENRSP_C_USER_CONTEXT)
                            user_ctx,
#endif
                            get_overlap_mat,
                            get_overlap_exp);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOverlapCreate");
    return QSUCCESS;
}
