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

   This file implements the function OpenRSPSetPerturbations().

   2015-06-29, Bin Gao:
   * first version
*/

#include "openrsp.h"

/*@% \brief sets all perturbations involved in response theory calculations
     \author Bin Gao
     \date 2015-06-29
     \param[OpenRSP:struct]{inout} open_rsp the context of response theory calculations
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
QErrorCode OpenRSPSetPerturbations(OpenRSP *open_rsp,
                                   const QInt num_pert,
                                   const QInt *pert_labels,
                                   const QInt *pert_max_orders,
                                   const QInt *pert_num_comps,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   QVoid *user_ctx,
#endif
                                   const GetPertCat get_pert_concatenation)
{
    QErrorCode ierr;  /* error information */
    /* creates the context of all perturbations involved in calculations */
    if (open_rsp->rsp_pert!=NULL) {
        ierr = RSPPertDestroy(open_rsp->rsp_pert);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertDestroy");
    }
    else {
        open_rsp->rsp_pert = (RSPPert *)malloc(sizeof(RSPPert));
        if (open_rsp->rsp_pert==NULL) {
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for rsp_pert");
        }
    }
    ierr = RSPPertCreate(open_rsp->rsp_pert,
                         num_pert,
                         pert_labels,
                         pert_max_orders,
                         pert_num_comps,
#if defined(OPENRSP_C_USER_CONTEXT)
                         user_ctx,
#endif
                         get_pert_concatenation);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertCreate");
    return QSUCCESS;
}
