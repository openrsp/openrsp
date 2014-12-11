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

   This file implements the function OpenRSPDestroy().

   2014-01-28, Bin Gao:
   * first version
*/

#include "openrsp.h"

/*@% \brief destroys the context of response theory calculations, should
         be called at the end
     \author Bin Gao
     \date 2014-01-28
     \param[OneRSP:struct]{inout} open_rsp the context of response theory calculations
     \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPDestroy(OpenRSP *open_rsp)
{
    QErrorCode ierr;  /* error information */
    open_rsp->assembled = QFALSE;
    if (open_rsp->elec_eom!=NULL) {
/*FIXME: to implement ierr = xxDestroy(open_rsp->elec_eom); */
        free(open_rsp->elec_eom);
        open_rsp->elec_eom = NULL;
    }
    /* destroys the context of response equation sovler */
    if (open_rsp->rsp_solver!=NULL) {
        ierr = RSPSolverDestroy(open_rsp->rsp_solver);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPSolverDestroy");
    }
#if defined(OPENRSP_PERTURBATION_FREE)
    /* destroys the context of perturbations */
    open_rsp->num_pert = 0;
    if (open_rsp->perturbations!=NULL) {
        free(open_rsp->perturbations);
        open_rsp->perturbations = NULL;
    }
    if (open_rsp->pert_max_orders!=NULL) {
        free(open_rsp->pert_max_orders);
        open_rsp->pert_max_orders = NULL;
    }
    if (open_rsp->size_ptr!=NULL) {
        free(open_rsp->size_ptr);
        open_rsp->size_ptr = NULL;
    }
    if (open_rsp->pert_sizes!=NULL) {
        free(open_rsp->pert_sizes);
        open_rsp->pert_sizes = NULL;
    }
    open_rsp->get_pert_comp = NULL;
    open_rsp->get_pert_rank = NULL;
#endif
    /* destroys the context of overlap integrals */
    if (open_rsp->overlap!=NULL) {
        ierr = RSPOverlapDestroy(open_rsp->overlap);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOverlapDestroy");
        free(open_rsp->overlap);
        open_rsp->overlap = NULL;
    }
    /* destroys the linked list of one-electron operators */
    if (open_rsp->one_oper!=NULL) {
        ierr = RSPOneOperDestroy(&open_rsp->one_oper);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOneOperDestroy");
    }
    /* destroys the linked list of two-electron operators */
    if (open_rsp->two_oper!=NULL) {
        ierr = RSPTwoOperDestroy(&open_rsp->two_oper);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPTwoOperDestroy");
    }
    /* destroys the linked list of exchange-correlation functionals */
    if (open_rsp->xc_fun!=NULL) {
    }
    /* destroys the context of (derivatives of) nuclear repulsion and nuclei-field interaction */
    if (open_rsp->nuc_contrib!=NULL) {
        free(open_rsp->nuc_contrib);
        open_rsp->nuc_contrib = NULL;
    }
    return QSUCCESS;
}
