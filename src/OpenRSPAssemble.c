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

   This file implements the function OpenRSPAssemble().

   2014-07-30, Bin Gao:
   * first version
*/

#include "openrsp.h"

/*@% \brief assembles the context of response theory calculations and checks
         its validity, should be called before any function OpenRSPGet...(),
         otherwise the results might be incorrect
     \author Bin Gao
     \date 2014-07-30
     \param[OneRSP:struct]{inout} open_rsp the context of response theory calculations
     \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPAssemble(OpenRSP *open_rsp)
{
    QErrorCode ierr;  /* error information */
    open_rsp->assembled = QFALSE;
/*FIXME: to implement ierr = xxAssemble(open_rsp->elec_eom); */
    /* assembles the context of response equation solver */
    if (open_rsp->rsp_solver!=NULL) {
        ierr = RSPSolverAssemble(open_rsp->rsp_solver);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPSolverAssemble");
    }
    else {
        QErrorExit(FILE_AND_LINE, "response equation solver should be set via OpenRSPSetSolver()");
    }
#if defined(OPENRSP_PERTURBATION_FREE)
    /* assembles the context of perturbations */
    if (open_rsp->num_pert<1 ||
        open_rsp->pert_labels==NULL ||
        open_rsp->pert_max_orders==NULL ||
        open_rsp->size_ptr==NULL ||
        open_rsp->pert_sizes==NULL ||
        open_rsp->get_pert_comp==NULL ||
        open_rsp->get_pert_rank==NULL) {
        QErrorExit(FILE_AND_LINE, "perturbations should be set via OpenRSPSetPerturbations()");
    }
#endif
    /* assembles the overlap integrals */
    if (open_rsp->overlap!=NULL) {
        ierr = RSPOverlapAssemble(open_rsp->overlap);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOverlapAssemble");
    }
    /* assembles the linked list of one-electron operators */
    if (open_rsp->one_oper!=NULL) {
        ierr = RSPOneOperAssemble(open_rsp->one_oper);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPOneOperAssemble");
    }
    /* assembles the linked list of two-electron operators */
    if (open_rsp->two_oper!=NULL) {
        ierr = RSPTwoOperAssemble(open_rsp->two_oper);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPTwoOperAssemble");
    }
    /* assembles the linked list of XC functionals */
    if (open_rsp->xc_fun!=NULL) {
        ierr = RSPXCFunAssemble(open_rsp->xc_fun);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPXCFunAssemble");
    }
    /* assembles the nuclear Hamiltonian */
    if (open_rsp->nuc_hamilton!=NULL) {
        ierr = RSPNucHamiltonAssemble(open_rsp->nuc_hamilton);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPNucHamiltonAssemble");
    }
    open_rsp->assembled = QTRUE;
    return QSUCCESS;
}
