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
     \param[OpenRSP:struct]{inout} open_rsp the context of response theory calculations
     \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPAssemble(OpenRSP *open_rsp)
{
    QErrorCode ierr;  /* error information */
    open_rsp->assembled = QFALSE;
/*FIXME: to implement ierr = xxAssemble(open_rsp->elec_eom); */
    /* assembles the context of perturbations */
    if (open_rsp->rsp_pert!=NULL) {
        ierr = RSPPertAssemble(open_rsp->rsp_pert);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPPertAssemble");
    }
#if defined(OPENRSP_PERTURBATION_FREE)
    else {
        QErrorExit(FILE_AND_LINE, "perturbations should be set via OpenRSPSetPerturbations()");
    }
#endif
    /* assembles the context of response equation solver */
    if (open_rsp->rsp_solver!=NULL) {
        ierr = RSPSolverAssemble(open_rsp->rsp_solver);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPSolverAssemble");
    }
    else {
        QErrorExit(FILE_AND_LINE, "linear response equation solver should be set via OpenRSPSetSolver()");
    }
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
