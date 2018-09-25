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

   This file implements the function RSPSolverCreate().

   2014-08-06, Bin Gao:
   * first version
*/

#include "eom/openrsp_solver.h"

/*% \brief creates the context of response equation solver, should be called at first
    \author Bin Gao
    \date 2014-08-06
    \param[RSPSolver:struct]{inout} rsp_solver the context of response equation solver
    \param[QVoid:void]{in} user_ctx user-defined callback function context
    \param[GetLinearRSPSolution:void]{in} get_linear_rsp_solution user specified function of
        linear response equation solver
    \return[QErrorCode:int] error information
*/
QErrorCode RSPSolverCreate(RSPSolver *rsp_solver,
#if defined(OPENRSP_C_USER_CONTEXT)
                           QVoid *user_ctx,
#endif
                           const GetLinearRSPSolution get_linear_rsp_solution)
{
#if defined(OPENRSP_C_USER_CONTEXT)
    rsp_solver->user_ctx = user_ctx;
#endif
    rsp_solver->get_linear_rsp_solution = get_linear_rsp_solution;
    return QSUCCESS;
}
