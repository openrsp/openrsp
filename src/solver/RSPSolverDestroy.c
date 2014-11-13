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

   This file implements the function RSPSolverDestroy().

   2014-08-06, Bin Gao:
   * first version
*/

#include "eom/openrsp_solver.h"

/*% \brief destroys the context of response equation solver, should be called at the end
    \author Bin Gao
    \date 2014-08-06
    \param[RSPSolver:struct]{inout} rsp_solver the context of response equation solver
    \return[QErrorCode:int] error information
*/
QErrorCode RSPSolverDestroy(RSPSolver *rsp_solver)
{
#if defined(OPENRSP_C_USER_CONTEXT)
    rsp_solver->user_ctx = NULL;
#endif
    rsp_solver->get_rsp_solution = NULL;
    return QSUCCESS;
}
