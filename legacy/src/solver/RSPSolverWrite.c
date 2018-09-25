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

   This file implements the function RSPSolverWrite().

   2014-08-06, Bin Gao:
   * first version
*/

#include "eom/openrsp_solver.h"

/*% \brief writes the context of response equation solver
    \author Bin Gao
    \date 2014-08-06
    \param[RSPSolver:struct]{in} rsp_solver the context of response equation solver
    \param[FILE]{inout} fp_solver file pointer
    \return[QErrorCode:int] error information
*/
QErrorCode RSPSolverWrite(const RSPSolver *rsp_solver, FILE *fp_solver)
{
#if defined(OPENRSP_C_USER_CONTEXT)
    if (rsp_solver->user_ctx!=NULL) {
        fprintf(fp_solver, "RSPSolverWrite>> user-defined function context given\n");
    }
#endif
    return QSUCCESS;
}
