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

   This file implements the function RSPSolverGetLinearRSPSolution().

   2014-08-06, Bin Gao:
   * first version
*/

#include "eom/openrsp_solver.h"

/*% \brief solves the linear response equation
    \author Bin Gao
    \date 2014-08-06
    \param[RSPSolver:struct]{in} rsp_solver the context of response equation solver
    \param[QInt:int]{in} num_freq_sums number of complex frequency sums
        on the left hand side of the linear response equation
    \param[QReal:real]{in} freq_sums the complex frequency sums on the left hand side
    \param[QInt:int]{in} size_pert size of perturbations acting on the
        time-dependent self-consistent-field (TDSCF) equation
    \param[QcMat:struct]{in} RHS_mat RHS matrices, size is \var{size_pert}*\var{num_freq_sums}
    \param[QcMat:struct]{out} rsp_param solved response parameters,
        size is \var{size_pert}*\var{num_freq_sums}
    \return[QErrorCode:int] error information
*/
QErrorCode RSPSolverGetLinearRSPSolution(const RSPSolver *rsp_solver,
                                         const QInt num_freq_sums,
                                         const QReal *freq_sums,
                                         const QInt size_pert,
                                         QcMat *RHS_mat[],
                                         QcMat *rsp_param[])
{
    rsp_solver->get_linear_rsp_solution(num_freq_sums,
                                        freq_sums,
                                        size_pert,
                                        RHS_mat,
#if defined(OPENRSP_C_USER_CONTEXT)
                                        rsp_solver->user_ctx,
#endif
                                        rsp_param);
    return QSUCCESS;
}
