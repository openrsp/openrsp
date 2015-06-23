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

   This file implements the function RSPSolverGetLinearRSPSolution().

   2014-08-06, Bin Gao:
   * first version
*/

#include "eom/openrsp_solver.h"

/*% \brief solves the linear response equation
    \author Bin Gao
    \date 2014-08-06
    \param[RSPSolver:struct]{in} rsp_solver the context of response equation solver
    \param[QInt:int]{in} size_pert size of perturbations acting on the
        time-dependent self-consistent-field (TDSCF) equation
    \param[QInt:int]{in} num_freq_sums number of complex frequency sums
        on the left hand side of the linear response equation
    \param[QReal:real]{in} freq_sums the complex frequency sums on the left hand side
    \param[QcMat:struct]{in} RHS_mat RHS matrices, size is \var{size_pert}*\var{num_freq_sums}
    \param[QcMat:struct]{out} rsp_param solved response parameters,
        size is \var{size_pert}*\var{num_freq_sums}
    \return[QErrorCode:int] error information
*/
QErrorCode RSPSolverGetLinearRSPSolution(const RSPSolver *rsp_solver,
                                         const QInt size_pert,
                                         const QInt num_freq_sums,
                                         const QReal *freq_sums,
                                         QcMat *RHS_mat[],
                                         QcMat *rsp_param[])
{
    rsp_solver->get_linear_rsp_solution(size_pert,
                                        num_freq_sums,
                                        freq_sums,
                                        RHS_mat,
#if defined(OPENRSP_C_USER_CONTEXT)
                                        rsp_solver->user_ctx,
#endif
                                        rsp_param);
    return QSUCCESS;
}
