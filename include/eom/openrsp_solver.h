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

   This is the header file of response equation solver.

   2014-08-06, Bin Gao:
   * first version
*/

#if !defined(OPENRSP_SOLVER_H)
#define OPENRSP_SOLVER_H

/* QcMatrix library */
#include "qcmatrix.h"

/* callback function to get the solutions of response equation */
typedef QVoid (*GetRSPSolution)(const QcMat*,
                                const QcMat*,
                                const QcMat*,
                                const QInt,
                                const QReal*,
                                const QInt,
                                QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                                QVoid*,
#endif
                                QcMat*[]);

/* context of response equation solver */
typedef struct {
#if defined(OPENRSP_C_USER_CONTEXT)
    QVoid *user_ctx;                 /* user-defined callback function context */
#endif
    GetRSPSolution get_rsp_solution; /* user specified function of response equation solver */
} RSPSolver;

/* functions related to the linked list of one-electron operators */
extern QErrorCode RSPSolverCreate(RSPSolver*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                  QVoid*,
#endif
                                  const GetRSPSolution);
extern QErrorCode RSPSolverAssemble(RSPSolver*);
extern QErrorCode RSPSolverWrite(const RSPSolver*,FILE*);
extern QErrorCode RSPSolverGetSolution(const RSPSolver*,
                                       const QcMat*,
                                       const QcMat*,
                                       const QcMat*,
                                       const QInt,
                                       const QReal*,
                                       const QInt,
                                       QcMat*[],
                                       QcMat*[]);
extern QErrorCode RSPSolverDestroy(RSPSolver*);

#endif
