/*
  OpenRSP: open-ended library for response theory
  Copyright 2015 Radovan Bast,
                 Daniel H. Friese,
                 Bin Gao,
                 Dan J. Jonsson,
                 Magnus Ringholm,
                 Kenneth Ruud,
                 Andreas Thorvaldsen

  OpenRSP is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.

  OpenRSP is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with OpenRSP. If not, see <http://www.gnu.org/licenses/>.


  <header name='RSPSolver.h' author='Bin Gao' date='2014-08-06'>
    The header file of linear response equation solver used inside OpenRSP
  </header>
*/

#if !defined(RSP_SOLVER_H)
#define RSP_SOLVER_H

#include "qcmatrix.h"

typedef QVoid (*GetLinearRSPSolution)(const QInt,
                                      const QReal*,
                                      const QInt,
                                      QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                                      QVoid*,
#endif
                                      QcMat*[]);

typedef struct {
#if defined(OPENRSP_C_USER_CONTEXT)
    QVoid *user_ctx;                               /* user-defined callback-function
                                                      context */
#endif
    GetLinearRSPSolution get_linear_rsp_solution;  /* user-specified function of
                                                      linear response equation solver */
} RSPSolver;

extern QErrorCode RSPSolverCreate(RSPSolver*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                  QVoid*,
#endif
                                  const GetLinearRSPSolution);
extern QErrorCode RSPSolverAssemble(RSPSolver*);
extern QErrorCode RSPSolverWrite(const RSPSolver*,FILE*);
extern QErrorCode RSPSolverGetLinearRSPSolution(const RSPSolver*,
                                                const QInt,
                                                const QReal*,
                                                const QInt,
                                                QcMat*[],
                                                QcMat*[]);
extern QErrorCode RSPSolverDestroy(RSPSolver*);

#endif
