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

   This is the header file of linear response equation solver.

   2014-08-06, Bin Gao:
   * first version
*/

#if !defined(OPENRSP_SOLVER_H)
#define OPENRSP_SOLVER_H

/* QcMatrix library */
#include "qcmatrix.h"

/* callback function of linear response equation solver */
typedef QVoid (*GetLinearRSPSolution)(const QInt,
                                      const QReal*,
                                      const QInt,
                                      QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                                      QVoid*,
#endif
                                      QcMat*[]);

/* context of linear response equation solver */
typedef struct {
#if defined(OPENRSP_C_USER_CONTEXT)
    QVoid *user_ctx;                              /* user-defined callback function context */
#endif
    GetLinearRSPSolution get_linear_rsp_solution; /* user specified function of linear response equation solver */
} RSPSolver;

/* functions related to the linear response equation solver */
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
