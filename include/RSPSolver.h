/*
  OpenRSP: open-ended library for response theory
  Copyright 2015 Radovan Bast,
                 Daniel H. Friese,
                 Bin Gao,
                 Dan J. Jonsson,
                 Magnus Ringholm,
                 Kenneth Ruud

  This source code form is subject to the terms of the
  GNU Lesser General Public License, version 2.1.
  If a copy of the GNU LGPL v2.1 was not distributed with this
  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.


  <header name='RSPSolver.h' author='Bin Gao' date='2014-08-06'>
    The header file of linear response equation solver used inside OpenRSP
  </header>
*/

#if !defined(RSP_SOLVER_H)
#define RSP_SOLVER_H

#include "qcmatrix.h"

typedef void (*GetLinearRSPSolution)(const QInt,
                                     const QInt*,
                                     const QInt*,
                                     const QReal*,
                                     QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                                     void*,
#endif
                                     QcMat*[]);

typedef struct {
#if defined(OPENRSP_C_USER_CONTEXT)
    void *user_ctx;                                /* user-defined callback-function
                                                      context */
#endif
    GetLinearRSPSolution get_linear_rsp_solution;  /* user-specified function of
                                                      linear response equation solver */
} RSPSolver;

extern QErrorCode RSPSolverCreate(RSPSolver*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                  void*,
#endif
                                  const GetLinearRSPSolution);
extern QErrorCode RSPSolverAssemble(RSPSolver*);
extern QErrorCode RSPSolverWrite(const RSPSolver*,FILE*);
extern QErrorCode RSPSolverGetLinearRSPSolution(const RSPSolver*,
                                                const QInt,
                                                const QInt*,
                                                const QInt*,
                                                const QReal*,
                                                QcMat*[],
                                                QcMat*[]);
extern QErrorCode RSPSolverDestroy(RSPSolver*);

#endif
