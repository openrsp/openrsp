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

   This is the header file of perturbations.

   2015-06-23, Bin Gao:
   * first version
*/

#if !defined(OPENRSP_PERTURBATION_H)
#define OPENRSP_PERTURBATION_H

/* QcMatrix library */
#include "qcmatrix.h"

/* callback function to get the components of a perturbation */
typedef QVoid (*GetPertComp)(const QInt,
                             const QInt,
                             const QInt,
#if defined(OPENRSP_C_USER_CONTEXT)
                             QVoid*,
#endif
                             QInt*,
                             QInt*,
                             QInt*);
/* callback function to get the rank of a perturbation */
typedef QVoid (*GetPertRank)(const QInt,
                             const QInt,
                             const QInt*,
                             const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                             QVoid*,
#endif
                             QInt*);

/* context of all perturbations involved in calculations */
typedef struct {
    /* perturbations */
    QInt num_pert;                 /* number of all perturbations involved in calculations */
    QInt *pert_labels;             /* labels of all perturbations */
    QInt *pert_max_orders;         /* maximum allowed orders of all perturbations */
    QInt *size_ptr;                /* pointer to the size of each perturbation */
    QInt *pert_sizes;              /* sizes of all perturbations up to their maximum orders */
#if defined(OPENRSP_C_USER_CONTEXT)
    QVoid *user_ctx;               /* user-defined callback function context */
#endif
    GetPertComp get_pert_comp;     /* user specified function for getting components of a perturbation */
    GetPertRank get_pert_rank;     /* user specified function for getting rank of a perturbation */
} RSPPert;

/* functions related to the perturbations */
extern QErrorCode RSPSolverCreate(RSPSolver*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                  QVoid*,
#endif
                                  const GetLinearRSPSolution);
extern QErrorCode RSPSolverAssemble(RSPSolver*);
extern QErrorCode RSPSolverWrite(const RSPSolver*,FILE*);
extern QErrorCode RSPSolverGetLinearRSPSolution(const RSPSolver*,
                                                const QInt,
                                                const QInt,
                                                const QReal*,
                                                QcMat*[],
                                                QcMat*[]);
extern QErrorCode RSPSolverDestroy(RSPSolver*);

#endif
