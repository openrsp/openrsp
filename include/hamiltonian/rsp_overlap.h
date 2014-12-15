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

   This is the header file of overlap integrals.

   2014-08-05, Bin Gao:
   * first version
*/

#if !defined(RSP_OVERLAP_H)
#define RSP_OVERLAP_H

/* QMatrix library */
#include "qmatrix.h"

/* callback functions to get the integral matrices and expectation values */
typedef QVoid (*GetOverlapMat)(const QInt,
                               const QInt*,
                               const QInt*,
                               const QInt,
                               const QInt*,
                               const QInt*,
                               const QInt,
                               const QInt*,
                               const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                               QVoid*,
#endif
                               const QInt,
                               QMat*[]);
typedef QVoid (*GetOverlapExp)(const QInt,
                               const QInt*,
                               const QInt*,
                               const QInt,
                               const QInt*,
                               const QInt*,
                               const QInt,
                               const QInt*,
                               const QInt*,
                               const QInt,
                               QMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                               QVoid*,
#endif
                               const QInt,
                               QReal*);

/* context of overlap integrals */
typedef struct {
    QInt num_pert;                  /* number of perturbations that the overlap integrals depend on */
    QInt *perturbations;            /* perturbations that the overlap integrals depend on */
    QInt *pert_max_orders;          /* maximum allowed orders of the perturbations */
#if defined(OPENRSP_C_USER_CONTEXT)
    QVoid *user_ctx;                /* user-defined callback function context */
#endif
    GetOverlapMat get_overlap_mat;  /* user specified function for getting integral matrices */
    GetOverlapExp get_overlap_exp;  /* user specified function for getting expectation values */
} RSPOverlap;

/* functions related to the overlap integrals */
extern QErrorCode RSPOverlapCreate(RSPOverlap*,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   QVoid*,
#endif
                                   const GetOverlapMat,
                                   const GetOverlapExp);
extern QErrorCode RSPOverlapAssemble(RSPOverlap*);
extern QErrorCode RSPOverlapWrite(const RSPOverlap*,FILE*);
extern QErrorCode RSPOverlapGetMat(const RSPOverlap*,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
                                   const QReal*,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
                                   const QReal*,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
                                   const QInt,
                                   QMat*[]);
extern QErrorCode RSPOverlapGetExp(const RSPOverlap*,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
                                   const QReal*,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
                                   const QReal*,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
                                   const QInt,
                                   QMat*[],
                                   const QInt,
                                   QReal*);
extern QErrorCode RSPOverlapDestroy(RSPOverlap*);

#endif
