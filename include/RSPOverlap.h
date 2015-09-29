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


  <header name='RSPOneOper.h' author='Bin Gao' date='2014-08-05'>
    The header file of overlap integrals used inside OpenRSP
  </header>
*/

#if !defined(RSP_OVERLAP_H)
#define RSP_OVERLAP_H

#include "qcmatrix.h"

typedef QVoid (*GetOverlapMat)(const QInt,
                               const QInt*,
                               const QInt,
                               const QInt*,
                               const QInt,
                               const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                               QVoid*,
#endif
                               const QInt,
                               QcMat*[]);
typedef QVoid (*GetOverlapExp)(const QInt,
                               const QInt*,
                               const QInt,
                               const QInt*,
                               const QInt,
                               const QInt*,
                               const QInt,
                               QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                               QVoid*,
#endif
                               const QInt,
                               QReal*);

typedef struct {
    QInt num_pert;                  /* number of different perturbation labels that
                                       can act as perturbations on the basis sets */
    QInt *pert_labels;              /* all the different perturbation labels */
    QInt *pert_max_orders;          /* maximum allowed order of each perturbation (label) */
#if defined(OPENRSP_C_USER_CONTEXT)
    QVoid *user_ctx;                /* user-defined callback function context */
#endif
    GetOverlapMat get_overlap_mat;  /* user specified function for getting integral matrices */
    GetOverlapExp get_overlap_exp;  /* user specified function for getting expectation values */
} RSPOverlap;

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
                                   const QInt,
                                   const QInt*,
                                   const QInt,
                                   const QInt*,
                                   const QInt,
                                   QcMat*[]);
extern QErrorCode RSPOverlapGetExp(const RSPOverlap*,
                                   const QInt,
                                   const QInt*,
                                   const QInt,
                                   const QInt*,
                                   const QInt,
                                   const QInt*,
                                   const QInt,
                                   QcMat*[],
                                   const QInt,
                                   QReal*);
extern QErrorCode RSPOverlapDestroy(RSPOverlap*);

#endif
