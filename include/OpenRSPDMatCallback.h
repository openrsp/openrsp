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


  This header file contains callback functions for AO-based response theory
  calculations.

  2015-10-16, Bin Gao:
  * first version
*/

#if !defined(OPENRSP_DMAT_CALLBACK_H)
#define OPENRSP_DMAT_CALLBACK_H

#include "OpenRSP.h"
#include "OpenRSPTestData.h"

#if defined(OPENRSP_C_USER_CONTEXT)
#define OVERLAP_CONTEXT "OVERLAP"
#define ONEHAM_CONTEXT "ONEHAM"
#define EXT_FIELD_CONTEXT "EXT_FIELD"
#define TWO_OPER_CONTEXT "NONLAO"
#define SOLVER_CONTEXT "SOLVER"
#endif

extern void get_overlap_mat(const QInt,
                            const QcPertInt*,
                            const QInt*,
                            const QInt,
                            const QcPertInt*,
                            const QInt*,
                            const QInt,
                            const QcPertInt*,
                            const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                            void*,
#endif
                            const QInt,
                            QcMat*[]);
extern void get_overlap_exp(const QInt,
                            const QcPertInt*,
                            const QInt*,
                            const QInt,
                            const QcPertInt*,
                            const QInt*,
                            const QInt,
                            const QcPertInt*,
                            const QInt*,
                            const QInt,
                            QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                            void*,
#endif
                            const QInt,
                            QReal*);
extern void get_one_oper_mat(const QInt,
                             const QcPertInt*,
                             const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                             void*,
#endif
                             const QInt,
                             QcMat*[]);
extern void get_one_oper_exp(const QInt,
                             const QcPertInt*,
                             const QInt*,
                             const QInt,
                             QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                             void*,
#endif
                             const QInt,
                             QReal*);
extern void get_two_oper_mat(const QInt,
                             const QcPertInt*,
                             const QInt*,
                             const QInt,
                             QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                             void*,
#endif
                             const QInt,
                             QcMat*[]);
extern void get_two_oper_exp(const QInt,
                             const QcPertInt*,
                             const QInt*,
                             const QInt,
                             const QInt*,
                             QcMat*[],
                             const QInt*,
                             QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                             void*,
#endif
                             const QInt,
                             QReal*);
extern void get_linear_rsp_solution(const QInt,
                                    const QInt*,
                                    const QInt*,
                                    const QReal*,
                                    QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                                    void*,
#endif
                                    QcMat*[]);

#endif
