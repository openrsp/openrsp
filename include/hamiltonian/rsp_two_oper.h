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

   This is the header file of two-electron operators.

   2014-08-05, Bin Gao:
   * first version
*/

#if !defined(RSP_TWO_OPER_H)
#define RSP_TWO_OPER_H

/* QMatrix library */
#include "qmatrix.h"

/* callback functions to get the integral matrices and expectation values */
typedef QVoid (*GetTwoOperMat)(const QInt,
                               const QInt*,
                               const QInt*,
                               const QInt,
                               QMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                               QVoid*,
#endif
                               const QInt,
                               QMat*[]);
typedef QVoid (*GetTwoOperExp)(const QInt,
                               const QInt*,
                               const QInt*,
                               const QInt,
                               QMat*[],
                               const QInt,
                               QMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                               QVoid*,
#endif
                               const QInt,
                               QReal*);

/* linked list of two-electron operators */
typedef struct RSPTwoOper RSPTwoOper;
struct RSPTwoOper {
    QInt num_pert;                   /* number of perturbations that the two-electron operator depends on */
    QInt *perturbations;             /* perturbations that the two-electron operator depends on */
    QInt *pert_max_orders;           /* maximum allowed orders of the perturbations */
#if defined(OPENRSP_C_USER_CONTEXT)
    QVoid *user_ctx;                 /* user-defined callback function context */
#endif
    GetTwoOperMat get_two_oper_mat;  /* user specified function for getting integral matrices */
    GetTwoOperExp get_two_oper_exp;  /* user specified function for getting expectation values */
    RSPTwoOper *next_oper;           /* pointer to the next two-electron operator */
};

/* functions related to the linked list of two-electron operators */
extern QErrorCode RSPTwoOperCreate(RSPTwoOper**,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   QVoid*,
#endif
                                   const GetTwoOperMat,
                                   const GetTwoOperExp);
extern QErrorCode RSPTwoOperAdd(RSPTwoOper*,
                                const QInt,
                                const QInt*,
                                const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                QVoid*,
#endif
                                const GetTwoOperMat,
                                const GetTwoOperExp);
extern QErrorCode RSPTwoOperAssemble(RSPTwoOper*);
extern QErrorCode RSPTwoOperWrite(RSPTwoOper*,FILE*);
extern QErrorCode RSPTwoOperGetMat(RSPTwoOper*,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
                                   const QInt,
                                   QMat*[],
                                   const QInt,
                                   QMat*[]);
extern QErrorCode RSPTwoOperGetExp(RSPTwoOper*,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
                                   const QInt,
                                   QMat*[],
                                   const QInt,
                                   QMat*[],
                                   const QInt,
                                   QReal*);
extern QErrorCode RSPTwoOperDestroy(RSPTwoOper**);

#endif