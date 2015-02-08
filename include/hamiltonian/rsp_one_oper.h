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

   This is the header file of one-electron operators.

   2014-07-30, Bin Gao:
   * first version
*/

#if !defined(RSP_ONE_OPER_H)
#define RSP_ONE_OPER_H

/* QcMatrix library */
#include "qcmatrix.h"

/* callback functions to get the integral matrices and expectation values */
typedef QVoid (*GetOneOperMat)(const QInt,
                               const QInt*,
                               const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                               QVoid*,
#endif
                               const QInt,
                               QcMat*[]);
typedef QVoid (*GetOneOperExp)(const QInt,
                               const QInt*,
                               const QInt*,
                               const QInt,
                               QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                               QVoid*,
#endif
                               const QInt,
                               QReal*);

/* linked list of one-electron operators */
typedef struct RSPOneOper RSPOneOper;
struct RSPOneOper {
    QInt num_pert;                   /* number of perturbations that the one-electron operator depends on */
    QInt *pert_labels;               /* labels of the perturbations */
    QInt *pert_max_orders;           /* maximum allowed orders of the perturbations */
#if defined(OPENRSP_C_USER_CONTEXT)
    QVoid *user_ctx;                 /* user-defined callback function context */
#endif
    GetOneOperMat get_one_oper_mat;  /* user specified function for getting integral matrices */
    GetOneOperExp get_one_oper_exp;  /* user specified function for getting expectation values */
    RSPOneOper *next_oper;           /* pointer to the next one-electron operator */
};

/* functions related to the linked list of one-electron operators */
extern QErrorCode RSPOneOperCreate(RSPOneOper**,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   QVoid*,
#endif
                                   const GetOneOperMat,
                                   const GetOneOperExp);
extern QErrorCode RSPOneOperAdd(RSPOneOper*,
                                const QInt,
                                const QInt*,
                                const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                QVoid*,
#endif
                                const GetOneOperMat,
                                const GetOneOperExp);
extern QErrorCode RSPOneOperAssemble(RSPOneOper*);
extern QErrorCode RSPOneOperWrite(RSPOneOper*,FILE*);
extern QErrorCode RSPOneOperGetMat(RSPOneOper*,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
                                   const QInt,
                                   QcMat*[]);
extern QErrorCode RSPOneOperGetExp(RSPOneOper*,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
                                   const QInt,
                                   QcMat*[],
                                   const QInt,
                                   QReal*);
extern QErrorCode RSPOneOperDestroy(RSPOneOper**);

#endif
