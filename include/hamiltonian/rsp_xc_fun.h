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

   This is the header file of exchange-correlation functionals.

   2014-08-06, Bin Gao:
   * first version
*/

#if !defined(RSP_XC_FUN_H)
#define RSP_XC_FUN_H

/* QMatrix library */
#include "qmatrix.h"

/* callback functions to get the integral matrices and expectation values */
typedef QVoid (*GetXCFunMat)(const QInt,
                             const QInt*,
                             const QInt*,
                             const QInt,
                             QMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                             QVoid*,
#endif
                             const QInt,
                             QMat*[]);
typedef QVoid (*GetXCFunExp)(const QInt,
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

/* linked list of exchange-corrrelation functionals */
typedef struct RSPXCFun RSPXCFun;
struct RSPXCFun {
    QInt num_pert;               /* number of perturbations that the exchange-corrrelation functional depends on */
    QInt *perturbations;         /* perturbations that the exchange-corrrelation functional depends on */
    QInt *pert_max_orders;       /* maximum allowed orders of the perturbations */
#if defined(OPENRSP_C_USER_CONTEXT)
    QVoid *user_ctx;             /* user-defined callback function context */
#endif
    GetXCFunMat get_xc_fun_mat;  /* user specified function for getting integral matrices */
    GetXCFunExp get_xc_fun_exp;  /* user specified function for getting expectation values */
    RSPXCFun *next_oper;         /* pointer to the next exchange-corrrelation functional */
};

/* functions related to the linked list of exchange-corrrelation functionals */
//extern QErrorCode RSPXCFunCreate(RSPXCFun**,
//                                 const QInt,
//                                 const QInt*,
//                                 const QInt*,
//#if defined(OPENRSP_C_USER_CONTEXT)
//                                 QVoid*,
//#endif
//                                 const GetXCFunMat,
//                                 const GetXCFunExp);
//extern QErrorCode RSPXCFunAdd(RSPXCFun*,
//                              const QInt,
//                              const QInt*,
//                              const QInt*,
//#if defined(OPENRSP_C_USER_CONTEXT)
//                              QVoid*,
//#endif
//                              const GetXCFunMat,
//                              const GetXCFunExp);
//extern QErrorCode RSPXCFunAssemble(RSPXCFun*);
//extern QErrorCode RSPXCFunWrite(const RSPXCFun*,FILE*);
//extern QErrorCode RSPXCFunGetMat(const RSPXCFun*,
//                                 const QInt,
//                                 const QInt*,
//                                 const QInt*,
//                                 const QInt,
//                                 QMat*[],
//                                 const QInt,
//                                 QMat*[]);
//extern QErrorCode RSPXCFunGetExp(const RSPXCFun*,
//                                 const QInt,
//                                 const QInt*,
//                                 const QInt*,
//                                 const QInt,
//                                 QMat*[],
//                                 const QInt,
//                                 QMat*[],
//                                 const QInt,
//                                 QReal*);
//extern QErrorCode RSPXCFunDestroy(RSPXCFun**);

#endif
