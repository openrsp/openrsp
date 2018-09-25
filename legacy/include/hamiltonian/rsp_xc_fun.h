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

   This is the header file of exchange-correlation (XC) functionals.

   2014-08-06, Bin Gao:
   * first version
*/

#if !defined(RSP_XC_FUN_H)
#define RSP_XC_FUN_H

/* QcMatrix library */
#include "qcmatrix.h"

/* callback functions to get the integral matrices and expectation values */
typedef QVoid (*GetXCFunMat)(const QInt,
                             const QInt*,
                             const QInt,
                             const QInt,
                             const QInt*,
                             const QInt,
                             QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                             QVoid*,
#endif
                             const QInt,
                             QcMat*[]);
typedef QVoid (*GetXCFunExp)(const QInt,
                             const QInt*,
                             const QInt,
                             const QInt,
                             const QInt*,
                             const QInt,
                             QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                             QVoid*,
#endif
                             const QInt,
                             QReal*);

/* linked list of exchange-corrrelation (XC) functionals */
typedef struct RSPXCFun RSPXCFun;
struct RSPXCFun {
    QInt num_pert;               /* number of different perturbation labels that
                                    can act as perturbations on the XC functional */
    QInt *pert_labels;           /* all the different perturbation labels */
    QInt *pert_max_orders;       /*  maximum allowed order of each perturbation (label) */
#if defined(OPENRSP_C_USER_CONTEXT)
    QVoid *user_ctx;             /* user-defined callback function context */
#endif
    GetXCFunMat get_xc_fun_mat;  /* user specified function for getting integral matrices */
    GetXCFunExp get_xc_fun_exp;  /* user specified function for getting expectation values */
    RSPXCFun *next_xc;           /* pointer to the next exchange-corrrelation functional */
};

/* functions related to the linked list of exchange-corrrelation functionals */
extern QErrorCode RSPXCFunCreate(RSPXCFun**,
                                 const QInt,
                                 const QInt*,
                                 const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                 QVoid*,
#endif
                                 const GetXCFunMat,
                                 const GetXCFunExp);
extern QErrorCode RSPXCFunAdd(RSPXCFun*,
                              const QInt,
                              const QInt*,
                              const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                              QVoid*,
#endif
                              const GetXCFunMat,
                              const GetXCFunExp);
extern QErrorCode RSPXCFunAssemble(RSPXCFun*);
extern QErrorCode RSPXCFunWrite(RSPXCFun*,FILE*);
extern QErrorCode RSPXCFunGetMat(RSPXCFun*,
                                 const QInt,
                                 const QInt*,
                                 const QInt,
                                 const QInt,
                                 const QInt*,
                                 const QInt,
                                 QcMat*[],
                                 const QInt,
                                 QcMat*[]);
extern QErrorCode RSPXCFunGetExp(RSPXCFun*,
                                 const QInt,
                                 const QInt*,
                                 const QInt,
                                 const QInt,
                                 const QInt*,
                                 const QInt,
                                 QcMat*[],
                                 const QInt,
                                 QReal*);
extern QErrorCode RSPXCFunDestroy(RSPXCFun**);

#endif
