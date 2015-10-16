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


  <header name='RSPTwoOper.h' author='Bin Gao' date='2014-08-05'>
    The header file of two-electron operators used inside OpenRSP
  </header>
*/

#if !defined(RSP_TWOOPER_H)
#define RSP_TWOOPER_H

#include "qcmatrix.h"
#include "RSPPerturbation.h"

typedef void (*GetTwoOperMat)(const QInt,
                              const QcPertInt*,
                              const QInt*,
                              const QInt,
                              QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                              void*,
#endif
                              const QInt,
                              QcMat*[]);
typedef void (*GetTwoOperExp)(const QInt,
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

typedef struct RSPTwoOper RSPTwoOper;
struct RSPTwoOper {
    QInt num_pert_lab;               /* number of different perturbation labels
                                        that can act as perturbations on the
                                        two-electron operator */
    QInt oper_num_pert;              /* number of perturbations on the
                                        two-electron operator, only used for
                                        callback functions */
    QInt *pert_max_orders;           /* allowed maximal order of a perturbation
                                        described by exactly one of these
                                        different labels */
    QInt *oper_pert_orders;          /* orders of perturbations on the
                                        two-electron operator, only used for
                                        callback functions */
    QcPertInt *pert_labels;          /* all the different perturbation labels */
    QcPertInt *oper_pert_labels;     /* labels of perturbations on the
                                        two-electron operator, only used for
                                        callback functions */
#if defined(OPENRSP_C_USER_CONTEXT)
    void *user_ctx;                  /* user-defined callback-function context */
#endif
    GetTwoOperMat get_two_oper_mat;  /* user-specified function for calculating
                                        integral matrices */
    GetTwoOperExp get_two_oper_exp;  /* user-specified function for calculating
                                        expectation values */
    RSPTwoOper *next_oper;           /* pointer to the next two-electron operator */
};

extern QErrorCode RSPTwoOperCreate(RSPTwoOper**,
                                   const QInt,
                                   const QcPertInt*,
                                   const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   void*,
#endif
                                   const GetTwoOperMat,
                                   const GetTwoOperExp);
extern QErrorCode RSPTwoOperAdd(RSPTwoOper*,
                                const QInt,
                                const QcPertInt*,
                                const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                void*,
#endif
                                const GetTwoOperMat,
                                const GetTwoOperExp);
extern QErrorCode RSPTwoOperAssemble(RSPTwoOper*,const RSPPert*);
extern QErrorCode RSPTwoOperWrite(RSPTwoOper*,FILE*);
extern QErrorCode RSPTwoOperGetMat(RSPTwoOper*,
                                   const QInt,
                                   const QcPertInt*,
                                   const QInt,
                                   QcMat*[],
                                   const QInt,
                                   QcMat*[]);
extern QErrorCode RSPTwoOperGetExp(RSPTwoOper*,
                                   const QInt,
                                   const QcPertInt*,
                                   const QInt,
                                   const QInt*,
                                   QcMat*[],
                                   const QInt*,
                                   QcMat*[],
                                   const QInt,
                                   QReal*);
extern QErrorCode RSPTwoOperDestroy(RSPTwoOper**);

#endif
