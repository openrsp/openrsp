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
    The header file of overlap operator used inside OpenRSP
  </header>
*/

#if !defined(RSP_OVERLAP_H)
#define RSP_OVERLAP_H

#include "qcmatrix.h"
#include "RSPPerturbation.h"

typedef void (*GetOverlapMat)(const QInt,
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
typedef void (*GetOverlapExp)(const QInt,
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

typedef struct {
    QInt num_pert_lab;              /* number of different perturbation labels
                                       that can act as perturbations on the
                                       overlap operator */
    QInt bra_num_pert;              /* number of perturbations on the bra center,
                                       only used for callback functions */
    QInt ket_num_pert;              /* number of perturbations on the ket center,
                                       only used for callback functions */
    QInt oper_num_pert;             /* number of perturbations on the overlap operator,
                                       only used for callback functions */
    QInt *pert_max_orders;          /* allowed maximal order of a perturbation
                                       described by exactly one of these
                                       different labels */
    QInt *bra_pert_orders;          /* orders of perturbations on the bra center,
                                       only used for callback functions */
    QInt *ket_pert_orders;          /* orders of perturbations on the ket center,
                                       only used for callback functions */
    QInt *oper_pert_orders;         /* orders of perturbations on the overlap operator,
                                       only used for callback functions */
    QcPertInt *pert_labels;         /* all the different perturbation labels */
    QcPertInt *bra_pert_labels;     /* labels of perturbations on the bra center,
                                       only used for callback functions */
    QcPertInt *ket_pert_labels;     /* labels of perturbations on the ket center,
                                       only used for callback functions */
    QcPertInt *oper_pert_labels;    /* labels of perturbations on the overlap operator,
                                       only used for callback functions */
#if defined(OPENRSP_C_USER_CONTEXT)
    void *user_ctx;                 /* user-defined callback-function context */
#endif
    GetOverlapMat get_overlap_mat;  /* user-specified function for calculating
                                       integral matrices */
    GetOverlapExp get_overlap_exp;  /* user-specified function for calculating
                                       expectation values */
} RSPOverlap;

extern QErrorCode RSPOverlapCreate(RSPOverlap*,
                                   const QInt,
                                   const QcPertInt*,
                                   const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   void*,
#endif
                                   const GetOverlapMat,
                                   const GetOverlapExp);
extern QErrorCode RSPOverlapAssemble(RSPOverlap*,const RSPPert*);
extern QErrorCode RSPOverlapWrite(const RSPOverlap*,FILE*);
extern QErrorCode RSPOverlapGetMat(RSPOverlap*,
                                   const QInt,
                                   const QcPertInt*,
                                   const QInt,
                                   const QcPertInt*,
                                   const QInt,
                                   const QcPertInt*,
                                   const QInt,
                                   QcMat*[]);
extern QErrorCode RSPOverlapGetExp(RSPOverlap*,
                                   const QInt,
                                   const QcPertInt*,
                                   const QInt,
                                   const QcPertInt*,
                                   const QInt,
                                   const QcPertInt*,
                                   const QInt,
                                   QcMat*[],
                                   const QInt,
                                   QReal*);
extern QErrorCode RSPOverlapDestroy(RSPOverlap*);

#endif
