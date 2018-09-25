/*
  OpenRSP: open-ended library for response theory
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


  <header name='RSPZeroOper.h' author='Bin Gao' date='2014-12-11'>
    The header file of zero-electron operator used inside OpenRSP
  </header>
*/

#if !defined(RSP_ELEC_FREE_OPER_H)
#define RSP_ELEC_FREE_OPER_H

#include "qcmatrix.h"
#include "RSPPerturbation.h"

typedef void (*GetZeroOperContrib)(const QInt,
                                   const QcPertInt*,
                                   const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   void*,
#endif
                                   const QInt,
                                   QReal*);

typedef struct RSPZeroOper RSPZeroOper;
struct RSPZeroOper {
    QInt num_pert_lab;              /* number of different perturbation labels
                                       that can act as perturbations on the
                                       zero-electron operator */
    QInt oper_num_pert;             /* number of perturbations on the
                                       zero-electron operator, only used for
                                       callback functions */
    QInt *pert_max_orders;          /* allowed maximal order of a perturbation
                                       described by exactly one of these
                                       different labels */
    QInt *oper_pert_orders;         /* orders of perturbations on the
                                       zero-electron operator, only used for
                                       callback functions */
    QcPertInt *pert_labels;         /* all the different perturbation labels */
    QcPertInt *oper_pert_labels;    /* labels of perturbations on the
                                       zero-electron operator, only used for
                                       callback functions */
#if defined(OPENRSP_C_USER_CONTEXT)
    void *user_ctx;                 /* user-defined callback-function context */
#endif
    GetZeroOperContrib get_zero_oper_contrib;  /* user-specified function for calculating
                                       contribution from the zero-electron operator */
    RSPZeroOper *next_oper;           /* pointer to the next zero-electron operator */
};

extern QErrorCode RSPZeroOperCreate(RSPZeroOper**,
                                    const QInt,
                                    const QcPertInt*,
                                    const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                    void*,
#endif
                                    const GetZeroOperContrib);
extern QErrorCode RSPZeroOperAdd(RSPZeroOper*,
                                 const QInt,
                                 const QcPertInt*,
                                 const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                 void*,
#endif
                                 const GetZeroOperContrib);
extern QErrorCode RSPZeroOperAssemble(RSPZeroOper*,const RSPPert*);
extern QErrorCode RSPZeroOperWrite(RSPZeroOper*,FILE*);
extern QErrorCode RSPZeroOperGetContribution(RSPZeroOper*,
                                             const QInt,
                                             const QcPertInt*,
                                             const QInt,
                                             QReal*);
extern QErrorCode RSPZeroOperDestroy(RSPZeroOper**);
/*FIXME: RSPZeroOperGetNumAtoms() to be removed after perturbation free scheme implemented*/
extern QErrorCode RSPZeroOperGetNumAtoms(const RSPZeroOper*,QInt*);

#endif
