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


  <header name='RSPPerturbation.h' author='Bin Gao' date='2015-06-23'>
    The header file of perturbations used inside OpenRSP
  </header>
*/

#if !defined(RSP_PERTURBATION_H)
#define RSP_PERTURBATION_H

/* QcMatrix library */
#include "qcmatrix.h"

/* <macrodef name='OPENRSP_PERT_LABEL_BIT'>
     Set <OPENRSP_PERT_LABEL_BIT>
   </macrodef>
   <constant name='OPENRSP_PERT_LABEL_BIT'>
     Number of bits in an object of <QcPertInt> type for a perturbation label
   </constant> */
#if !defined(OPENRSP_PERT_LABEL_BIT)
#define OPENRSP_PERT_LABEL_BIT 10
#endif

/* <datatype name='QcPertInt'>
     Data type of integers to represent perturbation labels
   </datatype>
   <constant name='QCPERTINT_MAX'>
     Maximal value of an object of the <QcPertInt> type
   </constant>
   <constant name='QCPERTINT_FMT'>
     Format string of <QcPertInt> type
   </constant> */
//typedef unsigned long QcPertInt;
//#define QCPERTINT_MAX ULONG_MAX
//#define QCPERTINT_FMT "lu"
//typedef unsigned int QcPertInt;
//#define QCPERTINT_MAX UINT_MAX
//#define QCPERTINT_FMT "u"
typedef QInt QcPertInt;
#define QCPERTINT_MAX INT_MAX
#define QCPERTINT_FMT QINT_FMT
/* <constant name='OPENRSP_PERT_LABEL_MAX'>
     Maximal value for perturbation labels
   </constant>
   <constant name='OPENRSP_PERT_ID_MAX'>
     Maximal value for internal perturbation identifier
   </constant> */
extern const QcPertInt OPENRSP_PERT_LABEL_MAX;
extern const QcPertInt OPENRSP_PERT_ID_MAX;

typedef void (*GetPertCat)(const QInt,
                           const QcPertInt,
                           const QInt,
                           const QInt,
                           const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                           void*,
#endif
                           QInt*);

typedef struct {
    QInt num_pert_lab;                  /* number of different perturbation
                                           labels $p$ */
    QInt *pert_max_orders;              /* $n_{1},n_{2},\cdots,n_{p}$ */
    QInt *ptr_ncomp;                    /* pointers to $[N_{j}^{k_{j}}]$
                                           for each $a_{j}$ */
    QInt *pert_num_comps;               /* $[N_{j}^{k_{j}}]$, where
                                           $1\le k_{j}\le n_{j}$ and $1\le j\le p$ */
    QcPertInt *pert_labels;             /* $a_{1},a_{2},\cdots,a_{p}$ */
#if defined(OPENRSP_C_USER_CONTEXT)
    void *user_ctx;                     /* user-defined callback function context */
#endif
    GetPertCat get_pert_concatenation;  /* user-specified function for getting
                                           the ranks of components of sub-perturbation
                                           tuples (with the same perturbation label)
                                           for given components of the corresponding
                                           concatenated perturbation tuple */
} RSPPert;

extern QErrorCode RSPPertCreate(RSPPert*,
                                const QInt,
                                const QcPertInt*,
                                const QInt*,
                                const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                void*,
#endif
                                const GetPertCat);
extern QErrorCode RSPPertAssemble(RSPPert*);
extern QErrorCode RSPPertWrite(const RSPPert*,FILE*);
extern QErrorCode RSPPertDestroy(RSPPert*);
extern QErrorCode RSPPertValidateLabelOrder(const RSPPert*,
                                            const QInt,
                                            const QcPertInt*,
                                            const QInt*);
extern QErrorCode RSPPertHostToInternal(const RSPPert*,
                                        const QInt,
                                        QcPertInt*,
                                        const QInt,
                                        QReal*);
extern QErrorCode RSPPertInternTupleToHostLabelOrder(const QInt,
                                                     const QcPertInt*,
                                                     const QInt,
                                                     const QcPertInt*,
                                                     const QInt*,
                                                     QInt*,
                                                     QcPertInt*,
                                                     QInt*);
extern QErrorCode RSPPertInternTupleToHostTuple(const QInt,
                                                const QcPertInt*,
                                                const QInt,
                                                const QcPertInt*,
                                                const QInt*,
                                                QInt*,
                                                QcPertInt*);
extern QErrorCode RSPPertGetConcatenation(const RSPPert*,
                                          const QcPertInt,
                                          const QInt,
                                          const QInt,
                                          const QInt,
                                          const QInt*,
                                          QInt*);

#endif
