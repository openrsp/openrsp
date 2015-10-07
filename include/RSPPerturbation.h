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


  <header name='RSPPerturbation.h' author='Bin Gao' date='2015-06-23'>
    The header file of perturbations used inside OpenRSP
  </header>
*/

#if !defined(RSP_PERTURBATION_H)
#define RSP_PERTURBATION_H

/* QcMatrix library */
#include "qcmatrix.h"

/* <macrodef name='OPENRSP_PERT_SHORT_INT'>
     Represent perturbation labels using unsigned short integers
   </macrodef> */
#if defined(OPENRSP_PERT_SHORT_INT)
/* <datatype name='QcPertInt'>
     Data type of integers to represent perturbation labels
   </datatype>
   <constant name='QCPERTINT_MAX'>
     Maximal value of an object of the <QcPertInt> type
   </constant>
   <constant name='QCPERTINT_FMT'>
     Format string of <QcPertInt> type
   </constant> */
typedef unsigned short QcPertInt;
#define QCPERTINT_MAX USHRT_MAX
#define QCPERTINT_FMT "hu"
/* <macrodef name='OPENRSP_PERT_INT'>
     Represent perturbation labels using unsigned integers
   </macrodef> */
#elif defined(OPENRSP_PERT_INT)
typedef unsigned int QcPertInt;
#define QCPERTINT_MAX UINT_MAX
#define QCPERTINT_FMT "u"
#else
typedef unsigned long QcPertInt;
#define QCPERTINT_MAX ULONG_MAX
#define QCPERTINT_FMT "lu"
#endif
/* <macrodef name='OPENRSP_PERT_LABEL_BIT'>
     Set <OPENRSP_PERT_LABEL_BIT>
   </macrodef>
   <constant name='OPENRSP_PERT_LABEL_BIT'>
     Number of bits in an object of <QcPertInt> type for a perturbation label
   </constant> */
#if !defined(OPENRSP_PERT_LABEL_BIT)
#define OPENRSP_PERT_LABEL_BIT 10
#endif
/* <constant name='OPENRSP_PERT_LABEL_MAX'>
     Maximal value for perturbation labels
   </constant>
   <constant name='OPENRSP_NUM_FREQ_MAX'>
     Maximal value for number of frequencies
   </constant> */
extern const QcPertInt OPENRSP_PERT_LABEL_MAX;
extern const QcPertInt OPENRSP_NUM_FREQ_MAX;

typedef QVoid (*GetPertCat)(const QcPertInt,
                            const QInt,
                            const QInt,
                            const QInt,
                            const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                            QVoid*,
#endif
                            QInt*);

typedef struct {
    QcPertInt num_pert;                 /* number of different perturbation labels $p$ */
    QcPertInt *pert_labels;             /* $a_{1},a_{2},\cdots,a_{p}$ */
    QInt *pert_max_orders;              /* $n_{1},n_{2},\cdots,n_{p}$ */
    QInt *ptr_ncomp;                    /* pointers to $[N_{j}^{k_{j}}]$
                                           for each $a_{j}$ */
    QInt *pert_num_comps;               /* $[N_{j}^{k_{j}}]$, where
                                           $1\le k_{j}\le n_{j}$ and $1\le j\le p$ */
#if defined(OPENRSP_C_USER_CONTEXT)     
    QVoid *user_ctx;                    /* user-defined callback function context */
#endif
    GetPertCat get_pert_concatenation;  /* user specified function for getting
                                           the ranks of components of sub-perturbation
                                           tuples (with the same perturbation label)
                                           for given components of the corresponding
                                           concatenated perturbation tuple */
} RSPPert;

extern QErrorCode RSPPertCreate(RSPPert*,
                                const QcPertInt,
                                const QcPertInt*,
                                const QInt*,
                                const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                QVoid*,
#endif
                                const GetPertCat);
extern QErrorCode RSPPertAssemble(RSPPert*);
extern QErrorCode RSPPertWrite(const RSPPert*,FILE*);
//extern QErrorCode RSPPertGetFromTuple(const RSPPert*,
//                                      const QInt,
//                                      const QcPertInt*,
//                                      const QcPertInt,
//                                      const QReal*,
//                                      QInt*);
extern QErrorCode RSPPertGetConcatenation(const RSPPert*,
                                          const QcPertInt,
                                          const QInt,
                                          const QInt,
                                          const QInt,
                                          const QInt*,
                                          QInt*);
extern QErrorCode RSPPertDestroy(RSPPert*);

#endif
