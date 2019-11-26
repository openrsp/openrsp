/*
  OpenRSP: open-ended library for response theory
  Copyright 2015 Radovan Bast,
                 Daniel H. Friese,
                 Bin Gao,
                 Dan J. Jonsson,
                 Magnus Ringholm,
                 Kenneth Ruud

  This source code form is subject to the terms of the
  GNU Lesser General Public License, version 2.1.
  If a copy of the GNU LGPL v2.1 was not distributed with this
  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.


  <header name='OpenRSP.h' author='Bin Gao' date='2014-01-27'>
    The header file of OpenRSP library for users
  </header>
*/

#if !defined(OPENRSP_H)
#define OPENRSP_H

/* host program perturbations */
#include "RSPPerturbation.h"
/* type of electronic wave function */
/*#include "RSPWaveFunction.h"*/
/* overlap integrals */
#include "RSPOverlap.h"
/* one-electron operators */
#include "RSPOneOper.h"
/* two-electron operators */
#include "RSPTwoOper.h"
/* exchange-correlation (XC) functionals */
#include "RSPXCFun.h"
/* zero-electron operators */
#include "RSPZeroOper.h"
/* linear response equation solver */
#include "RSPSolver.h"

typedef struct {
    QBool assembled;               /* indicates if the OpenRSP context assembled */
    RSPPert *rsp_pert;             /* host program perturbations */
    /*ElecWav *elec_wav;*/           /* implementation-specific data of (electronic) wave function */
    /*ElecWavType elec_wav_type;*/
    RSPOverlap *overlap;           /* overlap integrals */
    RSPOneOper *one_oper;          /* one-electron operators */
    RSPTwoOper *two_oper;          /* two-electron operators */
    RSPXCFun *xc_fun;              /* XC functionals */
    RSPZeroOper *zero_oper;        /* zero-electron operators */
    RSPSolver *rsp_solver;         /* linear response equation solver */
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
    QInt num_atoms;
} OpenRSP;

extern QErrorCode OpenRSPCreate(OpenRSP*,const QInt);
extern QErrorCode OpenRSPSetPerturbations(OpenRSP*,
                                          const QInt,
                                          const QcPertInt*,
                                          const QInt*,
                                          const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                          void*,
#endif
                                          const GetPertCat);
/*extern QErrorCode OpenRSPSetWaveFunction(OpenRSP*,const ElecWavType);*/
extern QErrorCode OpenRSPSetOverlap(OpenRSP*,
                                    const QInt,
                                    const QcPertInt*,
                                    const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                    void*,
#endif
                                    const GetOverlapMat,
                                    const GetOverlapExp);
extern QErrorCode OpenRSPAddOneOper(OpenRSP*,
                                    const QInt,
                                    const QcPertInt*,
                                    const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                    void*,
#endif
                                    const GetOneOperMat,
                                    const GetOneOperExp);
extern QErrorCode OpenRSPAddTwoOper(OpenRSP*,
                                    const QInt,
                                    const QcPertInt*,
                                    const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                    void*,
#endif
                                    const GetTwoOperMat,
                                    const GetTwoOperExp);
extern QErrorCode OpenRSPAddXCFun(OpenRSP*,
                                  const QInt,
                                  const QcPertInt*,
                                  const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                  void*,
#endif
                                  const GetXCFunMat,
                                  const GetXCFunExp);
extern QErrorCode OpenRSPAddZeroOper(OpenRSP*,
                                     const QInt,
                                     const QcPertInt*,
                                     const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                     void*,
#endif 
                                     const GetZeroOperContrib);
extern QErrorCode OpenRSPSetLinearRSPSolver(OpenRSP*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                            void*,
#endif
                                            const GetLinearRSPSolution);
extern QErrorCode OpenRSPAssemble(OpenRSP*);
extern QErrorCode OpenRSPWrite(const OpenRSP*,FILE*);
extern QErrorCode OpenRSPGetRSPFun(OpenRSP*,
                                   const QcMat*,
                                   const QcMat*,
                                   const QcMat*,
                                   const QInt,
                                   const QInt*,
                                   const QcPertInt*,
                                   const QInt*,
                                   const QReal*,
                                   const QInt*,
                                   const QInt,
                                   const QReal,
                                   const QInt,
                                   QReal*);
extern QErrorCode OpenRSPGetResidue(OpenRSP*,
                                    const QcMat*,
                                    const QcMat*,
                                    const QcMat*,
                                    const QInt,
                                    const QInt,
                                    const QReal*,
                                    QcMat*[],
                                    const QInt,
                                    const QInt*,
                                    const QcPertInt*,
                                    const QInt*,
                                    const QInt*,
                                    const QInt*,
                                    const QReal*,
                                    const QInt*,
                                    const QInt,
                                    const QReal,
                                    const QInt,
                                    QReal*);
extern QErrorCode OpenRSPDestroy(OpenRSP*);


#endif

