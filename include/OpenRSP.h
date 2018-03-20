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
/* nuclear Hamiltonian */
#include "RSPNucHamilton.h"
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
    RSPNucHamilton *nuc_hamilton;  /* nuclear Hamiltonian */
    RSPSolver *rsp_solver;         /* linear response equation solver */
} OpenRSP;

extern QErrorCode OpenRSPCreate(OpenRSP*);
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
extern QErrorCode OpenRSPSetNucHamilton(OpenRSP*,
                                        const QInt,
                                        const QcPertInt*,
                                        const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                        void*,
#endif 
                                        const GetNucContrib,
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
                                        const QInt);
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
                                    QReal*);
extern QErrorCode OpenRSPDestroy(OpenRSP*);

#endif

