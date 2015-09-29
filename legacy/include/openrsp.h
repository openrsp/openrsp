/* OpenRSP: open-ended library for response theory
   Copyright 2015 Radovan Bast,
                  Daniel H. Friese,
                  Bin Gao,
                  Dan J. Jonsson,
                  Magnus Ringholm,
                  Kenneth Ruud,
                  Andreas Thorvaldsen

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

   This is the header file of application program interface of OpenRSP.

   2014-01-27, Bin Gao:
   * first version
*/

#if !defined(OPENRSP_H)
#define OPENRSP_H

/* perturbations involved in calculations */
#include "eom/openrsp_perturbation.h"
/* electronic wave function */
#include "eom/openrsp_wav_fun.h"
/* linear response equation solver */
#include "eom/openrsp_solver.h"
/* overlap integrals due to perturbation dependent basis sets */
#include "hamiltonian/rsp_overlap.h"
/* one-electron operators */
#include "hamiltonian/rsp_one_oper.h"
/* two-electron operators */
#include "hamiltonian/rsp_two_oper.h"
/* exchange-correlation functionals */
#include "hamiltonian/rsp_xc_fun.h"
/* nuclear Hamiltonian contributions */
#include "hamiltonian/rsp_nuc_contrib.h"

/* context of response theory calculations */
typedef struct {
    QBool assembled;               /* indicates if the context of response theory calculations assembled */
    /*ElecWav *elec_wav;*/           /* implementation-specific data of (electronic) wave function */
    ElecWavType elec_wav_type;
    RSPPert *rsp_pert;             /* perturbations */
    RSPSolver *rsp_solver;         /* linear response equation solver */
    RSPOverlap *overlap;           /* overlap integrals */
    RSPOneOper *one_oper;          /* linked list of one-electron operators */
    RSPTwoOper *two_oper;          /* linked list of two-electron operators */
    RSPXCFun *xc_fun;              /* linked list of exchange-correlation functionals */
    RSPNucHamilton *nuc_hamilton;  /* nuclear Hamiltonian */
} OpenRSP;

/* APIs of OpenRSP */
extern QErrorCode OpenRSPCreate(OpenRSP*);
extern QErrorCode OpenRSPSetWaveFunction(OpenRSP*,const ElecWavType);
extern QErrorCode OpenRSPSetPerturbations(OpenRSP*,
                                          const QInt,
                                          const QInt*,
                                          const QInt*,
                                          const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                          QVoid*,
#endif
                                          const GetPertCat);
extern QErrorCode OpenRSPSetLinearRSPSolver(OpenRSP*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                            QVoid*,
#endif
                                            const GetLinearRSPSolution);
extern QErrorCode OpenRSPSetPDBS(OpenRSP*,
                                 const QInt,
                                 const QInt*,
                                 const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                 QVoid*,
#endif
                                 const GetOverlapMat,
                                 const GetOverlapExp);
extern QErrorCode OpenRSPAddOneOper(OpenRSP*,
                                    const QInt,
                                    const QInt*,
                                    const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                    QVoid*,
#endif
                                    const GetOneOperMat,
                                    const GetOneOperExp);
extern QErrorCode OpenRSPAddTwoOper(OpenRSP*,
                                    const QInt,
                                    const QInt*,
                                    const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                    QVoid*,
#endif
                                    const GetTwoOperMat,
                                    const GetTwoOperExp);
extern QErrorCode OpenRSPAddXCFun(OpenRSP*,
                                  const QInt,
                                  const QInt*,
                                  const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                  QVoid*,
#endif
                                  const GetXCFunMat,
                                  const GetXCFunExp);
extern QErrorCode OpenRSPSetNucContributions(OpenRSP*,
                                             const QInt,
                                             const QInt*,
                                             const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                             QVoid*,
#endif 
                                             const GetNucContrib,
/*FIXME: num_atoms to be removed after perturbation free scheme implemented*/
                                             const QInt);
extern QErrorCode OpenRSPAssemble(OpenRSP*);
extern QErrorCode OpenRSPWrite(const OpenRSP*,const QChar*);
extern QErrorCode OpenRSPGetRSPFun(OpenRSP*,
                                   const QcMat*,
                                   const QcMat*,
                                   const QcMat*,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
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
                                    const QInt*,
                                    const QInt*,
                                    const QInt*,
                                    const QInt*,
                                    const QReal*,
                                    const QInt*,
                                    const QInt,
                                    QReal*);
extern QErrorCode OpenRSPDestroy(OpenRSP*);

#endif
