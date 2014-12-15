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

   This is the header file of application program interface of OpenRSP.

   2014-01-27, Bin Gao:
   * first version
*/

#if !defined(OPENRSP_H)
#define OPENRSP_H

/* type of equation of motion (EOM) of electrons */
#include "eom/openrsp_elec_eom.h"
/* response equation solver */
#include "eom/openrsp_solver.h"
/* overlap integrals due to perturbation dependent basis sets */
#include "hamiltonian/rsp_overlap.h"
/* one-electron operators */
#include "hamiltonian/rsp_one_oper.h"
/* two-electron operators */
#include "hamiltonian/rsp_two_oper.h"
/* exchange-correlation functionals */
#include "hamiltonian/rsp_xc_fun.h"
/* (derivatives of) nuclear repulsion and nuclei-field interaction */
#include "hamiltonian/rsp_nuc_contrib.h"

#if defined(OPENRSP_PERTURBATION_FREE)
/* callback function to get the components of a perturbation */
typedef QVoid (*GetPertComp)(const QInt,
                             const QInt,
                             const QInt,
#if defined(OPENRSP_C_USER_CONTEXT)
                             QVoid*,
#endif
                             QInt*,
                             QInt*,
                             QInt*);
/* callback function to get the rank of a perturbation */
typedef QVoid (*GetPertRank)(const QInt,
                             const QInt,
                             const QInt*,
                             const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                             QVoid*,
#endif
                             QInt*);
#endif

/* context of response theory calculations */
typedef struct {
    QBool assembled;             /* indicates if the context of response theory calculations assembled */
#if defined(OPENRSP_PERTURBATION_FREE)
    /* perturbations */
    QInt num_pert;               /* number of all perturbations involved in calculations */
    QInt *perturbations;         /* all perturbations involved in calculations */
    QInt *pert_max_orders;       /* maximum allowed orders of all perturbations */
    QInt *size_ptr;              /* pointer to the size of each perturbation */
    QInt *pert_sizes;            /* sizes of all perturbations up to their maximum orders */
#if defined(OPENRSP_C_USER_CONTEXT)
    QVoid *user_ctx;             /* user-defined callback function context */
#endif
    GetPertComp get_pert_comp;   /* user specified function for getting components of a perturbation */
    GetPertRank get_pert_rank;   /* user specified function for getting rank of a perturbation */
#endif
    /* EOM and solver */
    //ElecEOM *elec_eom;           /* implementation-specific data of the EOM of electrons */
    ElecEOMType elec_EOM_type;
    RSPSolver *rsp_solver;       /* response equation solver */
    /* Hamiltonian */
    RSPOverlap *overlap;         /* overlap integrals */
    RSPOneOper *one_oper;        /* linked list of one-electron operators */
    RSPTwoOper *two_oper;        /* linked list of two-electron operators */
    RSPXCFun *xc_fun;            /* linked list of exchange-correlation functionals */
    RSPNucContrib *nuc_contrib;  /* (derivatives of) nuclear repulsion and nuclei-field interaction */
} OpenRSP;

/* APIs of OpenRSP */
extern QErrorCode OpenRSPCreate(OpenRSP*);
extern QErrorCode OpenRSPSetElecEOM(OpenRSP*,const ElecEOMType);
extern QErrorCode OpenRSPSetSolver(OpenRSP*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   QVoid*,
#endif
                                   const GetRSPSolution);
#if defined(OPENRSP_PERTURBATION_FREE)
extern QErrorCode OpenRSPSetPerturbations(OpenRSP*,
                                          const QInt,
                                          const QInt*,
                                          const QInt*,
                                          const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                          QVoid*,
#endif
                                          const GetPertComp,
                                          const GetPertRank);
#endif
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
//extern QErrorCode OpenRSPAddXCFun(OpenRSP*,);
extern QErrorCode OpenRSPSetAtoms(OpenRSP*,
                                  const QInt,
                                  const QReal*,
                                  const QReal*);
extern QErrorCode OpenRSPSetDipoleOrigin(OpenRSP*,const QReal[3]);
extern QErrorCode OpenRSPSetGaugeOrigin(OpenRSP*,const QReal[3]);
extern QErrorCode OpenRSPAssemble(OpenRSP*);
extern QErrorCode OpenRSPWrite(const OpenRSP*,const QChar*);
extern QErrorCode OpenRSPGetRSPFun(OpenRSP*,
                                   const QMat*,
                                   const QMat*,
                                   const QMat*,
                                   const QInt,
                                   const QInt*,
                                   const QInt*,
                                   const QReal*,
                                   const QInt[],
                                   const QInt,
                                   QReal*);
//extern QErrorCode OpenRSPGetResidue(OpenRSP*,);
extern QErrorCode OpenRSPDestroy(OpenRSP*);

#endif
