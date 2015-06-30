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

   This is the header file of perturbations.

   2015-06-23, Bin Gao:
   * first version
*/

#if !defined(OPENRSP_PERTURBATION_H)
#define OPENRSP_PERTURBATION_H

/* QcMatrix library */
#include "qcmatrix.h"

/* callback function to get the ranks of components of sub-perturbation
   tuples (with same perturbation label) for given components of the
   corresponding concatenated perturbation tuple */
typedef QVoid (*GetPertCat)(const QInt,
                            const QInt,
                            const QInt*,
                            const QInt,
                            const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                            QVoid*,
#endif
                            QInt*);

/* context of all perturbations involved in calculations */
typedef struct {
    QInt num_pert;                      /* number of all different perturbation labels involved */
    QInt *pert_labels;                  /* all different perturbation labels involved */
    QInt *pert_max_orders;              /* maximum allowed order of each perturbation (label) */
    QInt *ptr_ncomp;                    /* pointer to the numbers of components of each perturbation */
    QInt *pert_num_comps;               /* number of components of each perturbation (label) up to
                                           its maximum order */
#if defined(OPENRSP_C_USER_CONTEXT)     
    QVoid *user_ctx;                    /* user-defined callback function context */
#endif
    GetPertCat get_pert_concatenation;  /* user specified function for getting the ranks of
                                           components of sub-perturbation tuples (with same
                                           perturbation label) for given components of the
                                           corresponding concatenated perturbation tuple */
} RSPPert;

/* functions related to the perturbations */
extern QErrorCode RSPPertCreate(RSPPert*,
                                const QInt,
                                const QInt*,
                                const QInt*,
                                const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                QVoid*,
#endif
                                const GetPertCat);
extern QErrorCode RSPPertAssemble(RSPPert*);
extern QErrorCode RSPPertWrite(const RSPPert*,FILE*);
extern QErrorCode RSPPertGetNumComps(const RSPPert*,
                                     const QInt,
                                     const QInt*,
                                     const QInt,
                                     const QReal*,
                                     QInt*);
extern QErrorCode RSPPertGetConcatenation(const RSPPert*,
                                          const QInt,
                                          const QInt,
                                          const QInt,
                                          const QInt,
                                          const QInt*,
                                          QInt*);
extern QErrorCode RSPPertDestroy(RSPPert*);

#endif
