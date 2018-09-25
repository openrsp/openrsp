/* OpenRSP: open-ended library for response theory
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
extern QErrorCode RSPPertGetFromTuple(const RSPPert*,
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
