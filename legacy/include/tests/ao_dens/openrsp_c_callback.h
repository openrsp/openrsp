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

   This is the header file of callback functions for the atomic orbital
   density matrix-based response theory.

   2014-07-31, Bin Gao:
   * first version
*/

#if !defined(OPENRSP_C_AO_DENS_CALLBACK_H)
#define OPENRSP_C_AO_DENS_CALLBACK_H

/* QcMatrix library */
#include "qcmatrix.h"

/* callback function to get the solutions of linear response equation */
extern QVoid get_linear_rsp_solution(const QInt,
                                     const QReal*,
                                     const QInt,
                                     QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                                     QVoid*,
#endif
                                     QcMat*[]);

/* callback function for getting integral matrices of overlap integrals */
extern QVoid get_overlap_mat(const QInt,
                             const QInt*,
                             const QInt,
                             const QInt*,
                             const QInt,
                             const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                             QVoid*,
#endif
                             const QInt,
                             QcMat*[]);
/* callback function for getting expectation values of overlap integrals */
extern QVoid get_overlap_exp(const QInt,
                             const QInt*,
                             const QInt,
                             const QInt*,
                             const QInt,
                             const QInt*,
                             const QInt,
                             QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                             QVoid*,
#endif
                             const QInt,
                             QReal*);
/* callback function for getting integral matrices of one-electron operators */
extern QVoid get_one_oper_mat(const QInt,
                              const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                              QVoid*,
#endif
                              const QInt,
                              QcMat*[]);
/* callback function for getting expectation values of one-electron operators */
extern QVoid get_one_oper_exp(const QInt,
                              const QInt*,
                              const QInt,
                              QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                              QVoid*,
#endif
                              const QInt,
                              QReal*);

#endif
/* callback function for getting integral matrices of two-electron operators */
extern QVoid get_two_oper_mat(const QInt,
                              const QInt*,
                              const QInt,
                              QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                              QVoid*,
#endif
                              const QInt,
                              QcMat*[]);
/* callback function for getting expectation values of two-electron operators */
extern QVoid get_two_oper_exp(const QInt,
                              const QInt*,
                              const QInt,
                              const QInt*,
                              QcMat*[],
                              const QInt*,
                              QcMat*[],
#if defined(OPENRSP_C_USER_CONTEXT)
                              QVoid*,
#endif
                              const QInt,
                              QReal*);
