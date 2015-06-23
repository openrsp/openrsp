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
