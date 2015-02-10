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

   This is the header file of callback functions for the C tests.

   2014-07-31, Bin Gao:
   * first version
*/

#if !defined(OPENRSP_C_PERT_CALLBACK_H)
#define OPENRSP_C_PERT_CALLBACK_H

/* QcMatrix library */
#include "qcmatrix.h"

/* callback function for getting components of a perturbation */
extern QVoid get_pert_comp(const QInt,
                           const QInt,
                           const QInt,
#if defined(OPENRSP_C_USER_CONTEXT)
                           QVoid*,
#endif
                           QInt*,
                           QInt*,
                           QInt*);

/* callback function for getting rank of a perturbation */
extern QVoid get_pert_rank(const QInt,
                           const QInt,
                           const QInt*,
                           const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                           QVoid*,
#endif
                           QInt*);

#endif
