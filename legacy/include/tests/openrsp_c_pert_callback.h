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
