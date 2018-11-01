/*
  OpenRSP: open-ended library for response theory
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


  This header file contains callback function for zero-electron operator
  (nuclear Hamiltonian).

  2015-10-16, Bin Gao:
  * first version
*/

#if !defined(OPENRSP_ZERO_OPER_CALLBACK_H)
#define OPENRSP_ZERO_OPER_CALLBACK_H

#include "OpenRSP.h"

#if defined(OPENRSP_C_USER_CONTEXT)
#define ZERO_OPER_CONTEXT "ZERO_OPER"
#endif

extern void get_zero_oper_contrib(const QInt,
                                  const QcPertInt*,
                                  const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                  void*,
#endif
                                  const QInt,
                                  QReal*);

#endif
