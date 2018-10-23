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


  This is the header file of perturbations' callback function.

  2014-07-31, Bin Gao:
  * first version
*/

#if !defined(OPENRSP_PERT_CALLBACK_H)
#define OPENRSP_PERT_CALLBACK_H

#include "OpenRSP.h"

#if defined(OPENRSP_C_USER_CONTEXT)
#define PERT_CONTEXT "NRNZGEO"
#endif

extern void get_pert_concatenation(const QInt,
                                   const QcPertInt,
                                   const QInt,
                                   const QInt,
                                   const QInt*,
#if defined(OPENRSP_C_USER_CONTEXT)
                                   void*,
#endif
                                   QInt*);

#endif
