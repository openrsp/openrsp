/*
  OpenRSP: open-ended library for response theory
  Copyright 2015 Radovan Bast,
                 Daniel H. Friese,
                 Bin Gao,
                 Dan J. Jonsson,
                 Magnus Ringholm,
                 Kenneth Ruud

  This source code form is subject to the terms of the
  GNU Lesser General Public License, version 2.1.
  If a copy of the GNU LGPL v2.1 was not distributed with this
  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.


   This is the header file of unit testing of the AO density matrix-based
   response theory.

   2014-07-31, Bin Gao
   * first version
*/
#if !defined(OPENRSP_DMAT_TEST_H)
#define OPENRSP_DMAT_TEST_H

#include "OpenRSPPertCallback.h"
#include "OpenRSPDMatCallback.h"

extern QErrorCode OpenRSPDMatTest(OpenRSP*,FILE*);

#endif
