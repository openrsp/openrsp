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


  This is the header file of OpenRSP unit testing.

  2014-07-31, Bin Gao:
  * first version
*/

#if !defined(OPENRSP_TEST_H)
#define OPENRSP_TEST_H

#include <string.h>
#include "OpenRSP.h"

#include "OpenRSPTestPerturbations.h"
#include "OpenRSPPertCallback.h"
#include "OpenRSPZeroOperCallback.h"
#include "OpenRSPDMatTest.h"

#if !defined(OPENRSP_TEST_EXECUTABLE)
extern QErrorCode OpenRSPTest(FILE *fp_log);
#endif

#endif

