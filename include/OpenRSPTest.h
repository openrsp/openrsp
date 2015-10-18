/*
  OpenRSP: open-ended library for response theory
  Copyright 2015 Radovan Bast,
                 Daniel H. Friese,
                 Bin Gao,
                 Dan J. Jonsson,
                 Magnus Ringholm,
                 Kenneth Ruud,
                 Andreas Thorvaldsen

  OpenRSP is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.

  OpenRSP is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with OpenRSP. If not, see <http://www.gnu.org/licenses/>.


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
#include "OpenRSPNucHamiltonCallback.h"
#include "OpenRSPDensAOTest.h"

#if !defined(OPENRSP_TEST_EXECUTABLE)
extern QErrorCode OpenRSPTest(FILE *fp_log);
#endif

#endif

