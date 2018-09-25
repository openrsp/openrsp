/* OpenRSP: open-ended library for response theory
   Copyright 2014
  This source code form is subject to the terms of the
  GNU Lesser General Public License, version 2.1.
  If a copy of the GNU LGPL v2.1 was not distributed with this
  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

   This file initializes molecule H2O2.

   2015-02-10, Bin Gao:
   * first version
*/

#include "tests/openrsp_c_molecule.h"

const QReal ATOM_COORD[] = {0.00000000,1.40784586,-0.09885600,
                            0.00000000,-1.40784586,-0.09885600,
                            0.69081489,1.72614891,1.56881868,
                            -0.69081489,-1.72614891,1.56881868};
const QReal ATOM_CHARGE[] = {8.0, 8.0, 1.0, 1.0};
const QReal DIPOLE_ORIGIN[] = {0.0,0.0,0.0};
const QReal GAUGE_ORIGIN[] = {0.0,0.0,0.0};
