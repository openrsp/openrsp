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
