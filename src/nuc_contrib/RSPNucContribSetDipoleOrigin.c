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

   This file implements the function RSPNucContribSetDipoleOrigin().

   2014-12-15, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_nuc_contrib.h"

/*% \brief sets the coordinates of dipole origin
    \author Bin Gao
    \date 2014-12-15
    \param[RSPNucContrib:struct]{inout} nuc_contrib the context of nuclear contributions
    \param[QReal:real]{in} dipole_origin coordinates of dipole origin
    \return[QErrorCode:int] error information
*/
QErrorCode RSPNucContribSetDipoleOrigin(RSPNucContrib *nuc_contrib,
                                        const QReal dipole_origin[3])
{
    nuc_contrib->dipole_origin[0] = dipole_origin[0];
    nuc_contrib->dipole_origin[1] = dipole_origin[1];
    nuc_contrib->dipole_origin[2] = dipole_origin[2];
    return QSUCCESS;
}
