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

   This file implements the function RSPNucContribDestroy().

   2014-12-15, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_nuc_contrib.h"

/*% \brief destroys the context of nuclear contributions, should be called at the end
    \author Bin Gao
    \date 2014-12-15
    \param[RSPNucContrib:struct]{inout} nuc_contrib the context of nuclear contributions
    \return[QErrorCode:int] error information
*/
QErrorCode RSPNucContribDestroy(RSPNucContrib *nuc_contrib)
{
    free(nuc_contrib->atom_coord);
    nuc_contrib->atom_coord = NULL;
    free(nuc_contrib->atom_charge);
    nuc_contrib->atom_charge = NULL;
    return QSUCCESS;
}
