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

   This file implements the function RSPNucHamiltonDestroy().

   2015-02-12, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_nuc_hamiltonian.h"

/*% \brief destroys the context of nuclear Hamiltonian, should be called at the end
    \author Bin Gao
    \date 2015-02-12
    \param[RSPNucHamilton:struct]{inout} nuc_hamilton the context of nuclear Hamiltonian
    \return[QErrorCode:int] error information
*/
QErrorCode RSPNucHamiltonDestroy(RSPNucHamilton *nuc_hamilton)
{
    if (nuc_hamilton->atom_coord!=NULL) {
        free(nuc_hamilton->atom_coord);
        nuc_hamilton->atom_coord = NULL;
    }
    if (nuc_hamilton->atom_charge!=NULL) {
        free(nuc_hamilton->atom_charge);
        nuc_hamilton->atom_charge = NULL;
    }
    if (nuc_hamilton->dipole_origin!=NULL) {
        free(nuc_hamilton->dipole_origin);
        nuc_hamilton->dipole_origin = NULL;
    }
    if (nuc_hamilton->gauge_origin!=NULL) {
        free(nuc_hamilton->gauge_origin);
        nuc_hamilton->gauge_origin = NULL;
    }
    return QSUCCESS;
}
