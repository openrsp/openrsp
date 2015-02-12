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

   This file implements the function RSPNucHamiltonSetScalarPotential().

   2015-02-12, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_nuc_hamiltonian.h"

/*% \brief sets the terms in nuclear Hamiltonian due to the scalar potential
    \author Bin Gao
    \date 2015-02-12
    \param[RSPNucHamilton:struct]{inout} nuc_hamilton the context of nuclear Hamiltonian
    \param[QReal:real]{in} dipole_origin coordinates of dipole origin
    \return[QErrorCode:int] error information
*/
QErrorCode RSPNucHamiltonSetScalarPotential(RSPNucHamilton *nuc_hamilton,
                                            const QReal dipole_origin[3])
{
    if (nuc_hamilton->dipole_origin==NULL) {
        nuc_hamilton->dipole_origin = (QReal *)malloc(3*sizeof(QReal));
        if (nuc_hamilton->dipole_origin==NULL) {
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for dipole_origin");
        }
    }
    nuc_hamilton->dipole_origin[0] = dipole_origin[0];
    nuc_hamilton->dipole_origin[1] = dipole_origin[1];
    nuc_hamilton->dipole_origin[2] = dipole_origin[2];
    return QSUCCESS;
}
