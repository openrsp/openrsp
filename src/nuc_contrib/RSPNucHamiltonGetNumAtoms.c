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

   This file implements the function RSPNucHamiltonGetNumAtoms().

   2015-02-12, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_nuc_contrib.h"

/*% \brief gets the number of atoms
    \author Bin Gao
    \date 2015-02-12
    \param[RSPNucHamilton:struct]{in} nuc_hamilton the context of nuclear Hamiltonian
    \param[QInt:int]{out} num_atoms number of atoms
    \return[QErrorCode:int] error information
*/
QErrorCode RSPNucHamiltonGetNumAtoms(const RSPNucHamilton *nuc_hamilton,
                                     QInt *num_atoms)
{
    *num_atoms = nuc_hamilton->num_atoms;
    return QSUCCESS;
}
