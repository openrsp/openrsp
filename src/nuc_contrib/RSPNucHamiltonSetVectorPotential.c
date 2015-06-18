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

   This file implements the function RSPNucHamiltonSetVectorPotential().

   2015-02-12, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_nuc_contrib.h"

/*% \brief sets the terms in nuclear Hamiltonian due to the vector potential
    \author Bin Gao
    \date 2015-02-12
    \param[RSPNucHamilton:struct]{inout} nuc_hamilton the context of nuclear Hamiltonian
    \param[QReal:real]{in} gauge_origin coordinates of gauge origin
    \return[QErrorCode:int] error information
*/
QErrorCode RSPNucHamiltonSetVectorPotential(RSPNucHamilton *nuc_hamilton,
                                            const QReal gauge_origin[3])
{
    if (nuc_hamilton->gauge_origin==NULL) {
        nuc_hamilton->gauge_origin = (QReal *)malloc(3*sizeof(QReal));
        if (nuc_hamilton->gauge_origin==NULL) {
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for gauge_origin");
        }
    }
    nuc_hamilton->gauge_origin[0] = gauge_origin[0];
    nuc_hamilton->gauge_origin[1] = gauge_origin[1];
    nuc_hamilton->gauge_origin[2] = gauge_origin[2];
    return QSUCCESS;
}
