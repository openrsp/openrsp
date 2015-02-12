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

   This file implements the function OpenRSPSetDipoleOrigin().

   2015-02-12, Bin Gao:
   * first version
*/

#include "openrsp.h"

/*% \brief sets the terms in nuclear Hamiltonian due to the scalar potential
    \author Bin Gao
    \date 2015-02-12
    \param[OneRSP:struct]{inout} open_rsp the context of response theory calculations
    \param[QReal:real]{in} dipole_origin coordinates of dipole origin
    \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPSetNucScalarPotential(OpenRSP *open_rsp,
                                        const QReal dipole_origin[3])
{
    QErrorCode ierr;  /* error information */
    /* creates the context of nuclear Hamiltonian */
    if (open_rsp->nuc_hamilton==NULL) {
        open_rsp->nuc_hamilton = (RSPNucHamilton *)malloc(sizeof(RSPNucHamilton));
        if (open_rsp->nuc_hamilton==NULL) {
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for nuc_hamilton");
        }
        ierr = RSPNucHamiltonCreate(open_rsp->nuc_hamilton);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPNucHamiltonCreate");
    }
    ierr = RSPNucHamiltonSetScalarPotential(open_rsp->nuc_hamilton,
                                            dipole_origin);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPNucHamiltonSetScalarPotential");
    return QSUCCESS;
}
