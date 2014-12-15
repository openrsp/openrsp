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

   2014-12-15, Bin Gao:
   * first version
*/

#include "openrsp.h"

/*% \brief sets the coordinates of dipole origin
    \author Bin Gao
    \date 2014-12-15
    \param[OneRSP:struct]{inout} open_rsp the context of response theory calculations
    \param[QReal:real]{in} dipole_origin coordinates of dipole origin
    \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPSetDipoleOrigin(OpenRSP *open_rsp,
                                  const QReal dipole_origin[3])
{
    QErrorCode ierr;  /* error information */
    if (open_rsp->nuc_contrib==NULL) {
        QErrorExit(FILE_AND_LINE, "OpenRSPSetAtoms() should be called at first");
    }
    ierr = RSPNucContribSetDipoleOrigin(open_rsp->nuc_contrib, dipole_origin);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPNucContribSetDipoleOrigin");
    return QSUCCESS;
}
