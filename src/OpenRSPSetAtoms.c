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

   This file implements the function OpenRSPSetAtoms().

   2014-12-15, Bin Gao:
   * first version
*/

#include "openrsp.h"

/*% \brief sets the context of atoms for the nuclear contributions
    \author Bin Gao
    \date 2014-12-15
    \param[OneRSP:struct]{inout} open_rsp the context of response theory calculations
    \param[QInt:int]{in} num_atoms number of atoms
    \param[QReal:real]{in} atom_coord coordinates of atoms
    \param[QReal:real]{in} atom_charge charges of atoms
    \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPSetAtoms(OpenRSP *open_rsp,
                           const QInt num_atoms,
                           const QReal *atom_coord,
                           const QReal *atom_charge)
{
    QErrorCode ierr;  /* error information */
    /* creates the context of nuclear contributions */
    if (open_rsp->nuc_contrib!=NULL) {
        ierr = RSPNucContribDestroy(open_rsp->nuc_contrib);
        QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPNucContribDestroy");
    }
    else {
        open_rsp->nuc_contrib = (RSPNucContrib *)malloc(sizeof(RSPNucContrib));
        if (open_rsp->nuc_contrib==NULL) {
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for nuc_contrib");
        }
    }
    ierr = RSPNucContribCreate(open_rsp->nuc_contrib,
                               num_atoms,
                               atom_coord,
                               atom_charge);
    QErrorCheckCode(ierr, FILE_AND_LINE, "calling RSPNucContribCreate");
    return QSUCCESS;
}
