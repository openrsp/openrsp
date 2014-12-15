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

   This file implements the function RSPNucContribCreate().

   2014-12-15, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_nuc_contrib.h"

/*% \brief creates the context of nuclear contributions, should be called at first
    \author Bin Gao
    \date 2014-12-15
    \param[RSPNucContrib:struct]{inout} nuc_contrib the context of nuclear contributions
    \param[QInt:int]{in} num_atoms number of atoms
    \param[QReal:real]{in} atom_coord coordinates of atoms
    \param[QReal:real]{in} atom_charge charges of atoms
    \return[QErrorCode:int] error information
*/
QErrorCode RSPNucContribCreate(RSPNucContrib *nuc_contrib,
                               const QInt num_atoms,
                               const QInt *atom_coord,
                               const QInt *atom_charge)
{
    QInt iatom,ixyz;
    if (num_atoms>0) {
        nuc_contrib->num_atoms = num_atoms;
    }
    else {
        printf("RSPNucContribCreate>> number of atoms %"QINT_FMT"\n", num_atoms);
        QErrorExit(FILE_AND_LINE, "invalid number of atoms");
    }
    nuc_contrib->atom_coord = (QReal *)malloc(3*num_atoms*sizeof(QReal));
    if (nuc_contrib->atom_coord==NULL) {
        printf("RSPNucContribCreate>> number of atoms %"QINT_FMT"\n", num_atoms);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for atom_coord");
    }
    nuc_contrib->atom_charge = (QReal *)malloc(num_atoms*sizeof(QReal));
    if (nuc_contrib->atom_charge==NULL) {
        printf("RSPNucContribCreate>> number of atoms %"QINT_FMT"\n", num_atoms);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for atom_charge");
    }
    for (iatom=0,ixyz=0; iatom<num_atoms; iatom++) {
        nuc_contrib->atom_charge[iatom] = atom_charge[iatom];
        nuc_contrib->atom_coord[ixyz] = atom_coord[ixyz];  /* x */
        ixyz++;
        nuc_contrib->atom_coord[ixyz] = atom_coord[ixyz];  /* y */
        ixyz++;
        nuc_contrib->atom_coord[ixyz] = atom_coord[ixyz];  /* z */
        ixyz++;
    }
    nuc_contrib->dipole_origin[0] = 0.0;
    nuc_contrib->dipole_origin[1] = 0.0;
    nuc_contrib->dipole_origin[2] = 0.0;
    nuc_contrib->gauge_origin[0] = 0.0;
    nuc_contrib->gauge_origin[1] = 0.0;
    nuc_contrib->gauge_origin[2] = 0.0;
    return QSUCCESS;
}
