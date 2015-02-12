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

   This file implements the function RSPNucHamiltonSetGeoPerturbations().

   2015-02-12, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_nuc_hamiltonian.h"

/*% \brief sets the context of geometric perturbations for nuclear Hamiltonian
    \author Bin Gao
    \date 2015-02-12
    \param[RSPNucHamilton:struct]{inout} nuc_hamilton the context of nuclear Hamiltonian
    \param[QInt:int]{in} num_atoms number of atoms
    \param[QReal:real]{in} atom_coord coordinates of atoms
    \param[QReal:real]{in} atom_charge charges of atoms
    \return[QErrorCode:int] error information
*/
QErrorCode RSPNucHamiltonSetGeoPerturbations(RSPNucHamilton *nuc_hamilton,
                                             const QInt num_atoms,
                                             const QReal *atom_coord,
                                             const QReal *atom_charge)
{
    QInt iatom,ixyz;
    if (num_atoms>0) {
        nuc_hamilton->num_atoms = num_atoms;
    }
    else {
        printf("RSPNucHamiltonSetGeoPerturbations>> number of atoms %"QINT_FMT"\n",
               num_atoms);
        QErrorExit(FILE_AND_LINE, "invalid number of atoms");
    }
    if (nuc_hamilton->atom_coord!=NULL) {
        free(nuc_hamilton->atom_coord);
    }
    nuc_hamilton->atom_coord = (QReal *)malloc(3*num_atoms*sizeof(QReal));
    if (nuc_hamilton->atom_coord==NULL) {
        printf("RSPNucHamiltonSetGeoPerturbations>> number of atoms %"QINT_FMT"\n",
               num_atoms);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for atom_coord");
    }
    if (nuc_hamilton->atom_charge!=NULL) {
        free(nuc_hamilton->atom_charge);
    }
    nuc_hamilton->atom_charge = (QReal *)malloc(num_atoms*sizeof(QReal));
    if (nuc_hamilton->atom_charge==NULL) {
        printf("RSPNucHamiltonSetGeoPerturbations>> number of atoms %"QINT_FMT"\n",
               num_atoms);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for atom_charge");
    }
    for (iatom=0,ixyz=0; iatom<num_atoms; iatom++) {
        nuc_hamilton->atom_charge[iatom] = atom_charge[iatom];
        nuc_hamilton->atom_coord[ixyz] = atom_coord[ixyz];  /* x */
        ixyz++;
        nuc_hamilton->atom_coord[ixyz] = atom_coord[ixyz];  /* y */
        ixyz++;
        nuc_hamilton->atom_coord[ixyz] = atom_coord[ixyz];  /* z */
        ixyz++;
    }
    return QSUCCESS;
}
