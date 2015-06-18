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

   This file implements the function RSPNucHamiltonWrite().

   2015-02-12, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_nuc_contrib.h"

/*% \brief writes the context of nuclear Hamiltonian
    \author Bin Gao
    \date 2015-02-12
    \param[RSPNucHamilton:struct]{in} nuc_hamilton the context of nuclear Hamiltonian
    \param[FILE]{inout} fp_nuc file pointer
    \return[QErrorCode:int] error information
*/
QErrorCode RSPNucHamiltonWrite(const RSPNucHamilton *nuc_hamilton, FILE *fp_nuc)
{
    QInt iatom,ixyz;
    fprintf(fp_nuc,
            "RSPNucHamiltonWrite>> number of atoms %"QINT_FMT"\n",
            nuc_hamilton->num_atoms);
    fprintf(fp_nuc,
            "RSPNucHamiltonWrite>> atom    charge    coordinates\n");
    for (iatom=0,ixyz=0; iatom<nuc_hamilton->num_atoms; iatom++) {
        fprintf(fp_nuc,
                "RSPNucHamiltonWrite>> %"QINT_FMT"    %f    [%f, %f, %f]\n",
                iatom,
                nuc_hamilton->atom_charge[iatom],
                nuc_hamilton->atom_coord[ixyz],     /* x */
                nuc_hamilton->atom_coord[ixyz+1],   /* y */
                nuc_hamilton->atom_coord[ixyz+2]);  /* z */
       ixyz += 3;
    }
    if (nuc_hamilton->dipole_origin!=NULL) {
        fprintf(fp_nuc,
                "RSPNucHamiltonWrite>> dipole origin [%f, %f, %f]\n",
                nuc_hamilton->dipole_origin[0],
                nuc_hamilton->dipole_origin[1],
                nuc_hamilton->dipole_origin[2]);
    }
    if (nuc_hamilton->gauge_origin!=NULL) {
        fprintf(fp_nuc,
                "RSPNucHamiltonWrite>> gauge origin [%f, %f, %f]\n",
                nuc_hamilton->gauge_origin[0],
                nuc_hamilton->gauge_origin[1],
                nuc_hamilton->gauge_origin[2]);
    }
    return QSUCCESS;
}
