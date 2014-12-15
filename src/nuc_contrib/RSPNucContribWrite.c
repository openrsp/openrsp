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

   This file implements the function RSPNucContribWrite().

   2014-12-15, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_nuc_contrib.h"

/*% \brief writes the context of nuclear contributions
    \author Bin Gao
    \date 2014-12-15
    \param[RSPNucContrib:struct]{in} nuc_contrib the context of nuclear contributions
    \param[FILE]{inout} fp_nuc file pointer
    \return[QErrorCode:int] error information
*/
QErrorCode RSPNucContribWrite(const RSPNucContrib *nuc_contrib, FILE *fp_nuc)
{
    QInt iatom,ixyz;
    fprintf(fp_nuc,
            "RSPNucContribWrite>> number of atoms %"QINT_FMT"\n",
            nuc_contrib->num_atoms);
    fprintf(fp_nuc,
            "RSPNucContribWrite>> atom    charge    coordinates\n");
    for (iatom=0,ixyz=0; iatom<num_atoms; iatom++) {
        fprintf(fp_nuc,
                "RSPNucContribWrite>> %"QINT_FMT"    %f    [%f, %f, %f]\n",
                iatom,
                nuc_contrib->atom_charge[iatom],
                nuc_contrib->atom_coord[ixyz++],   /* x */
                nuc_contrib->atom_coord[ixyz++],   /* y */
                nuc_contrib->atom_coord[ixyz++]);  /* z */
    }
    fprintf(fp_nuc,
            "RSPNucContribWrite>> dipole origin [%f, %f, %f]\n",
            nuc_contrib->dipole_origin[0],
            nuc_contrib->dipole_origin[1],
            nuc_contrib->dipole_origin[2]);
    fprintf(fp_nuc,
            "RSPNucContribWrite>> gauge origin [%f, %f, %f]\n",
            nuc_contrib->gauge_origin[0],
            nuc_contrib->gauge_origin[1],
            nuc_contrib->gauge_origin[2]);
    return QSUCCESS;
}
