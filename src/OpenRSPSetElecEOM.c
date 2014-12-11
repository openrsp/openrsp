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

   This file implements the function OpenRSPSetElecEOM().

   2014-08-04, Bin Gao:
   * first version
*/

#include "openrsp.h"

/*@% \brief sets the equation of motion (EOM) of electrons
     \author Bin Gao
     \date 2014-08-04
     \param[OneRSP:struct]{inout} open_rsp the context of response theory calculations
     \param[ElecEOM:enum]{in} elec_EOM_type the type of EOM of electrons
     \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPSetElecEOM(OpenRSP *open_rsp, const ElecEOMType elec_EOM_type)
{
    switch (elec_EOM_type) {
    /* density matrix-based response theory */
    case ELEC_AO_D_MATRIX:
        break;
    /* molecular orbital (MO) coefficient matrix-based response theory */
    case ELEC_MO_C_MATRIX:
        break;
    /* couple cluster-based response theory */
    case ELEC_COUPLED_CLUSTER:
        break;
    default:
        printf("OpenRSPSetElecEOM>> type of EOM of electrons %d\n", elec_EOM_type);
        QErrorExit(FILE_AND_LINE, "invalid type of EOM of electrons");
    }
    open_rsp->elec_EOM_type = elec_EOM_type;
    return QSUCCESS;
}
