/* OpenRSP: open-ended library for response theory
   Copyright 2015 Radovan Bast,
                  Daniel H. Friese,
                  Bin Gao,
                  Dan J. Jonsson,
                  Magnus Ringholm,
                  Kenneth Ruud,
                  Andreas Thorvaldsen

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

   This file implements the function OpenRSPSetWaveFunction().

   2015-06-29, Bin Gao:
   * first version
*/

#include "openrsp.h"

/*@% \brief sets the type of (electronic) wave function
     \author Bin Gao
     \date 2014-08-04
     \param[OpenRSP:struct]{inout} open_rsp the context of response theory calculations
     \param[ElecWavType:enum]{in} elec_wav_type the type of (electronic) wave function
     \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPSetWaveFunction(OpenRSP *open_rsp, const ElecWavType elec_wav_type)
{
    switch (elec_wav_type) {
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
        printf("OpenRSPSetWaveFunction>> type of (electronic) wave function %d\n",
               elec_wav_type);
        QErrorExit(FILE_AND_LINE, "invalid type of (electronic) wave function");
    }
    open_rsp->elec_wav_type = elec_wav_type;
    return QSUCCESS;
}
