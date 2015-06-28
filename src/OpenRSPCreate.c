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

   This file implements the function OpenRSPCreate().

   2014-01-28, Bin Gao:
   * first version
*/

#include "openrsp.h"

/*@% \brief creates the context of response theory calculations,
         should be called at first
     \author Bin Gao
     \date 2014-01-28
     \param[OpenRSP:struct]{inout} open_rsp the context of response theory calculations
     \return[QErrorCode:int] error information
*/
QErrorCode OpenRSPCreate(OpenRSP *open_rsp)
{
    open_rsp->assembled = QFALSE;
    //open_rsp->elec_wav = NULL;
    open_rsp->elec_wav_type = ELEC_AO_D_MATRIX;
    open_rsp->rsp_pert = NULL;
    open_rsp->rsp_solver = NULL;
    open_rsp->overlap = NULL;
    open_rsp->one_oper = NULL;
    open_rsp->two_oper = NULL;
    open_rsp->xc_fun = NULL;
    open_rsp->nuc_hamilton = NULL;
    return QSUCCESS;
}
