/* OpenRSP: open-ended library for response theory
   Copyright 2015 Radovan Bast,
                  Daniel H. Friese,
                  Bin Gao,
                  Dan J. Jonsson,
                  Magnus Ringholm,
                  Kenneth Ruud,
                  Andreas Thorvaldsen
  This source code form is subject to the terms of the
  GNU Lesser General Public License, version 2.1.
  If a copy of the GNU LGPL v2.1 was not distributed with this
  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

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
