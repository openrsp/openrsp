/*
  OpenRSP: open-ended library for response theory
  Copyright 2015 Radovan Bast,
                 Daniel H. Friese,
                 Bin Gao,
                 Dan J. Jonsson,
                 Magnus Ringholm,
                 Kenneth Ruud,
                 Andreas Thorvaldsen

  OpenRSP is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.

  OpenRSP is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with OpenRSP. If not, see <http://www.gnu.org/licenses/>.

*/

#include "RSPXCFun.h"

/*% \brief destroys the linked list of XC functionals, should be called at the end
    \author Bin Gao
    \date 2015-06-23
    \param[RSPXCFun:struct]{inout} xc_fun the linked list of XC functionals
    \return[QErrorCode:int] error information
*/
QErrorCode RSPXCFunDestroy(RSPXCFun **xc_fun)
{
    RSPXCFun *cur_xc;   /* current XC functional */
    RSPXCFun *next_xc;  /* next XC functional */
    /* walks to the last XC functional */
    cur_xc = *xc_fun;
    while (cur_xc!=NULL) {
        cur_xc->num_pert = 0;
        free(cur_xc->pert_labels);
        cur_xc->pert_labels = NULL;
        free(cur_xc->pert_max_orders);
        cur_xc->pert_max_orders = NULL;
        cur_xc->user_ctx = NULL;
        cur_xc->get_xc_fun_mat = NULL;
        cur_xc->get_xc_fun_exp = NULL;
        next_xc = cur_xc->next_xc;
        free(cur_xc);
        cur_xc = NULL;
        cur_xc = next_xc;
    }
    return QSUCCESS;
}

