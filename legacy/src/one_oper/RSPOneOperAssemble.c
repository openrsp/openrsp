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

   This file implements the function RSPOneOperAssemble().

   2014-07-30, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_one_oper.h"

/*% \brief assembles the linked list of one-electron operators
    \author Bin Gao
    \date 2014-07-30
    \param[RSPOneOper:struct]{inout} one_oper the linked list of one-electron operators
    \return[QErrorCode:int] error information
*/
QErrorCode RSPOneOperAssemble(RSPOneOper *one_oper)
{
    QInt ioper;            /* incremental recorder over opertors */
    RSPOneOper *cur_oper;  /* current operator */
    /* walks to the last operator */
    ioper = 0;
    cur_oper = one_oper;
    do {
        /*FIXME: to implement */
        ioper++;
        cur_oper = cur_oper->next_oper;
    } while (cur_oper!=NULL);
    return QSUCCESS;
}
