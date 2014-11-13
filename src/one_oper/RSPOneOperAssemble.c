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
