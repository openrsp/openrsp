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

   This file implements the function RSPNucContribGet().

   2014-12-15, Bin Gao:
   * first version
*/

#include "hamiltonian/rsp_nuc_contrib.h"

/*% \brief gets the nuclear contributions
    \author Bin Gao
    \date 2014-12-15
    \param[RSPNucContrib:struct]{in} nuc_contrib the context of nuclear contributions
    \param[QInt:int]{in} num_pert number of perturbations
     \param[QInt:int]{in} pert_labels labels of the perturbations
    \param[QInt:int]{in} pert_orders orders of the perturbations
    \param[QInt:int]{in} size_contrib size of nuclear contributions
    \param[QReal:real]{inout} val_contrib the nuclear contributions
    \return[QErrorCode:int] error information
*/
QErrorCode RSPNucContribGet(const RSPNucContrib *nuc_contrib,
                            const QInt num_pert,
                            const QInt *pert_labels,
                            const QInt *pert_orders,
                            const QInt size_contrib,
                            QReal *val_contrib)
{
    return QSUCCESS;
}
