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

   This file implements the callback function get_pert_comp.

   2014-07-31, Bin Gao:
   * first version
*/

#include "tests/openrsp_c_pert_callback.h"

QVoid get_pert_comp(const QInt pert_label,
                    const QInt pert_order,
                    const QInt pert_rank,
#if defined(OPENRSP_C_USER_CONTEXT)
                    QVoid *user_ctx,
#endif
                    QInt *pert_num_comp,
                    QInt *pert_components,
                    QInt *pert_comp_orders)
{
}
