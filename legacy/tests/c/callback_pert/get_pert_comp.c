/* OpenRSP: open-ended library for response theory
   Copyright 2014
  This source code form is subject to the terms of the
  GNU Lesser General Public License, version 2.1.
  If a copy of the GNU LGPL v2.1 was not distributed with this
  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

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
