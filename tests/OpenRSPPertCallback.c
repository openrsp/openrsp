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

#include "OpenRSPPertCallback.h"

void get_pert_concatenation(const QInt pert_label,
                            const QcPertInt first_cat_comp,
                            const QInt num_cat_comps,
                            const QInt num_sub_tuples,
                            const QInt *len_sub_tuples,
#if defined(OPENRSP_C_USER_CONTEXT)
                            void *user_ctx,
#endif
                            QInt *rank_sub_comps)
{
}
