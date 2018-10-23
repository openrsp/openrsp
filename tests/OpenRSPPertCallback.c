/*
  OpenRSP: open-ended library for response theory
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
/*FIXME: to implement*/
}
