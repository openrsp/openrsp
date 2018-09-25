!!  OpenRSP: open-ended library for response theory
!!  Copyright 2014
!!
!!  This source code form is subject to the terms of the
!!  GNU Lesser General Public License, version 2.1.
!!  If a copy of the GNU LGPL v2.1 was not distributed with this
!!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!!
!!  This file implements the callback subroutine get_pert_rank_f.
!!
!!  2014-08-03, Bin Gao
!!  * first version

! configuration file of QcMatrix library
#include "qcmatrix_config.h"

#define OPENRSP_F_TEST_SRC "tests/f90/callback/get_pert_rank_f.F90"

    subroutine get_pert_rank_f(pert_label,       &
                               pert_num_comp,    &
                               pert_components,  &
                               pert_comp_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                               len_ctx,          &
                               user_ctx,         &
#endif
                               pert_rank)
        use qcmatrix_f, only: QINT
        implicit none
        integer(kind=QINT), intent(in) :: pert_label
        integer(kind=QINT), intent(in) :: pert_num_comp
        integer(kind=QINT), intent(in) :: pert_components(pert_num_comp)
        integer(kind=QINT), intent(in) :: pert_comp_orders(pert_num_comp)
#if defined(OPENRSP_F_USER_CONTEXT)
        integer(kind=QINT), intent(in) :: len_ctx
        character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
        integer(kind=QINT), intent(out) :: pert_rank
        return
    end subroutine get_pert_rank_f

#undef OPENRSP_F_TEST_SRC
