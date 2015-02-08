!!  OpenRSP: open-ended library for response theory
!!  Copyright 2014
!!
!!  OpenRSP is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  OpenRSP is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with OpenRSP. If not, see <http://www.gnu.org/licenses/>.
!!
!!  This file implements the callback subroutine get_pert_comp_f.
!!
!!  2014-08-03, Bin Gao
!!  * first version

! configuration file of QcMatrix library
#include "qcmatrix_config.h"

#define OPENRSP_F_TEST_SRC "tests/f90/callback/get_pert_comp_f.F90"

    subroutine get_pert_comp_f(pert_label,      &
                               pert_order,      &
                               pert_rank,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                               len_ctx,         &
                               user_ctx,        &
#endif
                               pert_num_comp,   &
                               pert_components, &
                               pert_comp_orders)
        use qcmatrix_f, only: QINT
        implicit none
        integer(kind=QINT), intent(in) :: pert_label
        integer(kind=QINT), intent(in) :: pert_order
        integer(kind=QINT), intent(in) :: pert_rank
#if defined(OPENRSP_F_USER_CONTEXT)
        integer(kind=QINT), intent(in) :: len_ctx
        character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
        integer(kind=QINT), intent(out) :: pert_num_comp
        integer(kind=QINT), intent(out) :: pert_components(pert_order)
        integer(kind=QINT), intent(out) :: pert_comp_orders(pert_order)
        return
    end subroutine get_pert_comp_f

#undef OPENRSP_F_TEST_SRC
