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
!!  This file implements the callback subroutine get_overlap_mat_f.
!!
!!  2014-08-03, Bin Gao
!!  * first version

! configuration file of QMatrix library
#include "qmatrix_config.h"

#define OPENRSP_F_TEST_SRC "tests/f90/callback/get_overlap_mat_f.F90"

    subroutine get_overlap_mat_f(bra_num_pert,      &
                                 bra_perturbations, &
                                 bra_pert_orders,   &
                                 ket_num_pert,      &
                                 ket_perturbations, &
                                 ket_pert_orders,   &
                                 num_pert,          &
                                 perturbations,     &
                                 pert_orders,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                                 len_ctx,           &
                                 user_ctx,          &
#endif
                                 num_int,           &
                                 val_int)
        use qmatrix, only: QINT,  &
                           QREAL, &
                           QMat,  &
                           QMatSetValues
        implicit none
        integer(kind=QINT), intent(in) :: bra_num_pert
        integer(kind=QINT), intent(in) :: bra_perturbations(bra_num_pert)
        integer(kind=QINT), intent(in) :: bra_pert_orders(bra_num_pert)
        integer(kind=QINT), intent(in) :: ket_num_pert
        integer(kind=QINT), intent(in) :: ket_perturbations(ket_num_pert)
        integer(kind=QINT), intent(in) :: ket_pert_orders(ket_num_pert)
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: perturbations(num_pert)
        integer(kind=QINT), intent(in) :: pert_orders(num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
        integer(kind=QINT), intent(in) :: len_ctx
        character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
        integer(kind=QINT), intent(in) :: num_int
        type(QMat), intent(inout) :: val_int(num_int)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: overlap_context(7) = (/"O","V","E","R","L","A","P"/)
#endif
#include "tests/openrsp_f_AO_dim.h90"
        integer(kind=QINT) ipert
        integer(kind=QINT) imat
        integer(kind=4) ierr
        write(6,100) "bra_num_pert", bra_num_pert
        do ipert = 1, bra_num_pert
            write(6,100) "bra_pert.",              &
                         bra_perturbations(ipert), &
                         bra_pert_orders(ipert)
        end do
        write(6,100) "ket_num_pert", ket_num_pert
        do ipert = 1, ket_num_pert
            write(6,100) "ket_pert.",              &
                         ket_perturbations(ipert), &
                         ket_pert_orders(ipert)
        end do
        write(6,100) "num_pert", num_pert
        do ipert = 1, num_pert
            write(6,100) "pert.",              &
                         perturbations(ipert), &
                         pert_orders(ipert)
        end do
#if defined(OPENRSP_F_USER_CONTEXT)
        write(6,"(10A)") "get_one_oper_mat_f>> label ", user_ctx
        if (all(user_ctx==overlap_context)) then
            write(6,100) "overlap integrals"
        else
            write(6,100) "unknown one-electron operator"
            stop
        end if
#endif
        return
100     format("get_overlap_mat_f>> ",A,2I6)
    end subroutine get_overlap_mat_f

#undef OPENRSP_F_TEST_SRC
