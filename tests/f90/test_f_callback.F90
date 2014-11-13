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
!!  This file implements the callback functions for the Fortran tests.
!!
!!  2014-08-03, Bin Gao
!!  * first version

! parameters for the test suite
#include "tests/openrsp_f_test.h"
! configuration file of QMatrix library
#include "qmatrix_config.h"

#define OPENRSP_F_TEST_SRC "tests/f90/test_f_callback.F90"

subroutine get_rsp_solution(ref_ham,       &
                            ref_state,     &
                            ref_overlap,   &
                            num_freq_sums, &
                            freq_sums,     &
                            size_pert,     &
                            num_RHS,       &
                            RHS_mat,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                            len_ctx,       &
                            user_ctx,      &
#endif
                            rsp_param)
    use qmatrix, only: QINT,QREAL,QMat
    implicit none
    type(QMat), intent(in) :: ref_ham
    type(QMat), intent(in) :: ref_state
    type(QMat), intent(in) :: ref_overlap
    integer(kind=QINT), intent(in) :: num_freq_sums
    real(kind=QREAL), intent(in) :: freq_sums(num_freq_sums)
    integer(kind=QINT), intent(in) :: size_pert
    integer(kind=QINT), intent(in) :: num_RHS
    type(QMat), intent(in) :: RHS_mat(size_pert*num_RHS)
#if defined(OPENRSP_F_USER_CONTEXT)
    integer(kind=QINT), intent(in) :: len_ctx
    character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
    type(QMat), intent(inout) :: rsp_param(size_pert*num_RHS)
    return
end subroutine get_rsp_solution

#if defined(OPENRSP_PERTURBATION_FREE)
subroutine get_pert_comp(perturbation,    &
                         pert_order,      &
                         pert_rank,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                         len_ctx,         &
                         user_ctx,        &
#endif
                         pert_num_comp,   &
                         pert_components, &
                         pert_comp_orders)
    use qmatrix, only: QINT
    implicit none
    integer(kind=QINT), intent(in) :: perturbation
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
end subroutine get_pert_comp

subroutine get_pert_rank(perturbation,     &
                         pert_num_comp,    &
                         pert_components,  &
                         pert_comp_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                         len_ctx,          &
                         user_ctx,         &
#endif
                         pert_rank)
    use qmatrix, only: QINT
    implicit none
    integer(kind=QINT), intent(in) :: perturbation
    integer(kind=QINT), intent(in) :: pert_num_comp
    integer(kind=QINT), intent(in) :: pert_components(pert_num_comp)
    integer(kind=QINT), intent(in) :: pert_comp_orders(pert_num_comp)
#if defined(OPENRSP_F_USER_CONTEXT)
    integer(kind=QINT), intent(in) :: len_ctx
    character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
    integer(kind=QINT), intent(out) :: pert_rank
    return
end subroutine get_pert_rank
#endif

subroutine get_overlap_mat(bra_num_pert,      &
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
    character(len=1) :: overlap_lab(7) = (/"O","V","E","R","L","A","P"/)
#endif
    integer(kind=QINT) ipert
    integer(kind=QINT) imat
    integer(kind=QINT) ierr
#if defined(ZERO_BASED_NUMBERING)
    integer(kind=QINT) :: idx_block_row=0
    integer(kind=QINT) :: idx_block_col=0
    integer(kind=QINT) :: idx_first_row=0
    integer(kind=QINT) :: idx_first_col=0
#else
    integer(kind=QINT) :: idx_block_row=1
    integer(kind=QINT) :: idx_block_col=1
    integer(kind=QINT) :: idx_first_row=1
    integer(kind=QINT) :: idx_first_col=1
#endif
    integer(kind=QINT) :: num_row_set=2
    integer(kind=QINT) :: num_col_set=2
    real(kind=QREAL) values_real(4)
    real(kind=QREAL) values_imag(4)
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
    write(6,"(10A)") "get_overlap_mat>> label ", user_ctx
    if (all(user_ctx==overlap_lab)) then
        write(6,100) "overlap integrals"
        values_real = (/0.1,0.2,0.3,0.4/)
        values_imag = (/1.1,1.2,1.3,1.4/)
    else
        write(6,100) "unknown one-electron operator"
        stop
    end if
#endif
    write(6,100) "num_int", num_int
    do imat = 1, num_int
        ierr = QMatSetValues(val_int(imat), &
                             idx_block_row, &
                             idx_block_col, &
                             idx_first_row, &
                             num_row_set,   &
                             idx_first_col, &
                             num_col_set,   &
                             values_real,   &
                             values_imag)
        call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
    end do
    return
100 format("get_overlap_mat>> ",A,2I6)
end subroutine get_overlap_mat

subroutine get_overlap_exp(bra_num_pert,      &
                           bra_perturbations, &
                           bra_pert_orders,   &
                           ket_num_pert,      &
                           ket_perturbations, &
                           ket_pert_orders,   &
                           num_pert,          &
                           perturbations,     &
                           pert_orders,       &
                           num_dens,          &
                           ao_dens,           &
#if defined(OPENRSP_F_USER_CONTEXT)
                           len_ctx,           &
                           user_ctx,          &
#endif
                           num_exp,           &
                           val_exp)
    use qmatrix, only: QINT,QREAL,QMat
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
    integer(kind=QINT), intent(in) :: num_dens
    type(QMat), intent(in) :: ao_dens(num_dens)
#if defined(OPENRSP_F_USER_CONTEXT)
    integer(kind=QINT), intent(in) :: len_ctx
    character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
    integer(kind=QINT), intent(in) :: num_exp
    real(kind=QREAL), intent(inout) :: val_exp(num_exp)
#if defined(OPENRSP_F_USER_CONTEXT)
    character(len=1) :: overlap_lab(7) = (/"O","V","E","R","L","A","P"/)
#endif
    integer(kind=QINT) ipert
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
    write(6,100) "num_dens", num_dens
#if defined(OPENRSP_F_USER_CONTEXT)
    write(6,"(10A)") "get_overlap_exp>> label ", user_ctx
    if (all(user_ctx==overlap_lab)) then
        write(6,100) "overlap integrals"
    else
        write(6,100) "unknown one-electron operator"
        stop
    end if
#endif
    return
100 format("get_overlap_exp>> ",A,2I6)
end subroutine get_overlap_exp

subroutine get_one_oper_mat(num_pert,      &
                            perturbations, &
                            pert_orders,   &
#if defined(OPENRSP_F_USER_CONTEXT)
                            len_ctx,       &
                            user_ctx,      &
#endif
                            num_int,       &
                            val_int)
    use qmatrix, only: QINT,  &
                       QREAL, &
                       QMat,  &
                       QMatSetValues
    implicit none
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
    character(len=1) :: oneham_lab(6) = (/"O","N","E","H","A","M"/)
    character(len=1) :: ext_field_lab(9) = (/"E","X","T","_","F","I", "E", "L", "D"/)
#endif
    integer(kind=QINT) ipert
    integer(kind=QINT) imat
    integer(kind=QINT) ierr
#if defined(ZERO_BASED_NUMBERING)
    integer(kind=QINT) :: idx_block_row=0
    integer(kind=QINT) :: idx_block_col=0
    integer(kind=QINT) :: idx_first_row=0
    integer(kind=QINT) :: idx_first_col=0
#else
    integer(kind=QINT) :: idx_block_row=1
    integer(kind=QINT) :: idx_block_col=1
    integer(kind=QINT) :: idx_first_row=1
    integer(kind=QINT) :: idx_first_col=1
#endif
    integer(kind=QINT) :: num_row_set=2
    integer(kind=QINT) :: num_col_set=2
    real(kind=QREAL) values_real(4)
    real(kind=QREAL) values_imag(4)
    write(6,100) "num_pert", num_pert
    do ipert = 1, num_pert
        write(6,100) "pert.",              &
                     perturbations(ipert), &
                     pert_orders(ipert)
    end do
#if defined(OPENRSP_F_USER_CONTEXT)
    write(6,"(10A)") "get_one_oper_mat>> label ", user_ctx
    if (all(user_ctx==oneham_lab)) then
        write(6,100) "one-electron Hamiltonian"
        values_real = (/0.1,0.2,0.3,0.4/)
        values_imag = (/1.1,1.2,1.3,1.4/)
    else if (all(user_ctx==ext_field_lab)) then
        write(6,100) "external field"
        values_real = (/2.1,2.2,2.3,2.4/)
        values_imag = (/3.1,3.2,3.3,3.4/)
    else
        write(6,100) "unknown one-electron operator"
        stop
    end if
#endif
    write(6,100) "num_int", num_int
    do imat = 1, num_int
        ierr = QMatSetValues(val_int(imat), &
                             idx_block_row, &
                             idx_block_col, &
                             idx_first_row, &
                             num_row_set,   &
                             idx_first_col, &
                             num_col_set,   &
                             values_real,   &
                             values_imag)
        call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
    end do
    return
100 format("get_one_oper_mat>> ",A,2I6)
end subroutine get_one_oper_mat

subroutine get_one_oper_exp(num_pert,      &
                            perturbations, &
                            pert_orders,   &
                            num_dens,      &
                            ao_dens,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                            len_ctx,       &
                            user_ctx,      &
#endif
                            num_exp,       &
                            val_exp)
    use qmatrix, only: QINT,QREAL,QMat
    implicit none
    integer(kind=QINT), intent(in) :: num_pert
    integer(kind=QINT), intent(in) :: perturbations(num_pert)
    integer(kind=QINT), intent(in) :: pert_orders(num_pert)
    integer(kind=QINT), intent(in) :: num_dens
    type(QMat), intent(in) :: ao_dens(num_dens)
#if defined(OPENRSP_F_USER_CONTEXT)
    integer(kind=QINT), intent(in) :: len_ctx
    character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
    integer(kind=QINT), intent(in) :: num_exp
    real(kind=QREAL), intent(inout) :: val_exp(num_exp)
#if defined(OPENRSP_F_USER_CONTEXT)
    character(len=1) :: oneham_lab(6) = (/"O","N","E","H","A","M"/)
    character(len=1) :: ext_field_lab(9) = (/"E","X","T","_","F","I", "E", "L", "D"/)
#endif
    integer(kind=QINT) ipert
    write(6,100) "num_pert", num_pert
    do ipert = 1, num_pert
        write(6,100) "pert.",              &
                     perturbations(ipert), &
                     pert_orders(ipert)
    end do
    write(6,100) "num_dens", num_dens
#if defined(OPENRSP_F_USER_CONTEXT)
    write(6,"(10A)") "get_one_oper_exp>> label ", user_ctx
    if (all(user_ctx==oneham_lab)) then
        write(6,100) "one-electron Hamiltonian"
    else if (all(user_ctx==ext_field_lab)) then
        write(6,100) "external field"
    else
        write(6,100) "unknown one-electron operator"
        stop
    end if
#endif
    return
100 format("get_one_oper_exp>> ",A,2I6)
end subroutine get_one_oper_exp

#undef OPENRSP_F_TEST_SRC
