!!  OpenRSP: open-ended library for response theory
!!  Copyright 2014
!!
!!  This source code form is subject to the terms of the
!!  GNU Lesser General Public License, version 2.1.
!!  If a copy of the GNU LGPL v2.1 was not distributed with this
!!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!!
!!  This file implements the callback subroutine get_overlap_exp_f.
!!
!!  2014-08-03, Bin Gao
!!  * first version

! configuration file of QcMatrix library
#include "qcmatrix_config.h"

#define OPENRSP_F_TEST_SRC "tests/f90/callback/get_overlap_exp_f.F90"

    subroutine get_overlap_exp_f(bra_num_pert,    &
                                 bra_pert_labels, &
                                 bra_pert_orders, &
                                 ket_num_pert,    &
                                 ket_pert_labels, &
                                 ket_pert_orders, &
                                 num_pert,        &
                                 pert_labels,     &
                                 pert_orders,     &
                                 num_dens,        &
                                 ao_dens,         &
#if defined(OPENRSP_F_USER_CONTEXT)
                                 len_ctx,         &
                                 user_ctx,        &
#endif
                                 num_exp,         &
                                 val_exp)
        use qcmatrix_f, only: QINT,                   &
                              QREAL,                  &
                              QSYMMAT,                &
                              QREALMAT,               &
                              MAT_NO_OPERATION,       &
                              QcMat,                  &
                              QcMatCreate_f,          &
                              QcMatBlockCreate_f,     &
                              QcMatSetSymType_f,      &
                              QcMatSetDataType_f,     &
                              QcMatSetDimMat_f,       &
                              QcMatAssemble_f,        &
                              QcMatSetValues_f,       &
                              QcMatGetMatProdTrace_f, &
                              QcMatDestroy_f
        implicit none
        integer(kind=QINT), intent(in) :: bra_num_pert
        integer(kind=QINT), intent(in) :: bra_pert_labels(bra_num_pert)
        integer(kind=QINT), intent(in) :: bra_pert_orders(bra_num_pert)
        integer(kind=QINT), intent(in) :: ket_num_pert
        integer(kind=QINT), intent(in) :: ket_pert_labels(ket_num_pert)
        integer(kind=QINT), intent(in) :: ket_pert_orders(ket_num_pert)
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_orders(num_pert)
        integer(kind=QINT), intent(in) :: num_dens
        type(QcMat), intent(in) :: ao_dens(num_dens)
#if defined(OPENRSP_F_USER_CONTEXT)
        integer(kind=QINT), intent(in) :: len_ctx
        character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
        integer(kind=QINT), intent(in) :: num_exp
        real(kind=QREAL), intent(inout) :: val_exp(2*num_exp)
! defined perturbations and their maximum orders
#include "tests/openrsp_f_perturbations.h90"
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: overlap_context(7) = (/"O","V","E","R","L","A","P"/)
#endif
#include "tests/ao_dens/openrsp_f_ao_dims.h90"
#include "tests/ao_dens/openrsp_f_ao_overlap.h90"
        type(QcMat) val_int(1)
        integer(kind=QINT) idens
        integer(kind=4) ierr
#if defined(OPENRSP_F_USER_CONTEXT)
        if (all(user_ctx/=overlap_context)) then
            write(6,100) "unknown one-electron operator"
            call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
        end if
#endif
        ! overlap integrals
        if (num_pert==0 .and. bra_num_pert==0 .and. ket_num_pert==0) then
            ierr = QcMatCreate_f(A=val_int(1))
            call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
            ierr = QcMatBlockCreate_f(A=val_int(1), dim_block=1_QINT)
            call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
            ierr = QcMatSetSymType_f(A=val_int(1), sym_type=QSYMMAT)
            call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
            ierr = QcMatSetDataType_f(A=val_int(1),                    &
                                      num_blocks=1_QINT,               &
                                      idx_block_row=(/IDX_BLOCK_ROW/), &
                                      idx_block_col=(/IDX_BLOCK_COL/), &
                                      data_type=(/QREALMAT/))
            call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
            ierr = QcMatSetDimMat_f(A=val_int(1),   &
                                    num_row=NUM_AO, &
                                    num_col=NUM_AO)
            call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
            ierr = QcMatAssemble_f(A=val_int(1))
            call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
            ierr = QcMatSetValues_f(A=val_int(1),                &
                                    idx_block_row=IDX_BLOCK_ROW, &
                                    idx_block_col=IDX_BLOCK_COL, &
                                    idx_first_row=IDX_FIRST_ROW, &
                                    num_row_set=NUM_AO,          &
                                    idx_first_col=IDX_FIRST_COL, &
                                    num_col_set=NUM_AO,          &
                                    values_real=values_overlap)
            call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
            do idens = 1, num_dens
                ierr = QcMatGetMatProdTrace_f(A=val_int(1),          &
                                              B=ao_dens(idens),      &
                                              op_B=MAT_NO_OPERATION, &
                                              num_blocks=1_QINT,     &
                                              trace=val_exp(2*idens-1:2*idens))
                call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
            end do
            ierr = QcMatDestroy_f(A=val_int(1))
            call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        else if (num_pert==1 .and. bra_num_pert==0 .and. ket_num_pert==0) then
            if (pert_labels(1)==PERT_GEOMETRIC) then
                write(6,100) "not implemented"
                call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
            ! zero integrals
            else if (pert_labels(1)==PERT_DIPOLE) then
                val_exp = 0
            else if (pert_labels(1)==PERT_MAGNETIC) then
                write(6,100) "not implemented"
                call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
            else
                write(6,100) "unknown perturbations"
                call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
            end if
        else if (num_pert==0 .and. bra_num_pert==1 .and. ket_num_pert==0) then
            if (bra_pert_labels(1)==PERT_GEOMETRIC) then
                write(6,100) "not implemented"
                call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
            ! zero integrals
            else if (bra_pert_labels(1)==PERT_DIPOLE) then
                val_exp = 0
            else if (bra_pert_labels(1)==PERT_MAGNETIC) then
                write(6,100) "not implemented"
                call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
            else
                write(6,100) "unknown perturbations"
                call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
            end if
        else if (num_pert==0 .and. bra_num_pert==0 .and. ket_num_pert==1) then
            if (ket_pert_labels(1)==PERT_GEOMETRIC) then
                write(6,100) "not implemented"
                call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
            ! zero integrals
            else if (ket_pert_labels(1)==PERT_DIPOLE) then
                val_exp = 0
            else if (ket_pert_labels(1)==PERT_MAGNETIC) then
                write(6,100) "not implemented"
                call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
            else
                write(6,100) "unknown perturbations"
                call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
            end if
        else if (num_pert==0 .and. bra_num_pert==1 .and. ket_num_pert==1) then
            ! zero integrals
            if (bra_pert_labels(1)==PERT_DIPOLE .or. ket_pert_labels(1)==PERT_DIPOLE) then
                val_exp = 0
            else
                write(6,100) "unknown perturbations"
                call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
            end if
        else
            write(6,100) "not implemented"
            call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
        end if
        return
100     format("get_overlap_exp_f>> ",A)
    end subroutine get_overlap_exp_f

#undef OPENRSP_F_TEST_SRC
