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
!!  This file implements the callback subroutine get_rsp_solution_f.
!!
!!  2014-08-03, Bin Gao
!!  * first version

! configuration file of QcMatrix library
#include "qcmatrix_config.h"

#define OPENRSP_F_TEST_SRC "tests/f90/callback/get_rsp_solution_f.F90"

    subroutine get_rsp_solution_f(ref_ham,       &
                                  ref_state,     &
                                  ref_overlap,   &
                                  num_freq_sums, &
                                  freq_sums,     &
                                  size_pert,     &
                                  RHS_mat,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                                  len_ctx,       &
                                  user_ctx,      &
#endif
                                  rsp_param)
        use qcmatrix_f, only: QINT,               &
                              QREAL,              &
                              QREALMAT,           &
                              QcMat,              &
                              QcMatIsAssembled_f, &
                              QcMatBlockCreate_f, &
                              QcMatSetDataType_f, &
                              QcMatSetDimMat_f,   &
                              QcMatAssemble_f,    &
                              QcMatSetValues_f
        implicit none
        type(QcMat), intent(in) :: ref_ham
        type(QcMat), intent(in) :: ref_state
        type(QcMat), intent(in) :: ref_overlap
        integer(kind=QINT), intent(in) :: num_freq_sums
        real(kind=QREAL), intent(in) :: freq_sums(num_freq_sums)
        integer(kind=QINT), intent(in) :: size_pert
        type(QcMat), intent(in) :: RHS_mat(size_pert*num_freq_sums)
#if defined(OPENRSP_F_USER_CONTEXT)
        integer(kind=QINT), intent(in) :: len_ctx
        character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
        type(QcMat), intent(inout) :: rsp_param(size_pert*num_freq_sums)
#include "tests/ao_dens/openrsp_f_ao_dims.h90"
#include "tests/ao_dens/openrsp_f_ao_alpha_rsp_param.h90"
        logical(kind=4) assembled
        integer(kind=4) ierr
        integer(kind=QINT), save :: id_rsp_param = 0
        id_rsp_param = id_rsp_param+1
        if (id_rsp_param>3) then
            write(6,100) "not implemented"
            call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
        end if
        ! checks if the matrix is assembled or not
        ierr = QcMatIsAssembled_f(A=rsp_param(1), assembled=assembled)
        call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        if (.not.assembled) then
            ierr = QcMatBlockCreate_f(A=rsp_param(1), dim_block=1_QINT)
            call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
            ierr = QcMatSetDataType_f(A=rsp_param(1),                  &
                                      num_blocks=1_QINT,               &
                                      idx_block_row=(/IDX_BLOCK_ROW/), &
                                      idx_block_col=(/IDX_BLOCK_COL/), &
                                      data_type=(/QREALMAT/))
            call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
            ierr = QcMatSetDimMat_f(A=rsp_param(1), &
                                    num_row=NUM_AO, &
                                    num_col=NUM_AO)
            call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
            ierr = QcMatAssemble_f(A=rsp_param(1))
            call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        end if
        ierr = QcMatSetValues_f(A=rsp_param(1),              &
                                idx_block_row=IDX_BLOCK_ROW, &
                                idx_block_col=IDX_BLOCK_COL, &
                                idx_first_row=IDX_FIRST_ROW, &
                                num_row_set=NUM_AO,          &
                                idx_first_col=IDX_FIRST_COL, &
                                num_col_set=NUM_AO,          &
                                values_real=alpha_rsp_param((id_rsp_param-1)*SIZE_AO_MAT+1:id_rsp_param*SIZE_AO_MAT))
        call QErrorCheckCode(6, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        return
100     format("get_rsp_solution_f>> ",A)
    end subroutine get_rsp_solution_f

#undef OPENRSP_F_TEST_SRC
