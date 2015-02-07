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
!!  This file implements the callback subroutine get_two_oper_mat_f.
!!
!!  2014-08-03, Bin Gao
!!  * first version

! configuration file of QMatrix library
#include "qmatrix_config.h"

#define OPENRSP_F_TEST_SRC "tests/f90/callback/get_two_oper_mat_f.F90"

    subroutine get_two_oper_mat_f(num_pert,    &
                                  pert_labels, &
                                  pert_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                  len_ctx,     &
                                  user_ctx,    &
#endif
                                  num_int,     &
                                  val_int)
        use qmatrix, only: QINT,            &
                           QREAL,           &
                           QSYMMAT,         &
                           QREALMAT,        &
                           QMat,            &
                           QMatIsAssembled, &
                           QMatBlockCreate, &
                           QMatSetSymType,  &
                           QMatSetDataType, &
                           QMatSetDimMat,   &
                           QMatAssemble,    &
                           QMatSetValues,   &
                           QMatZeroEntries
        implicit none
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_orders(num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
        integer(kind=QINT), intent(in) :: len_ctx
        character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
        integer(kind=QINT), intent(in) :: num_int
        type(QMat), intent(inout) :: val_int(num_int)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: oneham_context(6) = (/"O","N","E","H","A","M"/)
        character(len=1) :: ext_field_context(9) = (/"E","X","T","_","F","I", "E", "L", "D"/)
#endif
#include "tests/ao_dens/openrsp_f_ao_dims.h90"
#include "tests/ao_dens/openrsp_f_ao_diplen.h90"
        logical(kind=4) assembled
        integer(kind=QINT) imat
        integer(kind=4) ierr
#if defined(OPENRSP_F_USER_CONTEXT)
        if (all(user_ctx==oneham_context)) then
            write(6,100) "not implemented"
            call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
        else if (all(user_ctx==ext_field_context)) then
            ! electric fields
            if (num_pert==1 .and. pert_labels(1)==PERT_DIPOLE) then
                ! checks if the matrix is assembled or not
                do imat = 1, num_int
                    ierr = QMatIsAssembled(A=val_int(imat), assembled=assembled)
                    call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                    if (.not.assembled) then
                        ierr = QMatBlockCreate(A=val_int(imat), dim_block=1_QINT)
                        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                        ierr = QMatSetSymType(A=val_int(imat), sym_type=QSYMMAT)
                        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                        ierr = QMatSetDataType(A=val_int(imat),                 &
                                               num_blocks=1_QINT,               &
                                               idx_block_row=(/IDX_BLOCK_ROW/), &
                                               idx_block_col=(/IDX_BLOCK_COL/), &
                                               data_type=(/QREALMAT/))
                        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                        ierr = QMatSetDimMat(A=val_int(imat), dim_mat=NUM_AO)
                        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                        ierr = QMatAssemble(A=val_int(imat))
                        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                    end if
                end do
                ! dipole length integrals
                if (pert_orders(1)==1) then
                    do imat = 1, 3
                        ierr = QMatSetValues(A=val_int(imat),             &
                                             idx_block_row=IDX_BLOCK_ROW, &
                                             idx_block_col=IDX_BLOCK_COL, &
                                             idx_first_row=IDX_FIRST_ROW, &
                                             num_row_set=NUM_AO,          &
                                             idx_first_col=IDX_FIRST_COL, &
                                             num_col_set=NUM_AO,          &
                                             values_real=values_diplen((imat-1)*SIZE_AO_MAT+1:imat*SIZE_AO_MAT)
                        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                    end do
                ! zero integrals
                else
                    do imat = 1, num_int
                        ierr = QMatZeroEntries(A=val_int(imat))
                        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                    end do
                end if
            else
                write(6,100) "not implemented"
                call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
            end if
        else
            write(6,100) "unknown one-electron operator"
            call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
        end if
#else
        ! electric fields
        if (num_pert==1 .and. pert_labels(1)==PERT_DIPOLE) then
            ! checks if the matrix is assembled or not
            do imat = 1, num_int
                ierr = QMatIsAssembled(A=val_int(imat), assembled=assembled)
                call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                if (.not.assembled) then
                    ierr = QMatBlockCreate(A=val_int(imat), dim_block=1_QINT)
                    call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                    ierr = QMatSetSymType(A=val_int(imat), sym_type=QSYMMAT)
                    call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                    ierr = QMatSetDataType(A=val_int(imat),                 &
                                           num_blocks=1_QINT,               &
                                           idx_block_row=(/IDX_BLOCK_ROW/), &
                                           idx_block_col=(/IDX_BLOCK_COL/), &
                                           data_type=(/QREALMAT/))
                    call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                    ierr = QMatSetDimMat(A=val_int(imat), dim_mat=NUM_AO)
                    call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                    ierr = QMatAssemble(A=val_int(imat))
                    call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                end if
            end do
            ! dipole length integrals
            if (pert_orders(1)==1) then
                do imat = 1, 3
                    ierr = QMatSetValues(A=val_int(imat),             &
                                         idx_block_row=IDX_BLOCK_ROW, &
                                         idx_block_col=IDX_BLOCK_COL, &
                                         idx_first_row=IDX_FIRST_ROW, &
                                         num_row_set=NUM_AO,          &
                                         idx_first_col=IDX_FIRST_COL, &
                                         num_col_set=NUM_AO,          &
                                         values_real=values_diplen((imat-1)*SIZE_AO_MAT+1:imat*SIZE_AO_MAT)
                    call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                end do
            ! zero integrals
            else
                do imat = 1, num_int
                    ierr = QMatZeroEntries(A=val_int(imat))
                    call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
                end do
            end if
        else
            write(6,100) "not implemented"
            call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
        end if
#endif
        return
100     format("get_two_oper_mat_f>> ",A)
    end subroutine get_two_oper_mat_f

#undef OPENRSP_F_TEST_SRC
