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
!!  This file tests the density matrix-based response theory.
!!
!!  2014-08-03, Bin Gao
!!  * first version

! parameters for the test suite
#include "tests/openrsp_f_test.h"

#define OPENRSP_F_TEST_SRC "tests/f90/test_f_OpenRSP_AO.F90"

    subroutine test_f_OpenRSP_AO(open_rsp, io_log)
        use qmatrix, only: QINT,            &
                           QREAL,           &
                           QSYMMAT,         &
                           QREALMAT,        &
                           QMat,            &
                           QMatCreate,      &
                           QMatBlockCreate, &
                           QMatSetSymType,  &
                           QMatSetDataType, &
                           QMatSetDimMat,   &
                           QMatAssemble,    &
                           QMatSetValues,   &
                           QMatDestroy
        use openrsp_f, only: ELEC_AO_D_MATRIX,    &
                             OpenRSP,             &
                             OpenRSPSetElecEOM_f, &
                             OpenRSPSetSolver_f,  &
                             OpenRSPSetPDBS_f,    &
                             OpenRSPAddOneOper_f, &
                             OpenRSPAssemble_f,   &
                             OpenRSPWrite_f,      &
                             OpenRSPGetRSPFun_f
        implicit none
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=4), intent(in) :: io_log
! defined perturbations and their maximum orders
#include "tests/openrsp_f_perturbations.h90"
#if defined(OPENRSP_F_USER_CONTEXT)
        ! user defined context of linear response equation solver
        character(len=1) :: solver_context(6) = (/"S","O","L","V","E","R"/)
#endif
        ! callback subroutine of linear response equation solver
        external get_rsp_solution_f
        ! overlap integrals with London atomic orbitals
        integer(kind=QINT), parameter :: overlap_num_pert = 2
        integer(kind=QINT) :: overlap_perturbations(overlap_num_pert) = (/ &
            PERT_GEOMETRIC,PERT_MAGNETIC/)
        integer(kind=QINT) :: overlap_pert_orders(overlap_num_pert) = (/ &
            MAX_ORDER_GEOMETRIC,MAX_ORDER_MAGNETIC/)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: overlap_context(7) = (/"O","V","E","R","L","A","P"/)
#endif
        external get_overlap_mat_f
        external get_overlap_exp_f
        ! one-electron Hamiltonian
        integer(kind=QINT), parameter :: oneham_num_pert = 2
        integer(kind=QINT) :: oneham_perturbations(oneham_num_pert) = (/ &
            PERT_GEOMETRIC,PERT_MAGNETIC/)
        integer(kind=QINT) :: oneham_pert_orders(oneham_num_pert) = (/ &
            MAX_ORDER_GEOMETRIC,MAX_ORDER_MAGNETIC/)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: oneham_context(6) = (/"O","N","E","H","A","M"/)
#endif
        ! external field
        integer(kind=QINT), parameter :: ext_field_num_pert = 3
        integer(kind=QINT) :: ext_field_perturbations(ext_field_num_pert) = (/ &
            PERT_GEOMETRIC,PERT_DIPOLE,PERT_MAGNETIC/)
        integer(kind=QINT) :: ext_field_pert_orders(ext_field_num_pert) = (/ &
            MAX_ORDER_GEOMETRIC,MAX_ORDER_DIPOLE,MAX_ORDER_MAGNETIC/)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: ext_field_context(9) = (/"E","X","T","_","F","I", "E", "L", "D"/)
#endif
        external get_one_oper_mat_f
        external get_one_oper_exp_f
        ! referenced state
#include "tests/openrsp_f_AO_dim.h90"
#include "tests/openrsp_f_AO_state.h90"
        type(QMat) F_unpert
        type(QMat) D_unpert
        type(QMat) S_unpert
        ! polarizability
        integer(kind=QINT), parameter :: ALPHA_NUM_PERT = 1_QINT
        integer(kind=QINT), parameter :: ALPHA_PERTUBRATIONS(ALPHA_NUM_PERT) = (/PERT_DIPOLE/)
        integer(kind=QINT), parameter :: ALPHA_PERT_ORDER(ALPHA_NUM_PERT) = (/1_QINT/)
        real(kind=QREAL), parameter :: ALPHA_PERT_FREQ(2*ALPHA_NUM_PERT) = (/0.1_QREAL,0.0_QREAL/)
        ! kn rule and response functions
        integer(kind=QINT) kn_rule(2)
        integer(kind=QINT) size_rsp_fun
        real(kind=QREAL) rsp_fun(18)
        ! error information
        integer(kind=4) ierr

        ! sets the equation of motion of electrons
        ierr = OpenRSPSetElecEOM_f(open_rsp, ELEC_AO_D_MATRIX)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPSetElecEOM_f() passed"

        ! sets the context of linear response equation solver
        ierr = OpenRSPSetSolver_f(open_rsp,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                                  solver_context, &
#endif
                                  get_rsp_solution_f)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPSetSolver_f() passed"

        ! sets the context of perturbation dependent basis sets
        ierr = OpenRSPSetPDBS_f(open_rsp,              &
                                overlap_num_pert,      &
                                overlap_perturbations, &
                                overlap_pert_orders,   &
#if defined(OPENRSP_F_USER_CONTEXT)
                                overlap_context,       &
#endif
                                get_overlap_mat_f,     &
                                get_overlap_exp_f)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPSetPDBS_f() passed"

        ! adds one-electron Hamiltonian
        ierr = OpenRSPAddOneOper_f(open_rsp,             &
                                   oneham_num_pert,      &
                                   oneham_perturbations, &
                                   oneham_pert_orders,   &
#if defined(OPENRSP_F_USER_CONTEXT)
                                   oneham_context,       &
#endif
                                   get_one_oper_mat_f,   &
                                   get_one_oper_exp_f)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPAddOneOper_f(h) passed"

        ! adds external field
        ierr = OpenRSPAddOneOper_f(open_rsp,                &
                                   ext_field_num_pert,      &
                                   ext_field_perturbations, &
                                   ext_field_pert_orders,   &
#if defined(OPENRSP_F_USER_CONTEXT)
                                   ext_field_context,       &
#endif
                                   get_one_oper_mat_f,      &
                                   get_one_oper_exp_f)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPAddOneOper_f(V) passed"

        ! assembles the context of response theory calculations
        ierr = OpenRSPAssemble_f(open_rsp)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPAssemble_f() passed"

        ! writes the context of response theory calculations
        ierr = OpenRSPWrite_f(open_rsp, OPENRSP_F_LOG)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPWrite_f() passed"

        ! sets the unperturbed Fock matrix
        ierr = QMatCreate(A=F_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatBlockCreate(A=F_unpert, dim_block=1_QINT)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatSetSymType(A=F_unpert, sym_type=QSYMMAT)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatSetDataType(A=F_unpert,                      &
                               num_blocks=1_QINT,               &
                               idx_block_row=(/IDX_BLOCK_ROW/), &
                               idx_block_col=(/IDX_BLOCK_COL/), &
                               data_type=(/QREALMAT/))
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatSetDimMat(A=F_unpert, dim_mat=NUM_ROW_SET)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatAssemble(A=F_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatSetValues(A=F_unpert,                  &
                             idx_block_row=IDX_BLOCK_ROW, &
                             idx_block_col=IDX_BLOCK_COL, &
                             idx_first_row=IDX_FIRST_ROW, &
                             num_row_set=NUM_ROW_SET,     &
                             idx_first_col=IDX_FIRST_COL, &
                             num_col_set=NUM_COL_SET,     &
                             values_real=values_Fock)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ! sets the unperturbed density matrix
        ierr = QMatCreate(A=D_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatBlockCreate(A=D_unpert, dim_block=1_QINT)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatSetSymType(A=D_unpert, sym_type=QSYMMAT)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatSetDataType(A=D_unpert,                      &
                               num_blocks=1_QINT,               &
                               idx_block_row=(/IDX_BLOCK_ROW/), &
                               idx_block_col=(/IDX_BLOCK_COL/), &
                               data_type=(/QREALMAT/))
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatSetDimMat(A=D_unpert, dim_mat=NUM_ROW_SET)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatAssemble(A=D_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatSetValues(A=D_unpert,                  &
                             idx_block_row=IDX_BLOCK_ROW, &
                             idx_block_col=IDX_BLOCK_COL, &
                             idx_first_row=IDX_FIRST_ROW, &
                             num_row_set=NUM_ROW_SET,     &
                             idx_first_col=IDX_FIRST_COL, &
                             num_col_set=NUM_COL_SET,     &
                             values_real=values_density)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ! sets the unperturbed overlap integrals
        ierr = QMatCreate(A=S_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatBlockCreate(A=S_unpert, dim_block=1_QINT)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatSetSymType(A=S_unpert, sym_type=QSYMMAT)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatSetDataType(A=S_unpert,                      &
                               num_blocks=1_QINT,               &
                               idx_block_row=(/IDX_BLOCK_ROW/), &
                               idx_block_col=(/IDX_BLOCK_COL/), &
                               data_type=(/QREALMAT/))
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatSetDimMat(A=S_unpert, dim_mat=NUM_ROW_SET)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatAssemble(A=S_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatSetValues(A=S_unpert,                  &
                             idx_block_row=IDX_BLOCK_ROW, &
                             idx_block_col=IDX_BLOCK_COL, &
                             idx_first_row=IDX_FIRST_ROW, &
                             num_row_set=NUM_ROW_SET,     &
                             idx_first_col=IDX_FIRST_COL, &
                             num_col_set=NUM_COL_SET,     &
                             values_real=values_overlap)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)

        ! gets the polarizability
        kn_rule = (/0,1/)
        size_rsp_fun = 9
        ierr = OpenRSPGetRSPFun_f(open_rsp,            &
                                  F_unpert,            &
                                  D_unpert,            &
                                  S_unpert,            &
                                  ALPHA_NUM_PERT,      &
                                  ALPHA_PERTUBRATIONS, &
                                  ALPHA_PERT_ORDER,    &
                                  ALPHA_PERT_FREQ,     &
                                  kn_rule,             &
                                  size_rsp_fun,        &
                                  rsp_fun)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPGetRSPFun_f() passed"

        ! cleans
        ierr = QMatDestroy(A=F_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatDestroy(A=D_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QMatDestroy(A=S_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)

        return

100     format("test_f_OpenRSP_AO>> ",A)
    end subroutine test_f_OpenRSP_AO

#undef OPENRSP_F_TEST_SRC
