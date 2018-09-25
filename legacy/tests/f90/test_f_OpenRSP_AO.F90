!!  OpenRSP: open-ended library for response theory
!!  Copyright 2014
!!
!!  This source code form is subject to the terms of the
!!  GNU Lesser General Public License, version 2.1.
!!  If a copy of the GNU LGPL v2.1 was not distributed with this
!!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!!
!!  This file tests the atomic orbital density matrix-based response theory.
!!
!!  2014-08-03, Bin Gao
!!  * first version

! parameters for the test suite
#include "tests/openrsp_test_param.h"

#define OPENRSP_F_TEST_SRC "tests/f90/test_f_OpenRSP_AO.F90"

    subroutine test_f_OpenRSP_AO(open_rsp, io_log)
        use qcmatrix_f, only: QINT,               &
                              QREAL,              &
                              QSYMMAT,            &
                              QREALMAT,           &
                              QcMat,              &
                              QcMatCreate_f,      &
                              QcMatBlockCreate_f, &
                              QcMatSetSymType_f,  &
                              QcMatSetDataType_f, &
                              QcMatSetDimMat_f,   &
                              QcMatAssemble_f,    &
                              QcMatSetValues_f,   &
                              QcMatDestroy_f
        use openrsp_f, only: ELEC_AO_D_MATRIX,             &
                             OpenRSP,                      &
                             OpenRSPSetElecEOM_f,          &
                             OpenRSPSetLinearRSPSolver_f,  &
                             OpenRSPSetPDBS_f,             &
                             OpenRSPAddOneOper_f,          &
                             OpenRSPAddTwoOper_f,          &
                             OpenRSPAssemble_f,            &
                             OpenRSPWrite_f,               &
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
        external get_linear_rsp_solution_f
        ! overlap integrals with London atomic orbitals
        integer(kind=QINT), parameter :: overlap_num_pert = 2_QINT
        integer(kind=QINT) :: overlap_pert_labels(overlap_num_pert) = (/ &
            PERT_GEOMETRIC,PERT_MAGNETIC/)
        integer(kind=QINT) :: overlap_pert_orders(overlap_num_pert) = (/ &
            MAX_ORDER_GEOMETRIC,MAX_ORDER_MAGNETIC/)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: overlap_context(7) = (/"O","V","E","R","L","A","P"/)
#endif
        external get_overlap_mat_f
        external get_overlap_exp_f
        ! one-electron Hamiltonian
        integer(kind=QINT), parameter :: oneham_num_pert = 2_QINT
        integer(kind=QINT) :: oneham_pert_labels(oneham_num_pert) = (/ &
            PERT_GEOMETRIC,PERT_MAGNETIC/)
        integer(kind=QINT) :: oneham_pert_orders(oneham_num_pert) = (/ &
            MAX_ORDER_GEOMETRIC,MAX_ORDER_MAGNETIC/)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: oneham_context(6) = (/"O","N","E","H","A","M"/)
#endif
        ! external field
        integer(kind=QINT), parameter :: ext_field_num_pert = 3_QINT
        integer(kind=QINT) :: ext_field_pert_labels(ext_field_num_pert) = (/ &
            PERT_GEOMETRIC,PERT_DIPOLE,PERT_MAGNETIC/)
        integer(kind=QINT) :: ext_field_pert_orders(ext_field_num_pert) = (/ &
            MAX_ORDER_GEOMETRIC,MAX_ORDER_DIPOLE,MAX_ORDER_MAGNETIC/)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: ext_field_context(9) = (/"E","X","T","_","F","I", "E", "L", "D"/)
#endif
        external get_one_oper_mat_f
        external get_one_oper_exp_f
        ! two-electron Hamiltonian
        integer(kind=QINT), parameter :: twoel_num_pert = 2_QINT
        integer(kind=QINT) :: twoel_pert_labels(twoel_num_pert) = (/ &
            PERT_GEOMETRIC,PERT_MAGNETIC/)
        integer(kind=QINT) :: twoel_pert_orders(twoel_num_pert) = (/ &
            MAX_ORDER_GEOMETRIC,MAX_ORDER_MAGNETIC/)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: twoel_context(6) = (/"N","O","N","L","A","O"/)
#endif
        external get_two_oper_mat_f
        external get_two_oper_exp_f
        ! referenced state
#include "tests/ao_dens/openrsp_f_ao_dims.h90"
#include "tests/ao_dens/openrsp_f_ao_fock.h90"
#include "tests/ao_dens/openrsp_f_ao_density.h90"
#include "tests/ao_dens/openrsp_f_ao_overlap.h90"
        type(QcMat) F_unpert
        type(QcMat) D_unpert
        type(QcMat) S_unpert
        ! polarizability
        integer(kind=QINT), parameter :: ALPHA_NUM_PROPS = 1_QINT
        integer(kind=QINT), parameter :: ALPHA_NUM_PERT(ALPHA_NUM_PROPS) = (/2_QINT/)
        integer(kind=QINT), parameter :: ALPHA_PERT_LABELS(sum(ALPHA_NUM_PERT)) &
            = (/PERT_DIPOLE,PERT_DIPOLE/)
        integer(kind=QINT), parameter :: ALPHA_NUM_FREQS(ALPHA_NUM_PROPS) = (/1_QINT/)
        real(kind=QREAL), parameter :: ALPHA_PERT_FREQS(2*dot_product(ALPHA_NUM_FREQS,ALPHA_NUM_PERT)) &
            = (/-0.072_QREAL,0.0_QREAL,0.072_QREAL,0.0_QREAL/)
        integer(kind=QINT), parameter :: ALPHA_KN_RULES(ALPHA_NUM_PROPS) = (/0_QINT/)
        ! response functions
        integer(kind=QINT) size_rsp_funs
        real(kind=QREAL) rsp_funs(18)
        ! error information
        integer(kind=4) ierr

        ! sets the equation of motion of electrons
        ierr = OpenRSPSetElecEOM_f(open_rsp, ELEC_AO_D_MATRIX)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPSetElecEOM_f() passed"

        ! sets the context of linear response equation solver
        ierr = OpenRSPSetLinearRSPSolver_f(open_rsp,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                                           solver_context, &
#endif
                                           get_linear_rsp_solution_f)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPSetLinearRSPSolver_f() passed"

        ! sets the context of perturbation dependent basis sets
        ierr = OpenRSPSetPDBS_f(open_rsp,            &
                                overlap_num_pert,    &
                                overlap_pert_labels, &
                                overlap_pert_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                overlap_context,     &
#endif
                                get_overlap_mat_f,   &
                                get_overlap_exp_f)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPSetPDBS_f() passed"

        ! adds one-electron Hamiltonian
        ierr = OpenRSPAddOneOper_f(open_rsp,           &
                                   oneham_num_pert,    &
                                   oneham_pert_labels, &
                                   oneham_pert_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                   oneham_context,     &
#endif
                                   get_one_oper_mat_f, &
                                   get_one_oper_exp_f)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPAddOneOper_f(h) passed"

        ! adds external field
        ierr = OpenRSPAddOneOper_f(open_rsp,              &
                                   ext_field_num_pert,    &
                                   ext_field_pert_labels, &
                                   ext_field_pert_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                   ext_field_context,     &
#endif
                                   get_one_oper_mat_f,    &
                                   get_one_oper_exp_f)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPAddOneOper_f(V) passed"

        ! adds two-electron Hamiltonian
        ierr = OpenRSPAddTwoOper_f(open_rsp,           &
                                   twoel_num_pert,     &
                                   twoel_pert_labels,  &
                                   twoel_pert_orders,  &
#if defined(OPENRSP_F_USER_CONTEXT)
                                   twoel_context,      &
#endif
                                   get_two_oper_mat_f, &
                                   get_two_oper_exp_f)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPAddTwoOper_f() passed"

        ! assembles the context of response theory calculations
        ierr = OpenRSPAssemble_f(open_rsp)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPAssemble_f() passed"

        ! writes the context of response theory calculations
        ierr = OpenRSPWrite_f(open_rsp, OPENRSP_F_LOG)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPWrite_f() passed"

        ! sets the unperturbed Fock matrix
        ierr = QcMatCreate_f(A=F_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatBlockCreate_f(A=F_unpert, dim_block=1_QINT)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatSetSymType_f(A=F_unpert, sym_type=QSYMMAT)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatSetDataType_f(A=F_unpert,                      &
                                  num_blocks=1_QINT,               &
                                  idx_block_row=(/IDX_BLOCK_ROW/), &
                                  idx_block_col=(/IDX_BLOCK_COL/), &
                                  data_type=(/QREALMAT/))
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatSetDimMat_f(A=F_unpert, num_row=NUM_AO, num_col=NUM_AO)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatAssemble_f(A=F_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatSetValues_f(A=F_unpert,                  &
                                idx_block_row=IDX_BLOCK_ROW, &
                                idx_block_col=IDX_BLOCK_COL, &
                                idx_first_row=IDX_FIRST_ROW, &
                                num_row_set=NUM_AO,          &
                                idx_first_col=IDX_FIRST_COL, &
                                num_col_set=NUM_AO,          &
                                values_real=values_fock)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ! sets the unperturbed density matrix
        ierr = QcMatCreate_f(A=D_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatBlockCreate_f(A=D_unpert, dim_block=1_QINT)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatSetSymType_f(A=D_unpert, sym_type=QSYMMAT)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatSetDataType_f(A=D_unpert,                      &
                                  num_blocks=1_QINT,               &
                                  idx_block_row=(/IDX_BLOCK_ROW/), &
                                  idx_block_col=(/IDX_BLOCK_COL/), &
                                  data_type=(/QREALMAT/))
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatSetDimMat_f(A=D_unpert, num_row=NUM_AO, num_col=NUM_AO)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatAssemble_f(A=D_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatSetValues_f(A=D_unpert,                  &
                                idx_block_row=IDX_BLOCK_ROW, &
                                idx_block_col=IDX_BLOCK_COL, &
                                idx_first_row=IDX_FIRST_ROW, &
                                num_row_set=NUM_AO,          &
                                idx_first_col=IDX_FIRST_COL, &
                                num_col_set=NUM_AO,          &
                                values_real=values_density)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ! sets the unperturbed overlap integrals
        ierr = QcMatCreate_f(A=S_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatBlockCreate_f(A=S_unpert, dim_block=1_QINT)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatSetSymType_f(A=S_unpert, sym_type=QSYMMAT)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatSetDataType_f(A=S_unpert,                      &
                                  num_blocks=1_QINT,               &
                                  idx_block_row=(/IDX_BLOCK_ROW/), &
                                  idx_block_col=(/IDX_BLOCK_COL/), &
                                  data_type=(/QREALMAT/))
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatSetDimMat_f(A=S_unpert, num_row=NUM_AO, num_col=NUM_AO)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatAssemble_f(A=S_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatSetValues_f(A=S_unpert,                  &
                                idx_block_row=IDX_BLOCK_ROW, &
                                idx_block_col=IDX_BLOCK_COL, &
                                idx_first_row=IDX_FIRST_ROW, &
                                num_row_set=NUM_AO,          &
                                idx_first_col=IDX_FIRST_COL, &
                                num_col_set=NUM_AO,          &
                                values_real=values_overlap)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)

        ! gets the polarizability
        size_rsp_funs = 9_QINT
        ierr = OpenRSPGetRSPFun_f(open_rsp,          &
                                  F_unpert,          &
                                  D_unpert,          &
                                  S_unpert,          &
                                  ALPHA_NUM_PROPS,   &
                                  ALPHA_NUM_PERT,    &
                                  ALPHA_PERT_LABELS, &
                                  ALPHA_NUM_FREQS,   &
                                  ALPHA_PERT_FREQS,  &
                                  ALPHA_KN_RULES,    &
                                  size_rsp_funs,     &
                                  rsp_funs)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPGetRSPFun_f() passed"
        write(io_log,110) rsp_funs(1:2*size_rsp_funs)

        ! cleans
        ierr = QcMatDestroy_f(A=F_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatDestroy_f(A=D_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        ierr = QcMatDestroy_f(A=S_unpert)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)

        return

100     format("test_f_OpenRSP_AO>> ",A)
110     format(3(" (",F16.10,",",F16.10,")"))
    end subroutine test_f_OpenRSP_AO

#undef OPENRSP_F_TEST_SRC
