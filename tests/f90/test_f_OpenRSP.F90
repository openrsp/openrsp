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
!!  This file tests the APIs of OpenRSP library.
!!
!!  2014-08-03, Bin Gao
!!  * first version

! parameters for the test suite
#include "tests/openrsp_f_test.h"

#define OPENRSP_F_TEST_SRC "tests/f90/test_f_OpenRSP.F90"

#if defined(OPENRSP_TEST_EXECUTABLE)
    program test_f_OpenRSP
#else
    subroutine test_f_OpenRSP(io_log)
#endif
        use qmatrix, only: QINT,QREAL,QMat
        use openrsp_f, only: OpenRSP,                 &
                             OpenRSPCreate,           &
                             OpenRSPSetSolver,        &
#if defined(OPENRSP_PERTURBATION_FREE)
                               OpenRSPSetPerturbations, &
#endif
                               OpenRSPSetPDBS,          &
                               OpenRSPAddOneOper,       &
                               OpenRSPAssemble,         &
                               OpenRSPWrite,            &
                               OpenRSPGetRSPFun,        &
                               OpenRSPDestroy
        implicit none
        ! IO of standard output
#if defined(OPENRSP_TEST_EXECUTABLE)
        integer(kind=4), parameter :: io_log = 6
#else
        integer(kind=4), intent(in) :: io_log
#endif
        ! context of response theory calculations
        type(OpenRSP) open_rsp
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: solver_lab(6) = (/"S","O","L","V","E","R"/)
#endif
        external get_rsp_solution
        ! all perturbations involved in calculations
        integer(kind=QINT), parameter :: ALL_PERTURBATIONS(NUM_ALL_PERT) = (/ &
            PERT_DIPOLE,PERT_MAGNETIC,PERT_GEOMETRIC/)
        ! maximum allowed orders of all perturbations
        integer(kind=QINT), parameter :: ALL_PERT_MAX_ORDERS(NUM_ALL_PERT) = (/ &
            MAX_ORDER_DIPOLE,MAX_ORDER_MAGNETIC,MAX_ORDER_GEOMETRIC/)
        ! sizes of all perturbations up to their maximum orders
        integer(kind=QINT), parameter :: ALL_PERT_SIZES(15) = (/ &
            3,                                                   &  !electric dipole
            3,6,10,15,21,28,36,                                  &  !magnetic derivatives
            36,666,8436,82251,465552,1898232,6257232/)              !geometric derivatives (12 atoms)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: pert_lab(4) = (/"P","E","R","T"/)
#endif
#if defined(OPENRSP_PERTURBATION_FREE)
        external get_pert_comp
        external get_pert_rank
#endif
        ! overlap integrals with London atomic orbitals
        integer(kind=QINT), parameter :: overlap_num_pert = 2
        integer(kind=QINT) :: overlap_perturbations(overlap_num_pert) = (/ &
            PERT_MAGNETIC,PERT_GEOMETRIC/)
        integer(kind=QINT) :: overlap_pert_orders(overlap_num_pert) = (/ &
            MAX_ORDER_MAGNETIC,MAX_ORDER_GEOMETRIC/)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: overlap_lab(7) = (/"O","V","E","R","L","A","P"/)
#endif
        external get_overlap_mat
        external get_overlap_exp
        ! one-electron Hamiltonian
        integer(kind=QINT), parameter :: oneham_num_pert = 2
        integer(kind=QINT) :: oneham_perturbations(oneham_num_pert) = (/ &
            PERT_MAGNETIC,PERT_GEOMETRIC/)
        integer(kind=QINT) :: oneham_pert_orders(oneham_num_pert) = (/ &
            MAX_ORDER_MAGNETIC,MAX_ORDER_GEOMETRIC/)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: oneham_lab(6) = (/"O","N","E","H","A","M"/)
#endif
        ! external field
        integer(kind=QINT), parameter :: ext_field_num_pert = 3
        integer(kind=QINT) :: ext_field_perturbations(ext_field_num_pert) = (/ &
            PERT_DIPOLE,PERT_MAGNETIC,PERT_GEOMETRIC/)
        integer(kind=QINT) :: ext_field_pert_orders(ext_field_num_pert) = (/ &
            MAX_ORDER_DIPOLE,MAX_ORDER_MAGNETIC,MAX_ORDER_GEOMETRIC/)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: ext_field_lab(9) = (/"E","X","T","_","F","I", "E", "L", "D"/)
#endif
        external get_one_oper_mat
        external get_one_oper_exp
        ! referenced state
        type(QMat) ref_ham
        type(QMat) ref_state
        type(QMat) ref_overlap
!FIXME: to move to test_...
        integer(kind=QINT), parameter :: vib_alpha_num_pert = 2
        integer(kind=QINT) :: vib_alpha_perturbations(vib_alpha_num_pert) = (/PERT_DIPOLE,PERT_GEOMETRIC/)
        integer(kind=QINT) :: vib_alpha_pert_orders(vib_alpha_num_pert) = (/1,1/)
        real(kind=QREAL) :: vib_alpha_pert_freqs(vib_alpha_num_pert) = (/0.1,0.0/)
        integer(kind=QINT) :: kn_rule(3) = (/0,1,0/)
        integer(kind=QINT) :: size_rsp_fun = 9
        real(kind=QREAL) :: rsp_fun(9)
 
        ! error information
        integer(kind=4) ierr

        ierr = OpenRSPCreate(open_rsp)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPCreate() passed"

        ierr = OpenRSPSetSolver(open_rsp,   &
#if defined(OPENRSP_F_USER_CONTEXT)
                                solver_lab, &
#endif
                                get_rsp_solution)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPSetSolver() passed"

#if defined(OPENRSP_PERTURBATION_FREE)
        ierr = OpenRSPSetPerturbations(open_rsp,            &
                                       NUM_ALL_PERT,        &
                                       ALL_PERTURBATIONS,   &
                                       ALL_PERT_MAX_ORDERS, &
                                       ALL_PERT_SIZES,      &
#if defined(OPENRSP_F_USER_CONTEXT)
                                       pert_lab,            &
#endif
                                       get_pert_comp,       &
                                       get_pert_rank)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPSetPerturbations() passed"
#endif

        ierr = OpenRSPSetPDBS(open_rsp,              &
                              overlap_num_pert,      &
                              overlap_perturbations, &
                              overlap_pert_orders,   &
#if defined(OPENRSP_F_USER_CONTEXT)
                              overlap_lab,           &
#endif
                              get_overlap_mat,       &
                              get_overlap_exp)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPSetPDBS passed"

        ierr = OpenRSPAddOneOper(open_rsp,             &
                                 oneham_num_pert,      &
                                 oneham_perturbations, &
                                 oneham_pert_orders,   &
#if defined(OPENRSP_F_USER_CONTEXT)
                                 oneham_lab,           &
#endif
                                 get_one_oper_mat,     &
                                 get_one_oper_exp)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPAddOneOper(h) passed"

        ierr = OpenRSPAddOneOper(open_rsp,                &
                                 ext_field_num_pert,      &
                                 ext_field_perturbations, &
                                 ext_field_pert_orders,   &
#if defined(OPENRSP_F_USER_CONTEXT)
                                 ext_field_lab,           &
#endif
                                 get_one_oper_mat,        &
                                 get_one_oper_exp)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPAddOneOper(V) passed"

        ierr = OpenRSPAssemble(open_rsp)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPAssemble() passed"

        ierr = OpenRSPWrite(open_rsp, OPENRSP_F_LOG)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPWrite() passed"

        ierr = OpenRSPGetRSPFun(open_rsp,                &
                                ref_ham,                 &
                                ref_state,               &
                                ref_overlap,             &
                                vib_alpha_num_pert,      &
                                vib_alpha_perturbations, &
                                vib_alpha_pert_orders,   &
                                vib_alpha_pert_freqs,    &
                                kn_rule,                 &
                                size_rsp_fun,            &
                                rsp_fun)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPGetRSPFun() passed"

        ierr = OpenRSPDestroy(open_rsp)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPDestroy() passed"

100     format("test_f_OpenRSP>> ",A)
#if defined(OPENRSP_TEST_EXECUTABLE)
    end program test_f_OpenRSP
#else
    end subroutine test_f_OpenRSP
#endif

#undef OPENRSP_F_TEST_SRC
