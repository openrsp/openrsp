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
!!  This file contains the Fortran APIs of OpenRSP library.
!!
!!  2014-07-12, Bin Gao
!!  * first version

! basic data types
#include "api/qmatrix_c_type.h"

module openrsp_mod

    use, intrinsic :: iso_c_binding
    use qmatrix, only: QINT,  &
                       QREAL, &
                       QMat
    use rsp_solver_f, only: SolverFun_f,       &
                            RSPSolverCreate_f, &
                            RSPSolverDestroy_f
#if defined(OPENRSP_PERTURBATION_FREE)
    use rsp_pert_f, only: PertFun_f,       &
                          RSPPertCreate_f, &
                          RSPPertDestroy_f
#endif
    use rsp_overlap_f, only: OverlapFun_f,       &
                             RSPOverlapCreate_f, &
                             RSPOverlapDestroy_f
    use rsp_one_oper_f, only: OneOperFun_f,       &
                              RSPOneOperCreate_f, &
                              RSPOneOperDestroy_f
    use rsp_two_oper_f, only: TwoOperFun_f,       &
                              RSPTwoOperCreate_f, &
                              RSPTwoOperDestroy_f

    implicit none

    ! type of equation of motion (EOM) of electrons
    integer(kind=QINT), parameter, public :: ELEC_EOM_DMAT = 0
    integer(kind=QINT), parameter, public :: ELEC_EOM_CMAT = 1
    integer(kind=QINT), parameter, public :: ELEC_EOM_CC = 2

    ! linked list of context of callback subroutines of one-electron operators
    type, private :: OneOperList_f
        type(OneOperFun_f), pointer :: one_oper_fun => null()
        type(OneOperList_f), pointer :: next_one_oper => null()
    end type OneOperList_f

    ! linked list of context of callback subroutines of two-electron operators
    type, private :: TwoOperList_f
        type(TwoOperFun_f), pointer :: two_oper_fun => null()
        type(TwoOperList_f), pointer :: next_two_oper => null()
    end type TwoOperList_f

    ! OpenRSP type (inspired by http://wiki.rac.manchester.ac.uk/community/GPU/GpuFaq/FortranGPU)
    type, public :: OpenRSP
        private
        type(C_PTR) :: c_rsp = C_NULL_PTR
        type(SolverFun_f), pointer :: solver_fun => null()
#if defined(OPENRSP_PERTURBATION_FREE)
        type(PertFun_f), pointer :: pert_fun => null()
#endif
        type(OverlapFun_f), pointer :: overlap_fun => null()
        type(OneOperList_f), pointer :: list_one_oper => null()
        type(TwoOperList_f), pointer :: list_two_oper => null()
    end type OpenRSP

    ! functions provided by the Fortran APIs
    public :: OpenRSPCreate
    public :: OpenRSPSetElecEOM
    public :: OpenRSPSetSolver
#if defined(OPENRSP_PERTURBATION_FREE)
    public :: OpenRSPSetPerturbations
#endif
    public :: OpenRSPSetPDBS
    public :: OpenRSPAddOneOper
    public :: OpenRSPAddTwoOper
    !public :: OpenRSPAddXCFun
    !public :: OpenRSPAddNucFun
    public :: OpenRSPAssemble
    public :: OpenRSPWrite
    public :: OpenRSPGetRSPFun
    !public :: OpenRSPGetResidue
    public :: OpenRSPDestroy

    interface 
        integer(C_INT) function f03_api_OpenRSPCreate(open_rsp) &
            bind(C, name="f03_api_OpenRSPCreate")
            use, intrinsic :: iso_c_binding
            type(C_PTR), intent(inout) :: open_rsp
        end function f03_api_OpenRSPCreate
        integer(C_INT) function f03_api_OpenRSPSetElecEOM(open_rsp,      &
                                                          elec_eom_type) &
            bind(C, name="f03_api_OpenRSPSetElecEOM")
            use, intrinsic :: iso_c_binding
            type(C_PTR), intent(inout) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: elec_eom_type
        end function f03_api_OpenRSPSetElecEOM
        integer(C_INT) function f03_api_OpenRSPSetSolver(open_rsp,         &
                                                         user_ctx,         &
                                                         get_rsp_solution) &
            bind(C, name="f03_api_OpenRSPSetSolver")
            use, intrinsic :: iso_c_binding
            type(C_PTR), intent(inout) :: open_rsp
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_rsp_solution
        end function f03_api_OpenRSPSetSolver
#if defined(OPENRSP_PERTURBATION_FREE)
        integer(C_INT) function f03_api_OpenRSPSetPerturbations(open_rsp,        &
                                                                num_pert,        &
                                                                perturbations,   &
                                                                pert_max_orders, &
                                                                pert_sizes,      &
                                                                user_ctx,        &
                                                                get_pert_comp,   &
                                                                get_pert_rank)   &
            bind(C, name="f03_api_OpenRSPSetPerturbations")
            use, intrinsic :: iso_c_binding
            type(C_PTR), intent(inout) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: perturbations(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_sizes(sum(pert_max_orders))
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_pert_comp
            type(C_FUNPTR), value, intent(in) :: get_pert_rank
        end function f03_api_OpenRSPSetPerturbations
#endif
        integer(C_INT) function f03_api_OpenRSPSetPDBS(open_rsp,        &
                                                       num_pert,        &
                                                       perturbations,   &
                                                       pert_max_orders, &
                                                       user_ctx,        &
                                                       get_overlap_mat, &
                                                       get_overlap_exp) &
            bind(C, name="f03_api_OpenRSPSetPDBS")
            use, intrinsic :: iso_c_binding
            type(C_PTR), intent(inout) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: perturbations(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert)
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_overlap_mat
            type(C_FUNPTR), value, intent(in) :: get_overlap_exp
        end function f03_api_OpenRSPSetPDBS
        integer(C_INT) function f03_api_OpenRSPAddOneOper(open_rsp,         &
                                                          num_pert,         &
                                                          perturbations,    &
                                                          pert_max_orders,  &
                                                          user_ctx,         &
                                                          get_one_oper_mat, &
                                                          get_one_oper_exp) &
            bind(C, name="f03_api_OpenRSPAddOneOper")
            use, intrinsic :: iso_c_binding  
            type(C_PTR), intent(inout) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: perturbations(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert)
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_one_oper_mat
            type(C_FUNPTR), value, intent(in) :: get_one_oper_exp
        end function f03_api_OpenRSPAddOneOper
        integer(C_INT) function f03_api_OpenRSPAddTwoOper(open_rsp,         &
                                                          num_pert,         &
                                                          perturbations,    &
                                                          pert_max_orders,  &
                                                          user_ctx,         &
                                                          get_two_oper_mat, &
                                                          get_two_oper_exp) &
            bind(C, name="f03_api_OpenRSPAddTwoOper")
            use, intrinsic :: iso_c_binding
            type(C_PTR), intent(inout) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: perturbations(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert)
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_two_oper_mat
            type(C_FUNPTR), value, intent(in) :: get_two_oper_exp
        end function f03_api_OpenRSPAddTwoOper
        integer(C_INT) function f03_api_OpenRSPAssemble(open_rsp) &
            bind(C, name="f03_api_OpenRSPAssemble")
            use, intrinsic :: iso_c_binding
            type(C_PTR), intent(inout) :: open_rsp
        end function f03_api_OpenRSPAssemble
        integer(C_INT) function f03_api_OpenRSPWrite(open_rsp, file_name) &
            bind(C, name="f03_api_OpenRSPWrite")
            use, intrinsic :: iso_c_binding
            type(C_PTR), intent(in) :: open_rsp
            character(C_CHAR), intent(in) :: file_name(*)
        end function f03_api_OpenRSPWrite
        integer(C_INT) function f03_api_OpenRSPGetRSPFun(open_rsp,      &
                                                         ref_ham,       &
                                                         ref_state,     &
                                                         ref_overlap,   &
                                                         num_pert,      &
                                                         perturbations, &
                                                         pert_orders,   &
                                                         pert_freqs,    &
                                                         kn_rule,       &
                                                         size_rsp_fun,  &
                                                         rsp_fun)       &
            bind(C, name="f03_api_OpenRSPGetRSPFun")
            use, intrinsic :: iso_c_binding
            type(C_PTR), intent(inout) :: open_rsp
            type(C_PTR), value, intent(in) :: ref_ham
            type(C_PTR), value, intent(in) :: ref_state
            type(C_PTR), value, intent(in) :: ref_overlap
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: perturbations(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
            real(kind=C_QREAL), intent(in) :: pert_freqs(num_pert)
            integer(kind=C_QINT), intent(in) :: kn_rule(3)
            integer(kind=C_QINT), value, intent(in) :: size_rsp_fun
            real(kind=C_QREAL), intent(out) :: rsp_fun(size_rsp_fun)
        end function f03_api_OpenRSPGetRSPFun
        integer(C_INT) function f03_api_OpenRSPDestroy(open_rsp) &
            bind(C, name="f03_api_OpenRSPDestroy")
            use, intrinsic :: iso_c_binding
            type(C_PTR), intent(inout) :: open_rsp
        end function f03_api_OpenRSPDestroy
    end interface

    contains

    function OpenRSPCreate(open_rsp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        ierr = f03_api_OpenRSPCreate(open_rsp%c_rsp)
        nullify(open_rsp%solver_fun)
#if defined(OPENRSP_PERTURBATION_FREE)
        nullify(open_rsp%pert_fun)
#endif
        nullify(open_rsp%overlap_fun)
        nullify(open_rsp%list_one_oper)
        nullify(open_rsp%list_two_oper)
    end function OpenRSPCreate

    function OpenRSPSetElecEOM(open_rsp, elec_eom_type) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: elec_eom_type
        ierr = f03_api_OpenRSPSetElecEOM(open_rsp%c_rsp, elec_eom_type)
    end function OpenRSPSetElecEOM

    function OpenRSPSetSolver(open_rsp,        &
#if defined(OPENRSP_F_USER_CONTEXT)
                              user_ctx,        &
#endif
                              get_rsp_solution) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_rsp_solution(ref_ham,       &
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
                use qmatrix, only: QINT,QREAL,QMat
                type(QMat), intent(in) :: ref_ham
                type(QMat), intent(in) :: ref_state
                type(QMat), intent(in) :: ref_overlap
                integer(kind=QINT), intent(in) :: num_freq_sums
                real(kind=QREAL), intent(in) :: freq_sums(num_freq_sums)
                integer(kind=QINT), intent(in) :: size_pert
                type(QMat), intent(in) :: RHS_mat(size_pert*num_freq_sums)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                type(QMat), intent(inout) :: rsp_param(size_pert*num_freq_sums)
            end subroutine get_rsp_solution
            subroutine RSPSolverGetSolution_f(ref_ham,       &
                                              ref_state,     &
                                              ref_overlap,   &
                                              num_freq_sums, &
                                              freq_sums,     &
                                              size_pert,     &
                                              RHS_mat,       &
                                              user_ctx,      &
                                              rsp_param)     &
                bind(C, name="RSPSolverGetSolution_f")
                use, intrinsic :: iso_c_binding
                type(C_PTR), intent(in) :: ref_ham
                type(C_PTR), intent(in) :: ref_state
                type(C_PTR), intent(in) :: ref_overlap
                integer(kind=C_QINT), value, intent(in) :: num_freq_sums
                real(kind=C_QREAL), intent(in) :: freq_sums(num_freq_sums)
                integer(kind=C_QINT), value, intent(in) :: size_pert
                type(C_PTR), intent(in) :: RHS_mat(size_pert*num_freq_sums)
                type(C_PTR), value, intent(in) :: user_ctx
                type(C_PTR), intent(inout) :: rsp_param(size_pert*num_freq_sums)
            end subroutine RSPSolverGetSolution_f
        end interface
        if (associated(open_rsp%solver_fun)) then
            call RSPSolverDestroy_f(open_rsp%solver_fun)
        else
            allocate(open_rsp%solver_fun)
        end if
        ! adds context of callback function of response equation solver
        call RSPSolverCreate_f(open_rsp%solver_fun, &
#if defined(OPENRSP_F_USER_CONTEXT)
                               user_ctx,            &
#endif
                               get_rsp_solution)
        ierr = f03_api_OpenRSPSetSolver(open_rsp%c_rsp,             &
                                        c_loc(open_rsp%solver_fun), &
                                        c_funloc(RSPSolverGetSolution_f))
    end function OpenRSPSetSolver

#if defined(OPENRSP_PERTURBATION_FREE)
    function OpenRSPSetPerturbations(open_rsp,        &
                                     num_pert,        &
                                     perturbations,   &
                                     pert_max_orders, &
                                     pert_sizes,      &
#if defined(OPENRSP_F_USER_CONTEXT)
                                     user_ctx,        &
#endif
                                     get_pert_comp,   &
                                     get_pert_rank) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: perturbations(num_pert)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert)
        integer(kind=QINT), intent(in) :: pert_sizes(:)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
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
                integer(kind=QINT), intent(in) :: perturbation
                integer(kind=QINT), intent(in) :: pert_num_comp
                integer(kind=QINT), intent(in) :: pert_components(pert_num_comp)
                integer(kind=QINT), intent(in) :: pert_comp_orders(pert_num_comp)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(out) :: pert_rank
            end subroutine get_pert_rank
            subroutine RSPPertGetComp_f(perturbation,     &
                                        pert_order,       &
                                        pert_rank,        &
                                        user_ctx,         &
                                        pert_num_comp,    &
                                        pert_components,  &
                                        pert_comp_orders) &
                bind(C, name="RSPPertGetComp_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: perturbation
                integer(kind=C_QINT), value, intent(in) :: pert_order
                integer(kind=C_QINT), value, intent(in) :: pert_rank
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), intent(out) :: pert_num_comp
                integer(kind=C_QINT), intent(out) :: pert_components(pert_order)
                integer(kind=C_QINT), intent(out) :: pert_comp_orders(pert_order)
            end subroutine RSPPertGetComp_f
            subroutine RSPPertGetRank_f(perturbation,     &
                                        pert_num_comp,    &
                                        pert_components,  &
                                        pert_comp_orders, &
                                        user_ctx,         &
                                        pert_rank)        &
                bind(C, name="RSPPertGetRank_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: perturbation
                integer(kind=C_QINT), value, intent(in) :: pert_num_comp
                integer(kind=C_QINT), intent(in) :: pert_components(pert_num_comp)
                integer(kind=C_QINT), intent(in) :: pert_comp_orders(pert_num_comp)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), intent(out) :: pert_rank
            end subroutine RSPPertGetRank_f
        end interface
        if (associated(open_rsp%pert_fun)) then
            call RSPPertDestroy_f(open_rsp%pert_fun)
        else
            allocate(open_rsp%pert_fun)
        end if
        ! adds context of callback functions of perturbations
        call RSPPertCreate_f(open_rsp%pert_fun, &
#if defined(OPENRSP_F_USER_CONTEXT)
                             user_ctx,          &
#endif
                             get_pert_comp,     &
                             get_pert_rank)
        ierr = f03_api_OpenRSPSetPerturbations(open_rsp%c_rsp,             &
                                               num_pert,                   &
                                               perturbations,              &
                                               pert_max_orders,            &
                                               pert_sizes,                 &
                                               c_loc(open_rsp%pert_fun),   &
                                               c_funloc(RSPPertGetComp_f), &
                                               c_funloc(RSPPertGetRank_f))
    end function OpenRSPSetPerturbations
#endif

    function OpenRSPSetPDBS(open_rsp,        &
                            num_pert,        &
                            perturbations,   &
                            pert_max_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                            user_ctx,        &
#endif
                            get_overlap_mat, &
                            get_overlap_exp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: perturbations(num_pert)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
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
                use qmatrix, only: QINT,QREAL,QMat
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
            end subroutine get_overlap_exp
            subroutine RSPOverlapGetMat_f(bra_num_pert,      &
                                          bra_perturbations, &
                                          bra_pert_orders,   &
                                          ket_num_pert,      &
                                          ket_perturbations, &
                                          ket_pert_orders,   &
                                          num_pert,          &
                                          perturbations,     &
                                          pert_orders,       &
                                          user_ctx,          &
                                          num_int,           &
                                          val_int)           &
                bind(C, name="RSPOverlapGetMat_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: bra_num_pert
                integer(kind=C_QINT), intent(in) :: bra_perturbations(bra_num_pert)
                integer(kind=C_QINT), intent(in) :: bra_pert_orders(bra_num_pert)
                integer(kind=C_QINT), value, intent(in) :: ket_num_pert
                integer(kind=C_QINT), intent(in) :: ket_perturbations(ket_num_pert)
                integer(kind=C_QINT), intent(in) :: ket_pert_orders(ket_num_pert)
                integer(kind=C_QINT), value, intent(in) :: num_pert
                integer(kind=C_QINT), intent(in) :: perturbations(num_pert)
                integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_int
                type(C_PTR), intent(inout) :: val_int(num_int)
            end subroutine RSPOverlapGetMat_f
            subroutine RSPOverlapGetExp_f(bra_num_pert,      &
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
                                          user_ctx,          &
                                          num_exp,           &
                                          val_exp)           &
                bind(C, name="RSPOverlapGetExp_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: bra_num_pert
                integer(kind=C_QINT), intent(in) :: bra_perturbations(bra_num_pert)
                integer(kind=C_QINT), intent(in) :: bra_pert_orders(bra_num_pert)
                integer(kind=C_QINT), value, intent(in) :: ket_num_pert
                integer(kind=C_QINT), intent(in) :: ket_perturbations(ket_num_pert)
                integer(kind=C_QINT), intent(in) :: ket_pert_orders(ket_num_pert)
                integer(kind=C_QINT), value, intent(in) :: num_pert
                integer(kind=C_QINT), intent(in) :: perturbations(num_pert)
                integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
                integer(kind=C_QINT), value, intent(in) :: num_dens
                type(C_PTR), intent(in) :: ao_dens(num_dens)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_exp
                real(kind=C_QREAL), intent(inout) :: val_exp(num_exp)
            end subroutine RSPOverlapGetExp_f
        end interface
        if (associated(open_rsp%overlap_fun)) then
            call RSPOverlapDestroy_f(open_rsp%overlap_fun)
        else
            allocate(open_rsp%overlap_fun)
        end if
        ! adds context of callback functions of the overlap integrals
        call RSPOverlapCreate_f(open_rsp%overlap_fun, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                user_ctx,             &
#endif
                                get_overlap_mat,      &
                                get_overlap_exp)
        ierr = f03_api_OpenRSPSetPDBS(open_rsp%c_rsp,               &
                                      num_pert,                     &
                                      perturbations,                &
                                      pert_max_orders,              &
                                      c_loc(open_rsp%overlap_fun),  &
                                      c_funloc(RSPOverlapGetMat_f), &
                                      c_funloc(RSPOverlapGetExp_f))
    end function OpenRSPSetPDBS

    function OpenRSPAddOneOper(open_rsp,         &
                               num_pert,         &
                               perturbations,    &
                               pert_max_orders,  &
#if defined(OPENRSP_F_USER_CONTEXT)
                               user_ctx,         &
#endif
                               get_one_oper_mat, &
                               get_one_oper_exp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: perturbations(num_pert)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_one_oper_mat(num_pert,      &
                                        perturbations, &
                                        pert_orders,   &
#if defined(OPENRSP_F_USER_CONTEXT)
                                        len_ctx,       &
                                        user_ctx,      &
#endif
                                        num_int,       &
                                        val_int)
                use qmatrix, only: QINT,QREAL,QMat
                integer(kind=QINT), intent(in) :: num_pert
                integer(kind=QINT), intent(in) :: perturbations(num_pert)
                integer(kind=QINT), intent(in) :: pert_orders(num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_int
                type(QMat), intent(inout) :: val_int(num_int)
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
            end subroutine get_one_oper_exp
            subroutine RSPOneOperGetMat_f(num_pert,      &
                                          perturbations, &
                                          pert_orders,   &
                                          user_ctx,      &
                                          num_int,       &
                                          val_int)       &
                bind(C, name="RSPOneOperGetMat_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: num_pert
                integer(kind=C_QINT), intent(in) :: perturbations(num_pert)
                integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_int
                type(C_PTR), intent(inout) :: val_int(num_int)
            end subroutine RSPOneOperGetMat_f
            subroutine RSPOneOperGetExp_f(num_pert,      &
                                          perturbations, &
                                          pert_orders,   &
                                          num_dens,      &
                                          ao_dens,       &
                                          user_ctx,      &
                                          num_exp,       &
                                          val_exp)       &
                bind(C, name="RSPOneOperGetExp_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: num_pert
                integer(kind=C_QINT), intent(in) :: perturbations(num_pert)
                integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
                integer(kind=C_QINT), value, intent(in) :: num_dens
                type(C_PTR), intent(in) :: ao_dens(num_dens)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_exp
                real(kind=C_QREAL), intent(inout) :: val_exp(num_exp)
            end subroutine RSPOneOperGetExp_f
        end interface
        type(OneOperList_f), pointer :: cur_one_oper  !current one-electron operator
        ! inserts the context of callback functions to the tail of the linked list
        if (associated(open_rsp%list_one_oper)) then
            cur_one_oper => open_rsp%list_one_oper
            do while (associated(cur_one_oper%next_one_oper))
                cur_one_oper => cur_one_oper%next_one_oper
            end do
            allocate(cur_one_oper%next_one_oper)
            cur_one_oper => cur_one_oper%next_one_oper
        else
            allocate(open_rsp%list_one_oper)
            cur_one_oper => open_rsp%list_one_oper
        end if
        allocate(cur_one_oper%one_oper_fun)
        nullify(cur_one_oper%next_one_oper)
        ! adds context of callback functions of the new one-electron operator
        call RSPOneOperCreate_f(cur_one_oper%one_oper_fun, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                user_ctx,                  &
#endif
                                get_one_oper_mat,          &
                                get_one_oper_exp)
        ierr = f03_api_OpenRSPAddOneOper(open_rsp%c_rsp,                   &
                                         num_pert,                         &
                                         perturbations,                    &
                                         pert_max_orders,                  &
                                         c_loc(cur_one_oper%one_oper_fun), &
                                         c_funloc(RSPOneOperGetMat_f),     &
                                         c_funloc(RSPOneOperGetExp_f))
    end function OpenRSPAddOneOper

    function OpenRSPAddTwoOper(open_rsp,         &
                               num_pert,         &
                               perturbations,    &
                               pert_max_orders,  &
#if defined(OPENRSP_F_USER_CONTEXT)
                               user_ctx,         &
#endif
                               get_two_oper_mat, &
                               get_two_oper_exp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: perturbations(num_pert)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_two_oper_mat(num_pert,      &
                                        perturbations, &
                                        pert_orders,   &
                                        num_var_dens,  &
                                        var_ao_dens,   &
#if defined(OPENRSP_F_USER_CONTEXT)
                                        len_ctx,       &
                                        user_ctx,      &
#endif
                                        num_int,       &
                                        val_int)
                use qmatrix, only: QINT,QREAL,QMat
                integer(kind=QINT), intent(in) :: num_pert
                integer(kind=QINT), intent(in) :: perturbations(num_pert)
                integer(kind=QINT), intent(in) :: pert_orders(num_pert)
                integer(kind=QINT), intent(in) :: num_var_dens
                type(QMat), intent(in) :: var_ao_dens(num_var_dens)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_int
                type(QMat), intent(inout) :: val_int(num_int)
            end subroutine get_two_oper_mat
            subroutine get_two_oper_exp(num_pert,       &
                                        perturbations,  &
                                        pert_orders,    &
                                        num_var_dens,   &
                                        var_ao_dens,    &
                                        num_contr_dens, &
                                        contr_ao_dens,  &
#if defined(OPENRSP_F_USER_CONTEXT)
                                        len_ctx,        &
                                        user_ctx,       &
#endif
                                        num_exp,        &
                                        val_exp)
                use qmatrix, only: QINT,QREAL,QMat
                integer(kind=QINT), intent(in) :: num_pert
                integer(kind=QINT), intent(in) :: perturbations(num_pert)
                integer(kind=QINT), intent(in) :: pert_orders(num_pert)
                integer(kind=QINT), intent(in) :: num_var_dens
                type(QMat), intent(in) :: var_ao_dens(num_var_dens)
                integer(kind=QINT), intent(in) :: num_contr_dens
                type(QMat), intent(in) :: contr_ao_dens(num_contr_dens)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_exp
                real(kind=QREAL), intent(inout) :: val_exp(num_exp)
            end subroutine get_two_oper_exp
            subroutine RSPTwoOperGetMat_f(num_pert,      &
                                          perturbations, &
                                          pert_orders,   &
                                          num_var_dens,  &
                                          var_ao_dens,   &
                                          user_ctx,      &
                                          num_int,       &
                                          val_int)       &
                bind(C, name="RSPTwoOperGetMat_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: num_pert
                integer(kind=C_QINT), intent(in) :: perturbations(num_pert)
                integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
                integer(kind=C_QINT), value, intent(in) :: num_var_dens
                type(C_PTR), intent(in) :: var_ao_dens(num_var_dens)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_int
                type(C_PTR), intent(inout) :: val_int(num_int)
            end subroutine RSPTwoOperGetMat_f
            subroutine RSPTwoOperGetExp_f(num_pert,       &
                                          perturbations,  &
                                          pert_orders,    &
                                          num_var_dens,   &
                                          var_ao_dens,    &
                                          num_contr_dens, &
                                          contr_ao_dens,  &
                                          user_ctx,       &
                                          num_exp,        &
                                          val_exp)        &
                bind(C, name="RSPTwoOperGetExp_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: num_pert
                integer(kind=C_QINT), intent(in) :: perturbations(num_pert)
                integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
                integer(kind=C_QINT), value, intent(in) :: num_var_dens
                type(C_PTR), intent(in) :: var_ao_dens(num_var_dens)
                integer(kind=C_QINT), value, intent(in) :: num_contr_dens
                type(C_PTR), intent(in) :: contr_ao_dens(num_contr_dens)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_exp
                real(kind=C_QREAL), intent(inout) :: val_exp(num_exp)
            end subroutine RSPTwoOperGetExp_f
        end interface
        type(TwoOperList_f), pointer :: cur_two_oper  !current two-electron operator
        ! inserts the context of callback functions to the tail of the linked list
        if (associated(open_rsp%list_two_oper)) then
            cur_two_oper => open_rsp%list_two_oper
            do while (associated(cur_two_oper%next_two_oper))
                cur_two_oper => cur_two_oper%next_two_oper
            end do
            allocate(cur_two_oper%next_two_oper)
            cur_two_oper => cur_two_oper%next_two_oper
        else
            allocate(open_rsp%list_two_oper)
            cur_two_oper => open_rsp%list_two_oper
        end if
        allocate(cur_two_oper%two_oper_fun)
        nullify(cur_two_oper%next_two_oper)
        ! adds context of callback functions of the new two-electron operator
        call RSPTwoOperCreate_f(cur_two_oper%two_oper_fun, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                user_ctx,                  &
#endif
                                get_two_oper_mat,          &
                                get_two_oper_exp)
        ierr = f03_api_OpenRSPAddTwoOper(open_rsp%c_rsp,                   &
                                         num_pert,                         &
                                         perturbations,                    &
                                         pert_max_orders,                  &
                                         c_loc(cur_two_oper%two_oper_fun), &
                                         c_funloc(RSPTwoOperGetMat_f),     &
                                         c_funloc(RSPTwoOperGetExp_f))
    end function OpenRSPAddTwoOper

    function OpenRSPAssemble(open_rsp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        ierr = f03_api_OpenRSPAssemble(open_rsp%c_rsp)
    end function OpenRSPAssemble

    function OpenRSPWrite(open_rsp, file_name) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(in) :: open_rsp
        character*(*), intent(in) :: file_name
        ierr = f03_api_OpenRSPWrite(open_rsp%c_rsp, file_name//C_NULL_CHAR)
    end function OpenRSPWrite

    function OpenRSPGetRSPFun(open_rsp,      &
                              ref_ham,       &
                              ref_state,     &
                              ref_overlap,   &
                              num_pert,      &
                              perturbations, &
                              pert_orders,   &
                              pert_freqs,    &
                              kn_rule,       &
                              size_rsp_fun,  &
                              rsp_fun) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        type(QMat), target, intent(in) :: ref_ham
        type(QMat), target, intent(in) :: ref_state
        type(QMat), target, intent(in) :: ref_overlap
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: perturbations(num_pert)
        integer(kind=QINT), intent(in) :: pert_orders(num_pert)
        real(kind=QREAL), intent(in) :: pert_freqs(num_pert)
        integer(kind=QINT), intent(in) :: kn_rule(3)
        integer(kind=QINT), intent(in) :: size_rsp_fun
        real(kind=QREAL), intent(out) :: rsp_fun(size_rsp_fun)
        ierr = f03_api_OpenRSPGetRSPFun(open_rsp%c_rsp,     &
                                        c_loc(ref_ham),     &
                                        c_loc(ref_state),   &
                                        c_loc(ref_overlap), &
                                        num_pert,           &
                                        perturbations,      &
                                        pert_orders,        &
                                        pert_freqs,         &
                                        kn_rule,            &
                                        size_rsp_fun,       &
                                        rsp_fun)
    end function OpenRSPGetRSPFun

    function OpenRSPDestroy(open_rsp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        type(OneOperList_f), pointer :: cur_one_oper   !current one-electron operator
        type(OneOperList_f), pointer :: next_one_oper  !next one-electron operator
        type(TwoOperList_f), pointer :: cur_two_oper   !current two-electron operator
        type(TwoOperList_f), pointer :: next_two_oper  !next two-electron operator
        ierr = f03_api_OpenRSPDestroy(open_rsp%c_rsp)
        ! cleans up callback subroutine of response equation solver
        if (associated(open_rsp%solver_fun)) then
            call RSPSolverDestroy_f(open_rsp%solver_fun)
            deallocate(open_rsp%solver_fun)
            nullify(open_rsp%solver_fun)
        end if
#if defined(OPENRSP_PERTURBATION_FREE)
        ! cleans up callback subroutines of perturbations
        if (associated(open_rsp%pert_fun)) then
            call RSPPertDestroy_f(open_rsp%pert_fun)
            deallocate(open_rsp%pert_fun)
            nullify(open_rsp%pert_fun)
        end if
#endif
        ! cleans up callback subroutines of overlap integrals
        if (associated(open_rsp%overlap_fun)) then
            call RSPOverlapDestroy_f(open_rsp%overlap_fun)
            deallocate(open_rsp%overlap_fun)
            nullify(open_rsp%overlap_fun)
        end if
        ! cleans up the linked list of context of callback subroutines of one-electron operators
        cur_one_oper => open_rsp%list_one_oper
        do while (associated(cur_one_oper))
            next_one_oper => cur_one_oper%next_one_oper
            if (associated(cur_one_oper%one_oper_fun)) then
                call RSPOneOperDestroy_f(cur_one_oper%one_oper_fun)
                deallocate(cur_one_oper%one_oper_fun)
                nullify(cur_one_oper%one_oper_fun)
            end if
            deallocate(cur_one_oper)
            nullify(cur_one_oper)
            cur_one_oper => next_one_oper
        end do
        ! cleans up the linked list of context of callback subroutines of two-electron operators
        cur_two_oper => open_rsp%list_two_oper
        do while (associated(cur_two_oper))
            next_two_oper => cur_two_oper%next_two_oper
            if (associated(cur_two_oper%two_oper_fun)) then
                call RSPTwoOperDestroy_f(cur_two_oper%two_oper_fun)
                deallocate(cur_two_oper%two_oper_fun)
                nullify(cur_two_oper%two_oper_fun)
            end if
            deallocate(cur_two_oper)
            nullify(cur_two_oper)
            cur_two_oper => next_two_oper
        end do
    end function OpenRSPDestroy

end module openrsp_mod
