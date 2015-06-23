!!  OpenRSP: open-ended library for response theory
!!  Copyright 2015 Radovan Bast,
!!                 Daniel H. Friese,
!!                 Bin Gao,
!!                 Dan J. Jonsson,
!!                 Magnus Ringholm,
!!                 Kenneth Ruud,
!!                 Andreas Thorvaldsen
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
#include "api/qcmatrix_c_type.h"

module openrsp_f

    use, intrinsic :: iso_c_binding
    use qcmatrix_f, only: QINT,     &
                          QREAL,    &
                          QFAILURE, &
                          QcMat,    &
                          QcMat_C_LOC
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
    use rsp_xc_fun_f, only: XCFunFun_f,       &
                            RSPXCFunCreate_f, &
                            RSPXCFunDestroy_f
    use rsp_nuc_contrib_f, only: NucHamiltonFun_f,       &
                                 RSPNucHamiltonCreate_f, &
                                 RSPNucHamiltonDestroy_f

    implicit none

    ! type of equation of motion (EOM) of electrons
    integer(kind=QINT), parameter, public :: ELEC_AO_D_MATRIX = 0
    integer(kind=QINT), parameter, public :: ELEC_MO_C_MATRIX = 1
    integer(kind=QINT), parameter, public :: ELEC_COUPLED_CLUSTER = 2

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

    ! linked list of context of callback subroutines of XC functionals
    type, private :: XCFunList_f
        type(XCFunFun_f), pointer :: xc_fun_fun => null()
        type(XCFunList_f), pointer :: next_xc_fun => null()
    end type XCFunList_f

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
        type(XCFunList_f), pointer :: list_xc_fun => null()
        type(NucHamiltonFun_f), pointer :: nuc_contrib_fun => null()
    end type OpenRSP

    ! functions provided by the Fortran APIs
    public :: OpenRSPCreate_f
    public :: OpenRSPSetElecEOM_f
    public :: OpenRSPSetLinearRSPSolver_f
#if defined(OPENRSP_PERTURBATION_FREE)
    public :: OpenRSPSetPerturbations_f
#endif
    public :: OpenRSPSetPDBS_f
    public :: OpenRSPAddOneOper_f
    public :: OpenRSPAddTwoOper_f
    public :: OpenRSPAddXCFun_f
    public :: OpenRSPSetNucContributions_f
    public :: OpenRSPAssemble_f
    public :: OpenRSPWrite_f
    public :: OpenRSPGetRSPFun_f
    !public :: OpenRSPGetResidue_f
    public :: OpenRSPDestroy_f

    interface 
        integer(C_INT) function f_api_OpenRSPCreate(open_rsp) &
            bind(C, name="f_api_OpenRSPCreate")
            use, intrinsic :: iso_c_binding
            type(C_PTR), intent(inout) :: open_rsp
        end function f_api_OpenRSPCreate
        integer(C_INT) function f_api_OpenRSPSetElecEOM(open_rsp,      &
                                                        elec_EOM_type) &
            bind(C, name="f_api_OpenRSPSetElecEOM")
            use, intrinsic :: iso_c_binding
            type(C_PTR), intent(inout) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: elec_EOM_type
        end function f_api_OpenRSPSetElecEOM
        integer(C_INT) function OpenRSPSetLinearRSPSolver(open_rsp,         &
                                                          user_ctx,         &
                                                          get_linear_rsp_solution) &
            bind(C, name="OpenRSPSetLinearRSPSolver")
            use, intrinsic :: iso_c_binding
            type(C_PTR), value, intent(in) :: open_rsp
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_linear_rsp_solution
        end function OpenRSPSetLinearRSPSolver
#if defined(OPENRSP_PERTURBATION_FREE)
        integer(C_INT) function OpenRSPSetPerturbations(open_rsp,        &
                                                        num_pert,        &
                                                        pert_labels,     &
                                                        pert_max_orders, &
                                                        pert_sizes,      &
                                                        user_ctx,        &
                                                        get_pert_comp,   &
                                                        get_pert_rank)   &
            bind(C, name="OpenRSPSetPerturbations")
            use, intrinsic :: iso_c_binding
            type(C_PTR), value, intent(in) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: pert_labels(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_sizes(sum(pert_max_orders))
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_pert_comp
            type(C_FUNPTR), value, intent(in) :: get_pert_rank
        end function OpenRSPSetPerturbations
#endif
        integer(C_INT) function OpenRSPSetPDBS(open_rsp,        &
                                               num_pert,        &
                                               pert_labels,     &
                                               pert_max_orders, &
                                               user_ctx,        &
                                               get_overlap_mat, &
                                               get_overlap_exp) &
            bind(C, name="OpenRSPSetPDBS")
            use, intrinsic :: iso_c_binding
            type(C_PTR), value, intent(in) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: pert_labels(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert)
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_overlap_mat
            type(C_FUNPTR), value, intent(in) :: get_overlap_exp
        end function OpenRSPSetPDBS
        integer(C_INT) function OpenRSPAddOneOper(open_rsp,         &
                                                  num_pert,         &
                                                  pert_labels,      &
                                                  pert_max_orders,  &
                                                  user_ctx,         &
                                                  get_one_oper_mat, &
                                                  get_one_oper_exp) &
            bind(C, name="OpenRSPAddOneOper")
            use, intrinsic :: iso_c_binding  
            type(C_PTR), value, intent(in) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: pert_labels(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert)
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_one_oper_mat
            type(C_FUNPTR), value, intent(in) :: get_one_oper_exp
        end function OpenRSPAddOneOper
        integer(C_INT) function OpenRSPAddTwoOper(open_rsp,         &
                                                  num_pert,         &
                                                  pert_labels,      &
                                                  pert_max_orders,  &
                                                  user_ctx,         &
                                                  get_two_oper_mat, &
                                                  get_two_oper_exp) &
            bind(C, name="OpenRSPAddTwoOper")
            use, intrinsic :: iso_c_binding
            type(C_PTR), value, intent(in) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: pert_labels(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert)
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_two_oper_mat
            type(C_FUNPTR), value, intent(in) :: get_two_oper_exp
        end function OpenRSPAddTwoOper
        integer(C_INT) function OpenRSPAddXCFun(open_rsp,        &
                                                num_pert,        &
                                                pert_labels,     &
                                                pert_max_orders, &
                                                user_ctx,        &
                                                get_xc_fun_mat,  & 
                                                get_xc_fun_exp)  & 
            bind(C, name="OpenRSPAddXCFun")
            use, intrinsic :: iso_c_binding
            type(C_PTR), value, intent(in) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: pert_labels(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert)
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_xc_fun_mat
            type(C_FUNPTR), value, intent(in) :: get_xc_fun_exp
        end function OpenRSPAddXCFun
        integer(C_INT) function OpenRSPSetNucContributions(open_rsp,        &
                                                           num_pert,        &
                                                           pert_labels,     &
                                                           pert_max_orders, &
                                                           user_ctx,        &
                                                           get_nuc_contrib, &
                                                           num_atoms)       &
            bind(C, name="OpenRSPSetNucContributions")
            use, intrinsic :: iso_c_binding
            type(C_PTR), value, intent(in) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: pert_labels(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert)
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_nuc_contrib
            integer(kind=C_QINT), value, intent(in) :: num_atoms
        end function OpenRSPSetNucContributions
        integer(C_INT) function OpenRSPAssemble(open_rsp) &
            bind(C, name="OpenRSPAssemble")
            use, intrinsic :: iso_c_binding
            type(C_PTR), value, intent(in) :: open_rsp
        end function OpenRSPAssemble
        integer(C_INT) function OpenRSPWrite(open_rsp, file_name) &
            bind(C, name="OpenRSPWrite")
            use, intrinsic :: iso_c_binding
            type(C_PTR), value, intent(in) :: open_rsp
            character(C_CHAR), intent(in) :: file_name(*)
        end function OpenRSPWrite
        integer(C_INT) function OpenRSPGetRSPFun(open_rsp,         &
                                                 ref_ham,          &
                                                 ref_state,        &
                                                 ref_overlap,      &
                                                 num_props,        &
                                                 len_tuple,        &
                                                 pert_tuple,       &
                                                 num_freq_configs, &
                                                 pert_freqs,       &
                                                 kn_rules,         &
                                                 size_rsp_funs,    &
                                                 rsp_funs)         &
            bind(C, name="OpenRSPGetRSPFun")
            use, intrinsic :: iso_c_binding
            type(C_PTR), value, intent(in) :: open_rsp
            type(C_PTR), value, intent(in) :: ref_ham
            type(C_PTR), value, intent(in) :: ref_state
            type(C_PTR), value, intent(in) :: ref_overlap
            integer(kind=C_QINT), value, intent(in) :: num_props
            integer(kind=C_QINT), intent(in) :: len_tuple(num_props)
            integer(kind=C_QINT), intent(in) :: pert_tuple(sum(len_tuple))
            integer(kind=C_QINT), intent(in) :: num_freq_configs(num_props)
            real(kind=C_QREAL), intent(in) :: pert_freqs(2*dot_product(len_tuple,num_freq_configs))
            integer(kind=C_QINT), intent(in) :: kn_rules(num_props)
            integer(kind=C_QINT), value, intent(in) :: size_rsp_funs
            real(kind=C_QREAL), intent(out) :: rsp_funs(2*size_rsp_funs)
        end function OpenRSPGetRSPFun
        integer(C_INT) function f_api_OpenRSPDestroy(open_rsp) &
            bind(C, name="f_api_OpenRSPDestroy")
            use, intrinsic :: iso_c_binding
            type(C_PTR), intent(inout) :: open_rsp
        end function f_api_OpenRSPDestroy
    end interface

    contains

    function OpenRSPCreate_f(open_rsp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        ierr = f_api_OpenRSPCreate(open_rsp%c_rsp)
        nullify(open_rsp%solver_fun)
#if defined(OPENRSP_PERTURBATION_FREE)
        nullify(open_rsp%pert_fun)
#endif
        nullify(open_rsp%overlap_fun)
        nullify(open_rsp%list_one_oper)
        nullify(open_rsp%list_two_oper)
        nullify(open_rsp%list_xc_fun)
        nullify(open_rsp%nuc_contrib_fun)
    end function OpenRSPCreate_f

    function OpenRSPSetElecEOM_f(open_rsp, elec_EOM_type) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: elec_EOM_type
        ierr = f_api_OpenRSPSetElecEOM(open_rsp%c_rsp, elec_EOM_type)
    end function OpenRSPSetElecEOM_f

    function OpenRSPSetLinearRSPSolver_f(open_rsp,        &
#if defined(OPENRSP_F_USER_CONTEXT)
                                         user_ctx,        &
#endif
                                         get_linear_rsp_solution) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_linear_rsp_solution(size_pert,     &
                                               num_freq_sums, &
                                               freq_sums,     &
                                               RHS_mat,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                                               len_ctx,       &
                                               user_ctx,      &
#endif
                                               rsp_param)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                integer(kind=QINT), intent(in) :: size_pert
                integer(kind=QINT), intent(in) :: num_freq_sums
                real(kind=QREAL), intent(in) :: freq_sums(2*num_freq_sums)
                type(QcMat), intent(in) :: RHS_mat(size_pert*num_freq_sums)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                type(QcMat), intent(inout) :: rsp_param(size_pert*num_freq_sums)
            end subroutine get_linear_rsp_solution
            subroutine RSPSolverGetLinearRSPSolution_f(size_pert,     &
                                                       num_freq_sums, &
                                                       freq_sums,     &
                                                       RHS_mat,       &
                                                       user_ctx,      &
                                                       rsp_param)     &
                bind(C, name="RSPSolverGetLinearRSPSolution_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: size_pert
                integer(kind=C_QINT), value, intent(in) :: num_freq_sums
                real(kind=C_QREAL), intent(in) :: freq_sums(2*num_freq_sums)
                type(C_PTR), intent(in) :: RHS_mat(size_pert*num_freq_sums)
                type(C_PTR), value, intent(in) :: user_ctx
                type(C_PTR), intent(inout) :: rsp_param(size_pert*num_freq_sums)
            end subroutine RSPSolverGetLinearRSPSolution_f
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
                               get_linear_rsp_solution)
        ierr = OpenRSPSetLinearRSPSolver(open_rsp%c_rsp,             &
                                         c_loc(open_rsp%solver_fun), &
                                         c_funloc(RSPSolverGetLinearRSPSolution_f))
    end function OpenRSPSetLinearRSPSolver_f

#if defined(OPENRSP_PERTURBATION_FREE)
    function OpenRSPSetPerturbations_f(open_rsp,        &
                                       num_pert,        &
                                       pert_labels,     &
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
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert)
        integer(kind=QINT), intent(in) :: pert_sizes(:)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_pert_comp(pert_label,      &
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
            end subroutine get_pert_comp
            subroutine get_pert_rank(pert_label,       &
                                     pert_num_comp,    &
                                     pert_components,  &
                                     pert_comp_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                     len_ctx,          &
                                     user_ctx,         &
#endif
                                     pert_rank)
                use qcmatrix_f, only: QINT
                integer(kind=QINT), intent(in) :: pert_label
                integer(kind=QINT), intent(in) :: pert_num_comp
                integer(kind=QINT), intent(in) :: pert_components(pert_num_comp)
                integer(kind=QINT), intent(in) :: pert_comp_orders(pert_num_comp)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(out) :: pert_rank
            end subroutine get_pert_rank
            subroutine RSPPertGetComp_f(pert_label,       &
                                        pert_order,       &
                                        pert_rank,        &
                                        user_ctx,         &
                                        pert_num_comp,    &
                                        pert_components,  &
                                        pert_comp_orders) &
                bind(C, name="RSPPertGetComp_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: pert_label
                integer(kind=C_QINT), value, intent(in) :: pert_order
                integer(kind=C_QINT), value, intent(in) :: pert_rank
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), intent(out) :: pert_num_comp
                integer(kind=C_QINT), intent(out) :: pert_components(pert_order)
                integer(kind=C_QINT), intent(out) :: pert_comp_orders(pert_order)
            end subroutine RSPPertGetComp_f
            subroutine RSPPertGetRank_f(pert_label,       &
                                        pert_num_comp,    &
                                        pert_components,  &
                                        pert_comp_orders, &
                                        user_ctx,         &
                                        pert_rank)        &
                bind(C, name="RSPPertGetRank_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: pert_label
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
        ierr = OpenRSPSetPerturbations(open_rsp%c_rsp,             &
                                       num_pert,                   &
                                       pert_labels,                &
                                       pert_max_orders,            &
                                       pert_sizes,                 &
                                       c_loc(open_rsp%pert_fun),   &
                                       c_funloc(RSPPertGetComp_f), &
                                       c_funloc(RSPPertGetRank_f))
    end function OpenRSPSetPerturbations_f
#endif

    function OpenRSPSetPDBS_f(open_rsp,        &
                              num_pert,        &
                              pert_labels,     &
                              pert_max_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                              user_ctx,        &
#endif
                              get_overlap_mat, &
                              get_overlap_exp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_overlap_mat(bra_len_tuple,  &
                                       bra_pert_tuple, &
                                       ket_len_tuple,  &
                                       ket_pert_tuple, &
                                       len_tuple,      &
                                       pert_tuple,     &
#if defined(OPENRSP_F_USER_CONTEXT)
                                       len_ctx,        &
                                       user_ctx,       &
#endif
                                       num_int,        &
                                       val_int)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                integer(kind=QINT), intent(in) :: bra_len_tuple
                integer(kind=QINT), intent(in) :: bra_pert_tuple(bra_len_tuple)
                integer(kind=QINT), intent(in) :: ket_len_tuple
                integer(kind=QINT), intent(in) :: ket_pert_tuple(ket_len_tuple)
                integer(kind=QINT), intent(in) :: len_tuple
                integer(kind=QINT), intent(in) :: pert_tuple(len_tuple)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_int
                type(QcMat), intent(inout) :: val_int(num_int)
            end subroutine get_overlap_mat
            subroutine get_overlap_exp(bra_len_tuple,  &
                                       bra_pert_tuple, &
                                       ket_len_tuple,  &
                                       ket_pert_tuple, &
                                       len_tuple,      &
                                       pert_tuple,     &
                                       num_dmat,       &
                                       dens_mat,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                                       len_ctx,        &
                                       user_ctx,       &
#endif
                                       num_exp,        &
                                       val_exp)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                integer(kind=QINT), intent(in) :: bra_len_tuple
                integer(kind=QINT), intent(in) :: bra_pert_tuple(bra_len_tuple)
                integer(kind=QINT), intent(in) :: ket_len_tuple
                integer(kind=QINT), intent(in) :: ket_pert_tuple(ket_len_tuple)
                integer(kind=QINT), intent(in) :: len_tuple
                integer(kind=QINT), intent(in) :: pert_tuple(len_tuple)
                integer(kind=QINT), intent(in) :: num_dmat
                type(QcMat), intent(in) :: dens_mat(num_dmat)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_exp
                real(kind=QREAL), intent(inout) :: val_exp(num_exp)
            end subroutine get_overlap_exp
            subroutine RSPOverlapGetMat_f(bra_len_tuple,  &
                                          bra_pert_tuple, &
                                          ket_len_tuple,  &
                                          ket_pert_tuple, &
                                          len_tuple,      &
                                          pert_tuple,     &
                                          user_ctx,       &
                                          num_int,        &
                                          val_int)        &
                bind(C, name="RSPOverlapGetMat_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: bra_len_tuple
                integer(kind=C_QINT), intent(in) :: bra_pert_tuple(bra_len_tuple)
                integer(kind=C_QINT), value, intent(in) :: ket_len_tuple
                integer(kind=C_QINT), intent(in) :: ket_pert_tuple(ket_len_tuple)
                integer(kind=C_QINT), value, intent(in) :: len_tuple
                integer(kind=C_QINT), intent(in) :: pert_tuple(len_tuple)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_int
                type(C_PTR), intent(inout) :: val_int(num_int)
            end subroutine RSPOverlapGetMat_f
            subroutine RSPOverlapGetExp_f(bra_len_tuple,  &
                                          bra_pert_tuple, &
                                          ket_len_tuple,  &
                                          ket_pert_tuple, &
                                          len_tuple,      &
                                          pert_tuple,     &
                                          num_dmat,       &
                                          dens_mat,       &
                                          user_ctx,       &
                                          num_exp,        &
                                          val_exp)        &
                bind(C, name="RSPOverlapGetExp_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: bra_len_tuple
                integer(kind=C_QINT), intent(in) :: bra_pert_tuple(bra_len_tuple)
                integer(kind=C_QINT), value, intent(in) :: ket_len_tuple
                integer(kind=C_QINT), intent(in) :: ket_pert_tuple(ket_len_tuple)
                integer(kind=C_QINT), value, intent(in) :: len_tuple
                integer(kind=C_QINT), intent(in) :: pert_tuple(len_tuple)
                integer(kind=C_QINT), value, intent(in) :: num_dmat
                type(C_PTR), intent(in) :: dens_mat(num_dmat)
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
        ierr = OpenRSPSetPDBS(open_rsp%c_rsp,               &
                              num_pert,                     &
                              pert_labels,                  &
                              pert_max_orders,              &
                              c_loc(open_rsp%overlap_fun),  &
                              c_funloc(RSPOverlapGetMat_f), &
                              c_funloc(RSPOverlapGetExp_f))
    end function OpenRSPSetPDBS_f

    function OpenRSPAddOneOper_f(open_rsp,         &
                                 num_pert,         &
                                 pert_labels,      &
                                 pert_max_orders,  &
#if defined(OPENRSP_F_USER_CONTEXT)
                                 user_ctx,         &
#endif
                                 get_one_oper_mat, &
                                 get_one_oper_exp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_one_oper_mat(len_tuple,  &
                                        pert_tuple, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                        len_ctx,    &
                                        user_ctx,   &
#endif
                                        num_int,    &
                                        val_int)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                integer(kind=QINT), intent(in) :: len_tuple
                integer(kind=QINT), intent(in) :: pert_tuple(len_tuple)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_int
                type(QcMat), intent(inout) :: val_int(num_int)
            end subroutine get_one_oper_mat
            subroutine get_one_oper_exp(len_tuple,  &
                                        pert_tuple, &
                                        num_dmat,   &
                                        dens_mat,   &
#if defined(OPENRSP_F_USER_CONTEXT)
                                        len_ctx,    &
                                        user_ctx,   &
#endif
                                        num_exp,    &
                                        val_exp)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                integer(kind=QINT), intent(in) :: len_tuple
                integer(kind=QINT), intent(in) :: pert_tuple(len_tuple)
                integer(kind=QINT), intent(in) :: num_dmat
                type(QcMat), intent(in) :: dens_mat(num_dmat)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_exp
                real(kind=QREAL), intent(inout) :: val_exp(num_exp)
            end subroutine get_one_oper_exp
            subroutine RSPOneOperGetMat_f(len_tuple,  &
                                          pert_tuple, &
                                          user_ctx,   &
                                          num_int,    &
                                          val_int)    &
                bind(C, name="RSPOneOperGetMat_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: len_tuple
                integer(kind=C_QINT), intent(in) :: pert_tuple(len_tuple)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_int
                type(C_PTR), intent(inout) :: val_int(num_int)
            end subroutine RSPOneOperGetMat_f
            subroutine RSPOneOperGetExp_f(len_tuple,  &
                                          pert_tuple, &
                                          num_dmat,   &
                                          dens_mat,   &
                                          user_ctx,   &
                                          num_exp,    &
                                          val_exp)    &
                bind(C, name="RSPOneOperGetExp_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: len_tuple
                integer(kind=C_QINT), intent(in) :: pert_tuple(len_tuple)
                integer(kind=C_QINT), value, intent(in) :: num_dmat
                type(C_PTR), intent(in) :: dens_mat(num_dmat)
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
        ierr = OpenRSPAddOneOper(open_rsp%c_rsp,                   &
                                 num_pert,                         &
                                 pert_labels,                      &
                                 pert_max_orders,                  &
                                 c_loc(cur_one_oper%one_oper_fun), &
                                 c_funloc(RSPOneOperGetMat_f),     &
                                 c_funloc(RSPOneOperGetExp_f))
    end function OpenRSPAddOneOper_f

    function OpenRSPAddTwoOper_f(open_rsp,         &
                                 num_pert,         &
                                 pert_labels,      &
                                 pert_max_orders,  &
#if defined(OPENRSP_F_USER_CONTEXT)
                                 user_ctx,         &
#endif
                                 get_two_oper_mat, &
                                 get_two_oper_exp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_two_oper_mat(len_tuple,  &
                                        pert_tuple, &
                                        num_dmat,   &
                                        dens_mat,   &
#if defined(OPENRSP_F_USER_CONTEXT)
                                        len_ctx,    &
                                        user_ctx,   &
#endif
                                        num_int,    &
                                        val_int)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                integer(kind=QINT), intent(in) :: len_tuple
                integer(kind=QINT), intent(in) :: pert_tuple(len_tuple)
                integer(kind=QINT), intent(in) :: num_dmat
                type(QcMat), intent(in) :: dens_mat(num_dmat)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_int
                type(QcMat), intent(inout) :: val_int(num_int)
            end subroutine get_two_oper_mat
            subroutine get_two_oper_exp(len_tuple,      &
                                        pert_tuple,     &
                                        len_dmat_tuple, &
                                        num_LHS_dmat,   &
                                        LHS_dens_mat,   &
                                        num_RHS_dmat,   &
                                        RHS_dens_mat,   &
#if defined(OPENRSP_F_USER_CONTEXT)
                                        len_ctx,        &
                                        user_ctx,       &
#endif
                                        num_exp,        &
                                        val_exp)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                integer(kind=QINT), intent(in) :: len_tuple
                integer(kind=QINT), intent(in) :: pert_tuple(len_tuple)
                integer(kind=QINT), intent(in) :: len_dmat_tuple
                integer(kind=QINT), intent(in) :: num_LHS_dmat(len_dmat_tuple)
                type(QcMat), intent(in) :: LHS_dens_mat(sum(num_LHS_dmat))
                integer(kind=QINT), intent(in) :: num_RHS_dmat(len_dmat_tuple)
                type(QcMat), intent(in) :: RHS_dens_mat(sum(num_RHS_dmat))
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_exp
                real(kind=QREAL), intent(inout) :: val_exp(num_exp)
            end subroutine get_two_oper_exp
            subroutine RSPTwoOperGetMat_f(len_tuple,  &
                                          pert_tuple, &
                                          num_dmat,   &
                                          dens_mat,   &
                                          user_ctx,   &
                                          num_int,    &
                                          val_int)    &
                bind(C, name="RSPTwoOperGetMat_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: len_tuple
                integer(kind=C_QINT), intent(in) :: pert_tuple(len_tuple)
                integer(kind=C_QINT), value, intent(in) :: num_dmat
                type(C_PTR), intent(in) :: dens_mat(num_dmat)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_int
                type(C_PTR), intent(inout) :: val_int(num_int)
            end subroutine RSPTwoOperGetMat_f
            subroutine RSPTwoOperGetExp_f(len_tuple,      &
                                          pert_tuple,     &
                                          len_dmat_tuple, &
                                          num_LHS_dmat,   &
                                          LHS_dens_mat,   &
                                          num_RHS_dmat,   &
                                          RHS_dens_mat,   &
                                          user_ctx,       &
                                          num_exp,        &
                                          val_exp)        &
                bind(C, name="RSPTwoOperGetExp_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: len_tuple
                integer(kind=C_QINT), intent(in) :: pert_tuple(len_tuple)
                integer(kind=C_QINT), value, intent(in) :: len_dmat_tuple
                integer(kind=C_QINT), intent(in) :: num_LHS_dmat(len_dmat_tuple)
                type(C_PTR), intent(in) :: LHS_dens_mat(sum(num_LHS_dmat))
                integer(kind=C_QINT), intent(in) :: num_RHS_dmat(len_dmat_tuple)
                type(C_PTR), intent(in) :: RHS_dens_mat(sum(num_RHS_dmat))
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
        ierr = OpenRSPAddTwoOper(open_rsp%c_rsp,                   &
                                 num_pert,                         &
                                 pert_labels,                      &
                                 pert_max_orders,                  &
                                 c_loc(cur_two_oper%two_oper_fun), &
                                 c_funloc(RSPTwoOperGetMat_f),     &
                                 c_funloc(RSPTwoOperGetExp_f))
    end function OpenRSPAddTwoOper_f

    function OpenRSPAddXCFun_f(open_rsp,        &
                               num_pert,        &
                               pert_labels,     &
                               pert_max_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                               user_ctx,        &
#endif
                               get_xc_fun_mat,  &
                               get_xc_fun_exp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_xc_fun_mat(len_tuple,        &
                                      pert_tuple,       &
                                      num_freq_configs, &
                                      len_dmat_tuple,   &
                                      idx_dmat_tuple,   &
                                      num_dmat,         &
                                      dens_mat,         &
#if defined(OPENRSP_F_USER_CONTEXT)
                                      len_ctx,          &
                                      user_ctx,         &
#endif
                                      num_int,          &
                                      val_int)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                integer(kind=QINT), intent(in) :: len_tuple
                integer(kind=QINT), intent(in) :: pert_tuple(len_tuple)
                integer(kind=QINT), intent(in) :: num_freq_configs
                integer(kind=QINT), intent(in) :: len_dmat_tuple
                integer(kind=QINT), intent(in) :: idx_dmat_tuple(len_dmat_tuple)
                integer(kind=QINT), intent(in) :: num_dmat
                type(QcMat), intent(in) :: dens_mat(num_dmat)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_int
                type(QcMat), intent(inout) :: val_int(num_int)
            end subroutine get_xc_fun_mat
            subroutine get_xc_fun_exp(len_tuple,        &
                                      pert_tuple,       &
                                      num_freq_configs, &
                                      len_dmat_tuple,   &
                                      idx_dmat_tuple,   &
                                      num_dmat,         &
                                      dens_mat,         &
#if defined(OPENRSP_F_USER_CONTEXT)
                                      len_ctx,          &
                                      user_ctx,         &
#endif
                                      num_exp,          &
                                      val_exp)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                integer(kind=QINT), intent(in) :: len_tuple
                integer(kind=QINT), intent(in) :: pert_tuple(len_tuple)
                integer(kind=QINT), intent(in) :: num_freq_configs
                integer(kind=QINT), intent(in) :: len_dmat_tuple
                integer(kind=QINT), intent(in) :: idx_dmat_tuple(len_dmat_tuple)
                integer(kind=QINT), intent(in) :: num_dmat
                type(QcMat), intent(in) :: dens_mat(num_dmat)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_exp
                real(kind=QREAL), intent(inout) :: val_exp(num_exp)
            end subroutine get_xc_fun_exp
            subroutine RSPXCFunGetMat_f(len_tuple,        &
                                        pert_tuple,       &
                                        num_freq_configs, &
                                        len_dmat_tuple,   &
                                        idx_dmat_tuple,   &
                                        num_dmat,         &
                                        dens_mat,         &
                                        user_ctx,         &
                                        num_int,          &
                                        val_int)          &
                bind(C, name="RSPXCFunGetMat_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: len_tuple
                integer(kind=C_QINT), intent(in) :: pert_tuple(len_tuple)
                integer(kind=C_QINT), value, intent(in) :: num_freq_configs
                integer(kind=C_QINT), value, intent(in) :: len_dmat_tuple
                integer(kind=C_QINT), intent(in) :: idx_dmat_tuple(len_dmat_tuple)
                integer(kind=C_QINT), value, intent(in) :: num_dmat
                type(C_PTR), intent(in) :: dens_mat(num_dmat)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_int
                type(C_PTR), intent(inout) :: val_int(num_int)
            end subroutine RSPXCFunGetMat_f
            subroutine RSPXCFunGetExp_f(len_tuple,        &
                                        pert_tuple,       &
                                        num_freq_configs, &
                                        len_dmat_tuple,   &
                                        idx_dmat_tuple,   &
                                        num_dmat,         &
                                        dens_mat,         &
                                        user_ctx,         &
                                        num_exp,          &
                                        val_exp)          &
                bind(C, name="RSPXCFunGetExp_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: len_tuple
                integer(kind=C_QINT), intent(in) :: pert_tuple(len_tuple)
                integer(kind=C_QINT), value, intent(in) :: num_freq_configs
                integer(kind=C_QINT), value, intent(in) :: len_dmat_tuple
                integer(kind=C_QINT), intent(in) :: idx_dmat_tuple(len_dmat_tuple)
                integer(kind=C_QINT), value, intent(in) :: num_dmat
                type(C_PTR), intent(in) :: dens_mat(num_dmat)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_exp
                real(kind=C_QREAL), intent(inout) :: val_exp(num_exp)
            end subroutine RSPXCFunGetExp_f
        end interface
        type(XCFunList_f), pointer :: cur_xc_fun  !current XC functional
        ! inserts the context of callback functions to the tail of the linked list
        if (associated(open_rsp%list_xc_fun)) then
            cur_xc_fun => open_rsp%list_xc_fun
            do while (associated(cur_xc_fun%next_xc_fun))
                cur_xc_fun => cur_xc_fun%next_xc_fun
            end do
            allocate(cur_xc_fun%next_xc_fun)
            cur_xc_fun => cur_xc_fun%next_xc_fun
        else
            allocate(open_rsp%list_xc_fun)
            cur_xc_fun => open_rsp%list_xc_fun
        end if
        allocate(cur_xc_fun%xc_fun_fun)
        nullify(cur_xc_fun%next_xc_fun)
        ! adds context of callback functions of the new XC functional
        call RSPXCFunCreate_f(cur_xc_fun%xc_fun_fun, &
#if defined(OPENRSP_F_USER_CONTEXT)
                              user_ctx,              &
#endif
                              get_xc_fun_mat,        &
                              get_xc_fun_exp)
        ierr = OpenRSPAddXCFun(open_rsp%c_rsp,               &
                               num_pert,                     &
                               pert_labels,                  &
                               pert_max_orders,              &
                               c_loc(cur_xc_fun%xc_fun_fun), &
                               c_funloc(RSPXCFunGetMat_f),   &
                               c_funloc(RSPXCFunGetExp_f))
    end function OpenRSPAddXCFun_f

    function OpenRSPSetNucContributions_f(open_rsp,        &
                                          num_pert,        &
                                          pert_labels,     &
                                          pert_max_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                          user_ctx,        &
#endif
                                          get_nuc_contrib, &
                                          num_atoms) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        integer(kind=QINT), intent(in) :: num_atoms
        interface
            subroutine get_nuc_contrib(len_tuple,  &
                                       pert_tuple, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                       len_ctx,    &
                                       user_ctx,   &
#endif
                                       size_pert,  &
                                       val_nuc)
                use qcmatrix_f, only: QINT,QREAL
                integer(kind=QINT), intent(in) :: len_tuple
                integer(kind=QINT), intent(in) :: pert_tuple(len_tuple)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: size_pert
                real(kind=QREAL), intent(inout) :: val_nuc(size_pert)
            end subroutine get_nuc_contrib
            subroutine RSPNucHamiltonGetContributions_f(len_tuple,  &
                                                        pert_tuple, &
                                                        user_ctx,   &
                                                        size_pert,  &
                                                        val_nuc)    &
                bind(C, name="RSPNucHamiltonGetContributions_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: len_tuple
                integer(kind=C_QINT), intent(in) :: pert_tuple(len_tuple)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: size_pert
                real(kind=C_QREAL), intent(inout) :: val_nuc(size_pert)
            end subroutine RSPNucHamiltonGetContributions_f
        end interface
        if (associated(open_rsp%nuc_contrib_fun)) then
            call RSPNucHamiltonDestroy_f(open_rsp%nuc_contrib_fun)
        else
            allocate(open_rsp%nuc_contrib_fun)
        end if
        ! adds context of callback function of the nuclear Hamiltonian
        call RSPNucHamiltonCreate_f(open_rsp%nuc_contrib_fun, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                    user_ctx,                 &
#endif
                                    get_nuc_contrib)
        ierr = OpenRSPSetNucContributions(open_rsp%c_rsp,                             &
                                          num_pert,                                   &
                                          pert_labels,                                &
                                          pert_max_orders,                            &
                                          c_loc(open_rsp%nuc_contrib_fun),            &
                                          c_funloc(RSPNucHamiltonGetContributions_f), &
                                          num_atoms)
    end function OpenRSPSetNucContributions_f

    function OpenRSPAssemble_f(open_rsp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        ierr = OpenRSPAssemble(open_rsp%c_rsp)
    end function OpenRSPAssemble_f

    function OpenRSPWrite_f(open_rsp, file_name) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(in) :: open_rsp
        character*(*), intent(in) :: file_name
        ierr = OpenRSPWrite(open_rsp%c_rsp, file_name//C_NULL_CHAR)
    end function OpenRSPWrite_f

    function OpenRSPGetRSPFun_f(open_rsp,         &
                                ref_ham,          &
                                ref_state,        &
                                ref_overlap,      &
                                num_props,        &
                                len_tuple,        &
                                pert_tuple,       &
                                num_freq_configs, &
                                pert_freqs,       &
                                kn_rules,         &
                                size_rsp_funs,    &
                                rsp_funs) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(in) :: open_rsp
        type(QcMat), target, intent(in) :: ref_ham
        type(QcMat), target, intent(in) :: ref_state
        type(QcMat), target, intent(in) :: ref_overlap
        integer(kind=QINT), intent(in) :: num_props
        integer(kind=QINT), intent(in) :: len_tuple(num_props)
        integer(kind=QINT), intent(in) :: pert_tuple(sum(len_tuple))
        integer(kind=QINT), intent(in) :: num_freq_configs(num_props)
        real(kind=QREAL), intent(in) :: pert_freqs(2*dot_product(len_tuple,num_freq_configs))
        integer(kind=QINT), intent(in) :: kn_rules(num_props)
        integer(kind=QINT), intent(in) :: size_rsp_funs
        real(kind=QREAL), intent(out) :: rsp_funs(2*size_rsp_funs)
        type(C_PTR) c_ref_ham(1)
        type(C_PTR) c_ref_state(1)
        type(C_PTR) c_ref_overlap(1)
        ierr = QcMat_C_LOC((/ref_ham/), c_ref_ham)
        if (ierr==QFAILURE) return
        ierr = QcMat_C_LOC((/ref_state/), c_ref_state)
        if (ierr==QFAILURE) return
        ierr = QcMat_C_LOC((/ref_overlap/), c_ref_overlap)
        if (ierr==QFAILURE) return
        ierr = OpenRSPGetRSPFun(open_rsp%c_rsp,    &
                                c_ref_ham(1),      &
                                c_ref_state(1),    &
                                c_ref_overlap(1),  &
                                num_props,         &
                                len_tuple,         &
                                pert_tuple,        &
                                num_freq_configs,  &
                                pert_freqs,        &
                                kn_rules,          &
                                size_rsp_funs,     &
                                rsp_funs)
        c_ref_ham(1) = C_NULL_PTR
        c_ref_state(1) = C_NULL_PTR
        c_ref_overlap(1) = C_NULL_PTR
    end function OpenRSPGetRSPFun_f

    function OpenRSPDestroy_f(open_rsp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        type(OneOperList_f), pointer :: cur_one_oper   !current one-electron operator
        type(OneOperList_f), pointer :: next_one_oper  !next one-electron operator
        type(TwoOperList_f), pointer :: cur_two_oper   !current two-electron operator
        type(TwoOperList_f), pointer :: next_two_oper  !next two-electron operator
        type(XCFunList_f), pointer :: cur_xc_fun       !current XC functional
        type(XCFunList_f), pointer :: next_xc_fun      !next XC functional
        ierr = f_api_OpenRSPDestroy(open_rsp%c_rsp)
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
        ! cleans up the linked list of context of callback subroutines of XC functionals
        cur_xc_fun => open_rsp%list_xc_fun
        do while (associated(cur_xc_fun))
            next_xc_fun => cur_xc_fun%next_xc_fun
            if (associated(cur_xc_fun%xc_fun_fun)) then
                call RSPXCFunDestroy_f(cur_xc_fun%xc_fun_fun)
                deallocate(cur_xc_fun%xc_fun_fun)
                nullify(cur_xc_fun%xc_fun_fun)
            end if
            deallocate(cur_xc_fun)
            nullify(cur_xc_fun)
            cur_xc_fun => next_xc_fun
        end do
        ! cleans up callback subroutine of nuclear Hamiltonian
        if (associated(open_rsp%nuc_contrib_fun)) then
            call RSPNucHamiltonDestroy_f(open_rsp%nuc_contrib_fun)
            deallocate(open_rsp%nuc_contrib_fun)
            nullify(open_rsp%nuc_contrib_fun)
        end if
    end function OpenRSPDestroy_f

end module openrsp_f
