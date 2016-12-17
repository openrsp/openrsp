!!
!! <QCLANG='Fortran'>
!! <para>
!!   Fortran users should use the module <OpenRSP_f> in their codes to access
!!   the functionalities of OpenRSP. We have used the same name for Fortran
!!   data types and constants, for instance <OpenRSP>; macro definitions are
!!   also controlled by the same names, such as <OPENRSP_USER_CONTEXT>; however
!!   all Fortran modules and functions are appended by <c>_f</c>.
!! </para>
!!
!! OpenRSP: open-ended library for response theory
!! Copyright 2015 Radovan Bast,
!!                Daniel H. Friese,
!!                Bin Gao,
!!                Dan J. Jonsson,
!!                Magnus Ringholm,
!!                Kenneth Ruud,
!!                Andreas Thorvaldsen
!!
!! OpenRSP is free software: you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as
!! published by the Free Software Foundation, either version 3 of
!! the License, or (at your option) any later version.
!!
!! OpenRSP is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!! GNU Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public
!! License along with OpenRSP. If not, see <http://www.gnu.org/licenses/>.

!! <module name='OpenRSP_f' author='Bin Gao' date='2014-07-12'>
!!   The module file of OpenRSP library for Fortran users
!! </module>

! basic data types
#include "api/qcmatrix_c_type.h"

module OpenRSP_f

    use, intrinsic :: iso_c_binding
    use qcmatrix_f, only: QINT,     &
                          QREAL,    &
                          QFAILURE, &
                          QcMat,    &
                          QcMat_C_LOC
    use RSPSolver_f, only: SolverFun_f,       &
                           RSPSolverCreate_f, &
                           RSPSolverDestroy_f
    use RSPPerturbation_f, only: QcPertInt,       &
                                 PertFun_f,       &
                                 RSPPertCreate_f, &
                                 RSPPertDestroy_f
    use RSPOverlap_f, only: OverlapFun_f,       &
                            RSPOverlapCreate_f, &
                            RSPOverlapDestroy_f
    use RSPOneOper_f, only: OneOperFun_f,       &
                            RSPOneOperCreate_f, &
                            RSPOneOperDestroy_f
    use RSPTwoOper_f, only: TwoOperFun_f,       &
                            RSPTwoOperCreate_f, &
                            RSPTwoOperDestroy_f
    use RSPXCFun_f, only: XCFunFun_f,       &
                          RSPXCFunCreate_f, &
                          RSPXCFunDestroy_f
    use RSPNucHamilton_f, only: NucHamiltonFun_f,       &
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
        type(XCFunFun_f), pointer :: xcfun_fun => null()
        type(XCFunList_f), pointer :: next_xc_fun => null()
    end type XCFunList_f

    ! OpenRSP type (inspired by http://wiki.rac.manchester.ac.uk/community/GPU/GpuFaq/FortranGPU)
    type, public :: OpenRSP
        private
        type(C_PTR) :: c_rsp = C_NULL_PTR
        type(SolverFun_f), pointer :: solver_fun => null()
        type(PertFun_f), pointer :: pert_fun => null()
        type(OverlapFun_f), pointer :: overlap_fun => null()
        type(OneOperList_f), pointer :: list_one_oper => null()
        type(TwoOperList_f), pointer :: list_two_oper => null()
        type(XCFunList_f), pointer :: list_xc_fun => null()
        type(NucHamiltonFun_f), pointer :: nuc_hamilton_fun => null()
    end type OpenRSP

    ! functions provided by the Fortran APIs
    public :: OpenRSPCreate_f
    !public :: OpenRSPSetElecEOM_f
    public :: OpenRSPSetLinearRSPSolver_f
    public :: OpenRSPSetPerturbations_f
    public :: OpenRSPSetOverlap_f
    public :: OpenRSPAddOneOper_f
    public :: OpenRSPAddTwoOper_f
    public :: OpenRSPAddXCFun_f
    public :: OpenRSPSetNucHamilton_f
    public :: OpenRSPAssemble_f
    public :: OpenRSPWrite_f
    public :: OpenRSPGetRSPFun_f
    public :: OpenRSPGetResidue_f
    public :: OpenRSPDestroy_f

    interface OpenRSPWrite_f
        module procedure OpenRSPWritebyFileName_f
        module procedure OpenRSPWritebyUnit_f
    end interface OpenRSPWrite_f

    interface 
        integer(C_INT) function OpenRSPCreateFortranAdapter(open_rsp) &
            bind(C, name="OpenRSPCreateFortranAdapter")
            use, intrinsic :: iso_c_binding
            type(C_PTR), intent(inout) :: open_rsp
        end function OpenRSPCreateFortranAdapter
        !integer(C_INT) function f_api_OpenRSPSetElecEOM(open_rsp,      &
        !                                                elec_EOM_type) &
        !    bind(C, name="f_api_OpenRSPSetElecEOM")
        !    use, intrinsic :: iso_c_binding
        !    type(C_PTR), intent(inout) :: open_rsp
        !    integer(kind=C_QINT), value, intent(in) :: elec_EOM_type
        !end function f_api_OpenRSPSetElecEOM
        integer(C_INT) function OpenRSPSetLinearRSPSolver(open_rsp,         &
                                                          user_ctx,         &
                                                          get_linear_rsp_solution) &
            bind(C, name="OpenRSPSetLinearRSPSolver")
            use, intrinsic :: iso_c_binding
            type(C_PTR), value, intent(in) :: open_rsp
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_linear_rsp_solution
        end function OpenRSPSetLinearRSPSolver
        integer(C_INT) function OpenRSPSetPerturbations(open_rsp,               &
                                                        num_pert_lab,           &
                                                        pert_labels,            &
                                                        pert_max_orders,        &
                                                        pert_num_comps,         &
                                                        user_ctx,               &
                                                        get_pert_concatenation) &
            bind(C, name="OpenRSPSetPerturbations")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            type(C_PTR), value, intent(in) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert_lab
            integer(kind=C_QCPERTINT), intent(in) :: pert_labels(num_pert_lab)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert_lab)
            integer(kind=C_QINT), intent(in) :: pert_num_comps(sum(pert_max_orders))
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_pert_concatenation
        end function OpenRSPSetPerturbations
        integer(C_INT) function OpenRSPSetOverlap(open_rsp,        &
                                                  num_pert_lab,    &
                                                  pert_labels,     &
                                                  pert_max_orders, &
                                                  user_ctx,        &
                                                  get_overlap_mat, &
                                                  get_overlap_exp) &
            bind(C, name="OpenRSPSetOverlap")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            type(C_PTR), value, intent(in) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert_lab
            integer(kind=C_QCPERTINT), intent(in) :: pert_labels(num_pert_lab)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert_lab)
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_overlap_mat
            type(C_FUNPTR), value, intent(in) :: get_overlap_exp
        end function OpenRSPSetOverlap
        integer(C_INT) function OpenRSPAddOneOper(open_rsp,         &
                                                  num_pert_lab,     &
                                                  pert_labels,      &
                                                  pert_max_orders,  &
                                                  user_ctx,         &
                                                  get_one_oper_mat, &
                                                  get_one_oper_exp) &
            bind(C, name="OpenRSPAddOneOper")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            type(C_PTR), value, intent(in) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert_lab
            integer(kind=C_QCPERTINT), intent(in) :: pert_labels(num_pert_lab)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert_lab)
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_one_oper_mat
            type(C_FUNPTR), value, intent(in) :: get_one_oper_exp
        end function OpenRSPAddOneOper
        integer(C_INT) function OpenRSPAddTwoOper(open_rsp,         &
                                                  num_pert_lab,     &
                                                  pert_labels,      &
                                                  pert_max_orders,  &
                                                  user_ctx,         &
                                                  get_two_oper_mat, &
                                                  get_two_oper_exp) &
            bind(C, name="OpenRSPAddTwoOper")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            type(C_PTR), value, intent(in) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert_lab
            integer(kind=C_QCPERTINT), intent(in) :: pert_labels(num_pert_lab)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert_lab)
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_two_oper_mat
            type(C_FUNPTR), value, intent(in) :: get_two_oper_exp
        end function OpenRSPAddTwoOper
        integer(C_INT) function OpenRSPAddXCFun(open_rsp,        &
                                                num_pert_lab,    &
                                                pert_labels,     &
                                                pert_max_orders, &
                                                user_ctx,        &
                                                get_xc_fun_mat,  & 
                                                get_xc_fun_exp)  & 
            bind(C, name="OpenRSPAddXCFun")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            type(C_PTR), value, intent(in) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert_lab
            integer(kind=C_QCPERTINT), intent(in) :: pert_labels(num_pert_lab)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert_lab)
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_xc_fun_mat
            type(C_FUNPTR), value, intent(in) :: get_xc_fun_exp
        end function OpenRSPAddXCFun
        integer(C_INT) function OpenRSPSetNucHamilton(open_rsp,        &
                                                      num_pert_lab,    &
                                                      pert_labels,     &
                                                      pert_max_orders, &
                                                      user_ctx,        &
                                                      get_nuc_contrib, &
                                                      num_atoms)       &
            bind(C, name="OpenRSPSetNucHamilton")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            type(C_PTR), value, intent(in) :: open_rsp
            integer(kind=C_QINT), value, intent(in) :: num_pert_lab
            integer(kind=C_QCPERTINT), intent(in) :: pert_labels(num_pert_lab)
            integer(kind=C_QINT), intent(in) :: pert_max_orders(num_pert_lab)
            type(C_PTR), value, intent(in) :: user_ctx
            type(C_FUNPTR), value, intent(in) :: get_nuc_contrib
            integer(kind=C_QINT), value, intent(in) :: num_atoms
        end function OpenRSPSetNucHamilton
        integer(C_INT) function OpenRSPAssemble(open_rsp) &
            bind(C, name="OpenRSPAssemble")
            use, intrinsic :: iso_c_binding
            type(C_PTR), value, intent(in) :: open_rsp
        end function OpenRSPAssemble
        integer(C_INT) function OpenRSPWriteFortranAdapter(open_rsp,  &
                                                           file_name) &
            bind(C, name="OpenRSPWriteFortranAdapter")
            use, intrinsic :: iso_c_binding
            type(C_PTR), value, intent(in) :: open_rsp
            character(C_CHAR), intent(in) :: file_name(*)
        end function OpenRSPWriteFortranAdapter
        integer(C_INT) function OpenRSPWriteStdOutFortranAdapter(open_rsp) &
            bind(C, name="OpenRSPWriteStdOutFortranAdapter")
            use, intrinsic :: iso_c_binding
            type(C_PTR), value, intent(in) :: open_rsp
        end function OpenRSPWriteStdOutFortranAdapter
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
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            type(C_PTR), value, intent(in) :: open_rsp
            type(C_PTR), value, intent(in) :: ref_ham
            type(C_PTR), value, intent(in) :: ref_state
            type(C_PTR), value, intent(in) :: ref_overlap
            integer(kind=C_QINT), value, intent(in) :: num_props
            integer(kind=C_QINT), intent(in) :: len_tuple(num_props)
            integer(kind=C_QCPERTINT), intent(in) :: pert_tuple(sum(len_tuple))
            integer(kind=C_QINT), intent(in) :: num_freq_configs(num_props)
            real(kind=C_QREAL), intent(in) :: &
                pert_freqs(2*dot_product(len_tuple,num_freq_configs))
            integer(kind=C_QINT), intent(in) :: kn_rules(num_props)
            integer(kind=C_QINT), value, intent(in) :: size_rsp_funs
            real(kind=C_QREAL), intent(out) :: rsp_funs(2*size_rsp_funs)
        end function OpenRSPGetRSPFun
        integer(C_INT) function OpenRSPGetResidue(open_rsp,         &
                                                  ref_ham,          &
                                                  ref_state,        &
                                                  ref_overlap,      &
                                                  order_residue,    &
                                                  num_excit,        &
                                                  excit_energy,     &
                                                  eigen_vector,     &
                                                  num_props,        &
                                                  len_tuple,        &
                                                  pert_tuple,       &
                                                  residue_num_pert, &
                                                  residue_idx_pert, &
                                                  num_freq_configs, &
                                                  pert_freqs,       &
                                                  kn_rules,         &
                                                  size_residues,    &
                                                  residues)         &
            bind(C, name="OpenRSPGetResidue")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            type(C_PTR), value, intent(in) :: open_rsp
            type(C_PTR), value, intent(in) :: ref_ham
            type(C_PTR), value, intent(in) :: ref_state
            type(C_PTR), value, intent(in) :: ref_overlap
            integer(kind=C_QINT), value, intent(in) :: order_residue
            integer(kind=C_QINT), value, intent(in) :: num_excit
            real(kind=C_QREAL), intent(in) :: excit_energy(order_residue*num_excit)
            type(C_PTR), intent(in) :: eigen_vector(order_residue*num_excit)
            integer(kind=C_QINT), value, intent(in) :: num_props
            integer(kind=C_QINT), intent(in) :: len_tuple(num_props)
            integer(kind=C_QCPERTINT), intent(in) :: pert_tuple(sum(len_tuple))
            integer(kind=C_QINT), intent(in) :: residue_num_pert(order_residue*num_props)
            integer(kind=C_QINT), intent(in) :: residue_idx_pert(sum(residue_num_pert))
            integer(kind=C_QINT), intent(in) :: num_freq_configs(num_props)
            real(kind=C_QREAL), intent(in) :: &
                pert_freqs(2*dot_product(len_tuple,num_freq_configs)*num_excit)
            integer(kind=C_QINT), intent(in) :: kn_rules(num_props)
            integer(kind=C_QINT), value, intent(in) :: size_residues
            real(kind=C_QREAL), intent(out) :: residues(2*size_residues)
        end function OpenRSPGetResidue
        integer(C_INT) function OpenRSPDestroyFortranAdapter(open_rsp) &
            bind(C, name="OpenRSPDestroyFortranAdapter")
            use, intrinsic :: iso_c_binding
            type(C_PTR), intent(inout) :: open_rsp
        end function OpenRSPDestroyFortranAdapter
    end interface

    contains

    function OpenRSPCreate_f(open_rsp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        ierr = OpenRSPCreateFortranAdapter(open_rsp%c_rsp)
        nullify(open_rsp%solver_fun)
        nullify(open_rsp%pert_fun)
        nullify(open_rsp%overlap_fun)
        nullify(open_rsp%list_one_oper)
        nullify(open_rsp%list_two_oper)
        nullify(open_rsp%list_xc_fun)
        nullify(open_rsp%nuc_hamilton_fun)
    end function OpenRSPCreate_f

    !function OpenRSPSetElecEOM_f(open_rsp, elec_EOM_type) result(ierr)
    !    integer(kind=4) :: ierr
    !    type(OpenRSP), intent(inout) :: open_rsp
    !    integer(kind=QINT), intent(in) :: elec_EOM_type
    !    ierr = f_api_OpenRSPSetElecEOM(open_rsp%c_rsp, elec_EOM_type)
    !end function OpenRSPSetElecEOM_f

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
            subroutine get_linear_rsp_solution(num_pert,      &
                                               num_comps,     &
                                               num_freq_sums, &
                                               freq_sums,     &
                                               RHS_mat,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                                               len_ctx,       &
                                               user_ctx,      &
#endif
                                               rsp_param)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                integer(kind=QINT), intent(in) :: num_pert
                integer(kind=QINT), intent(in) :: num_comps(num_pert)
                integer(kind=QINT), intent(in) :: num_freq_sums(num_pert)
                real(kind=QREAL), intent(in) :: freq_sums(2*sum(num_freq_sums))
                type(QcMat), intent(in) :: RHS_mat(dot_product(num_comps,num_freq_sums))
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                type(QcMat), intent(inout) :: rsp_param(dot_product(num_comps,num_freq_sums))
            end subroutine get_linear_rsp_solution
            subroutine RSPSolverGetLinearRSPSolution_f(num_pert,      &
                                                       num_comps,     &
                                                       num_freq_sums, &
                                                       freq_sums,     &
                                                       RHS_mat,       &
                                                       user_ctx,      &
                                                       rsp_param)     &
                bind(C, name="RSPSolverGetLinearRSPSolution_f")
                use, intrinsic :: iso_c_binding
                integer(kind=C_QINT), value, intent(in) :: num_pert
                integer(kind=C_QINT), intent(in) :: num_comps(num_pert)
                integer(kind=C_QINT), intent(in) :: num_freq_sums(num_pert)
                real(kind=C_QREAL), intent(in) :: freq_sums(2*sum(num_freq_sums))
                type(C_PTR), intent(in) :: RHS_mat(dot_product(num_comps,num_freq_sums))
                type(C_PTR), value, intent(in) :: user_ctx
                type(C_PTR), intent(in) :: rsp_param(dot_product(num_comps,num_freq_sums))
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

    function OpenRSPSetPerturbations_f(open_rsp,        &
                                       num_pert_lab,    &
                                       pert_labels,     &
                                       pert_max_orders, &
                                       pert_num_comps,  &
#if defined(OPENRSP_F_USER_CONTEXT)
                                       user_ctx,        &
#endif
                                       get_pert_concatenation) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: num_pert_lab
        integer(kind=QcPertInt), intent(in) :: pert_labels(num_pert_lab)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert_lab)
        integer(kind=QINT), intent(in) :: pert_num_comps(:)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_pert_concatenation(pert_label,     &
                                              first_cat_comp, &
                                              num_cat_comps,  &
                                              num_sub_tuples, &
                                              len_sub_tuples, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                              len_ctx,        &
                                              user_ctx,       &
#endif
                                              rank_sub_comps)
                use qcmatrix_f, only: QINT
                use RSPPertBasicTypes_f, only: QcPertInt
                integer(kind=QcPertInt), intent(in) :: pert_label
                integer(kind=QINT), intent(in) :: first_cat_comp
                integer(kind=QINT), intent(in) :: num_cat_comps
                integer(kind=QINT), intent(in) :: num_sub_tuples
                integer(kind=QINT), intent(in) :: len_sub_tuples(num_sub_tuples)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(out) :: rank_sub_comps(num_sub_tuples*num_cat_comps)
            end subroutine get_pert_concatenation
            subroutine RSPPertGetConcatenation_f(pert_label,     &
                                                 first_cat_comp, &
                                                 num_cat_comps,  &
                                                 num_sub_tuples, &
                                                 len_sub_tuples, &
                                                 user_ctx,       &
                                                 rank_sub_comps) &
                bind(C, name="RSPPertGetConcatenation_f")
                use, intrinsic :: iso_c_binding
                use RSPPertBasicTypes_f, only: C_QCPERTINT
                integer(kind=C_QCPERTINT), value, intent(in) :: pert_label
                integer(kind=C_QINT), value, intent(in) :: first_cat_comp
                integer(kind=C_QINT), value, intent(in) :: num_cat_comps
                integer(kind=C_QINT), value, intent(in) :: num_sub_tuples
                integer(kind=C_QINT), intent(in) :: len_sub_tuples(num_sub_tuples)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), intent(out) :: rank_sub_comps(num_sub_tuples*num_cat_comps)
            end subroutine RSPPertGetConcatenation_f
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
                             get_pert_concatenation)
        ierr = OpenRSPSetPerturbations(open_rsp%c_rsp,             &
                                       num_pert_lab,               &
                                       pert_labels,                &
                                       pert_max_orders,            &
                                       pert_num_comps,             &
                                       c_loc(open_rsp%pert_fun),   &
                                       c_funloc(RSPPertGetConcatenation_f))
    end function OpenRSPSetPerturbations_f

    function OpenRSPSetOverlap_f(open_rsp,        &
                                 num_pert_lab,    &
                                 pert_labels,     &
                                 pert_max_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                 user_ctx,        &
#endif
                                 get_overlap_mat, &
                                 get_overlap_exp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: num_pert_lab
        integer(kind=QcPertInt), intent(in) :: pert_labels(num_pert_lab)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert_lab)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_overlap_mat(bra_num_pert,     &
                                       bra_pert_labels,  &
                                       bra_pert_orders,  &
                                       ket_num_pert,     &
                                       ket_pert_labels,  &
                                       ket_pert_orders,  &
                                       oper_num_pert,    &
                                       oper_pert_labels, &
                                       oper_pert_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                       len_ctx,          &
                                       user_ctx,         &
#endif
                                       num_int,          &
                                       val_int)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                use RSPPertBasicTypes_f, only: QcPertInt
                integer(kind=QINT), intent(in) :: bra_num_pert
                integer(kind=QcPertInt), intent(in) :: bra_pert_labels(bra_num_pert)
                integer(kind=QINT), intent(in) :: bra_pert_orders(bra_num_pert)
                integer(kind=QINT), intent(in) :: ket_num_pert
                integer(kind=QcPertInt), intent(in) :: ket_pert_labels(ket_num_pert)
                integer(kind=QINT), intent(in) :: ket_pert_orders(ket_num_pert)
                integer(kind=QINT), intent(in) :: oper_num_pert
                integer(kind=QcPertInt), intent(in) :: oper_pert_labels(oper_num_pert)
                integer(kind=QINT), intent(in) :: oper_pert_orders(oper_num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_int
                type(QcMat), intent(inout) :: val_int(num_int)
            end subroutine get_overlap_mat
            subroutine get_overlap_exp(bra_num_pert,     &
                                       bra_pert_labels,  &
                                       bra_pert_orders,  &
                                       ket_num_pert,     &
                                       ket_pert_labels,  &
                                       ket_pert_orders,  &
                                       oper_num_pert,    &
                                       oper_pert_labels, &
                                       oper_pert_orders, &
                                       num_dmat,         &
                                       dens_mat,         &
#if defined(OPENRSP_F_USER_CONTEXT)
                                       len_ctx,          &
                                       user_ctx,         &
#endif
                                       num_exp,          &
                                       val_exp)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                use RSPPertBasicTypes_f, only: QcPertInt
                integer(kind=QINT), intent(in) :: bra_num_pert
                integer(kind=QcPertInt), intent(in) :: bra_pert_labels(bra_num_pert)
                integer(kind=QINT), intent(in) :: bra_pert_orders(bra_num_pert)
                integer(kind=QINT), intent(in) :: ket_num_pert
                integer(kind=QcPertInt), intent(in) :: ket_pert_labels(ket_num_pert)
                integer(kind=QINT), intent(in) :: ket_pert_orders(ket_num_pert)
                integer(kind=QINT), intent(in) :: oper_num_pert
                integer(kind=QcPertInt), intent(in) :: oper_pert_labels(oper_num_pert)
                integer(kind=QINT), intent(in) :: oper_pert_orders(oper_num_pert)
                integer(kind=QINT), intent(in) :: num_dmat
                type(QcMat), intent(in) :: dens_mat(num_dmat)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_exp
                real(kind=QREAL), intent(inout) :: val_exp(2*num_exp)
            end subroutine get_overlap_exp
            subroutine RSPOverlapGetMat_f(bra_num_pert,     &
                                          bra_pert_labels,  &
                                          bra_pert_orders,  &
                                          ket_num_pert,     &
                                          ket_pert_labels,  &
                                          ket_pert_orders,  &
                                          oper_num_pert,    &
                                          oper_pert_labels, &
                                          oper_pert_orders, &
                                          user_ctx,         &
                                          num_int,          &
                                          val_int)          &
                bind(C, name="RSPOverlapGetMat_f")
                use, intrinsic :: iso_c_binding
                use RSPPertBasicTypes_f, only: C_QCPERTINT
                integer(kind=C_QINT), value, intent(in) :: bra_num_pert
                integer(kind=C_QCPERTINT), intent(in) :: bra_pert_labels(bra_num_pert)
                integer(kind=C_QINT), intent(in) :: bra_pert_orders(bra_num_pert)
                integer(kind=C_QINT), value, intent(in) :: ket_num_pert
                integer(kind=C_QCPERTINT), intent(in) :: ket_pert_labels(ket_num_pert)
                integer(kind=C_QINT), intent(in) :: ket_pert_orders(ket_num_pert)
                integer(kind=C_QINT), value, intent(in) :: oper_num_pert
                integer(kind=C_QCPERTINT), intent(in) :: oper_pert_labels(oper_num_pert)
                integer(kind=C_QINT), intent(in) :: oper_pert_orders(oper_num_pert)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_int
                type(C_PTR), intent(inout) :: val_int(num_int)
            end subroutine RSPOverlapGetMat_f
            subroutine RSPOverlapGetExp_f(bra_num_pert,     &
                                          bra_pert_labels,  &
                                          bra_pert_orders,  &
                                          ket_num_pert,     &
                                          ket_pert_labels,  &
                                          ket_pert_orders,  &
                                          oper_num_pert,    &
                                          oper_pert_labels, &
                                          oper_pert_orders, &
                                          num_dmat,         &
                                          dens_mat,         &
                                          user_ctx,         &
                                          num_exp,          &
                                          val_exp)          &
                bind(C, name="RSPOverlapGetExp_f")
                use, intrinsic :: iso_c_binding
                use RSPPertBasicTypes_f, only: C_QCPERTINT
                integer(kind=C_QINT), value, intent(in) :: bra_num_pert
                integer(kind=C_QCPERTINT), intent(in) :: bra_pert_labels(bra_num_pert)
                integer(kind=C_QINT), intent(in) :: bra_pert_orders(bra_num_pert)
                integer(kind=C_QINT), value, intent(in) :: ket_num_pert
                integer(kind=C_QCPERTINT), intent(in) :: ket_pert_labels(ket_num_pert)
                integer(kind=C_QINT), intent(in) :: ket_pert_orders(ket_num_pert)
                integer(kind=C_QINT), value, intent(in) :: oper_num_pert
                integer(kind=C_QCPERTINT), intent(in) :: oper_pert_labels(oper_num_pert)
                integer(kind=C_QINT), intent(in) :: oper_pert_orders(oper_num_pert)
                integer(kind=C_QINT), value, intent(in) :: num_dmat
                type(C_PTR), intent(in) :: dens_mat(num_dmat)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_exp
                real(kind=C_QREAL), intent(inout) :: val_exp(2*num_exp)
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
        ierr = OpenRSPSetOverlap(open_rsp%c_rsp,               &
                                 num_pert_lab,                 &
                                 pert_labels,                  &
                                 pert_max_orders,              &
                                 c_loc(open_rsp%overlap_fun),  &
                                 c_funloc(RSPOverlapGetMat_f), &
                                 c_funloc(RSPOverlapGetExp_f))
    end function OpenRSPSetOverlap_f

    function OpenRSPAddOneOper_f(open_rsp,         &
                                 num_pert_lab,     &
                                 pert_labels,      &
                                 pert_max_orders,  &
#if defined(OPENRSP_F_USER_CONTEXT)
                                 user_ctx,         &
#endif
                                 get_one_oper_mat, &
                                 get_one_oper_exp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: num_pert_lab
        integer(kind=QcPertInt), intent(in) :: pert_labels(num_pert_lab)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert_lab)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_one_oper_mat(oper_num_pert,    &
                                        oper_pert_labels, &
                                        oper_pert_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                        len_ctx,          &
                                        user_ctx,         &
#endif
                                        num_int,          &
                                        val_int)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                use RSPPertBasicTypes_f, only: QcPertInt
                integer(kind=QINT), intent(in) :: oper_num_pert
                integer(kind=QcPertInt), intent(in) :: oper_pert_labels(oper_num_pert)
                integer(kind=QINT), intent(in) :: oper_pert_orders(oper_num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_int
                type(QcMat), intent(inout) :: val_int(num_int)
            end subroutine get_one_oper_mat
            subroutine get_one_oper_exp(oper_num_pert,    &
                                        oper_pert_labels, &
                                        oper_pert_orders, &
                                        num_dmat,         &
                                        dens_mat,         &
#if defined(OPENRSP_F_USER_CONTEXT)
                                        len_ctx,          &
                                        user_ctx,         &
#endif
                                        num_exp,          &
                                        val_exp)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                use RSPPertBasicTypes_f, only: QcPertInt
                integer(kind=QINT), intent(in) :: oper_num_pert
                integer(kind=QcPertInt), intent(in) :: oper_pert_labels(oper_num_pert)
                integer(kind=QINT), intent(in) :: oper_pert_orders(oper_num_pert)
                integer(kind=QINT), intent(in) :: num_dmat
                type(QcMat), intent(in) :: dens_mat(num_dmat)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_exp
                real(kind=QREAL), intent(inout) :: val_exp(2*num_exp)
            end subroutine get_one_oper_exp
            subroutine RSPOneOperGetMat_f(oper_num_pert,    &
                                          oper_pert_labels, &
                                          oper_pert_orders, &
                                          user_ctx,         &
                                          num_int,          &
                                          val_int)          &
                bind(C, name="RSPOneOperGetMat_f")
                use, intrinsic :: iso_c_binding
                use RSPPertBasicTypes_f, only: C_QCPERTINT
                integer(kind=C_QINT), value, intent(in) :: oper_num_pert
                integer(kind=C_QCPERTINT), intent(in) :: oper_pert_labels(oper_num_pert)
                integer(kind=C_QINT), intent(in) :: oper_pert_orders(oper_num_pert)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_int
                type(C_PTR), intent(inout) :: val_int(num_int)
            end subroutine RSPOneOperGetMat_f
            subroutine RSPOneOperGetExp_f(oper_num_pert,    &
                                          oper_pert_labels, &
                                          oper_pert_orders, &
                                          num_dmat,         &
                                          dens_mat,         &
                                          user_ctx,         &
                                          num_exp,          &
                                          val_exp)          &
                bind(C, name="RSPOneOperGetExp_f")
                use, intrinsic :: iso_c_binding
                use RSPPertBasicTypes_f, only: C_QCPERTINT
                integer(kind=C_QINT), value, intent(in) :: oper_num_pert
                integer(kind=C_QCPERTINT), intent(in) :: oper_pert_labels(oper_num_pert)
                integer(kind=C_QINT), intent(in) :: oper_pert_orders(oper_num_pert)
                integer(kind=C_QINT), value, intent(in) :: num_dmat
                type(C_PTR), intent(in) :: dens_mat(num_dmat)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_exp
                real(kind=C_QREAL), intent(inout) :: val_exp(2*num_exp)
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
                                 num_pert_lab,                     &
                                 pert_labels,                      &
                                 pert_max_orders,                  &
                                 c_loc(cur_one_oper%one_oper_fun), &
                                 c_funloc(RSPOneOperGetMat_f),     &
                                 c_funloc(RSPOneOperGetExp_f))
    end function OpenRSPAddOneOper_f

    function OpenRSPAddTwoOper_f(open_rsp,         &
                                 num_pert_lab,     &
                                 pert_labels,      &
                                 pert_max_orders,  &
#if defined(OPENRSP_F_USER_CONTEXT)
                                 user_ctx,         &
#endif
                                 get_two_oper_mat, &
                                 get_two_oper_exp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: num_pert_lab
        integer(kind=QcPertInt), intent(in) :: pert_labels(num_pert_lab)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert_lab)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_two_oper_mat(oper_num_pert,    &
                                        oper_pert_labels, &
                                        oper_pert_orders, &
                                        num_dmat,         &
                                        dens_mat,         &
#if defined(OPENRSP_F_USER_CONTEXT)
                                        len_ctx,          &
                                        user_ctx,         &
#endif
                                        num_int,          &
                                        val_int)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                use RSPPertBasicTypes_f, only: QcPertInt
                integer(kind=QINT), intent(in) :: oper_num_pert
                integer(kind=QcPertInt), intent(in) :: oper_pert_labels(oper_num_pert)
                integer(kind=QINT), intent(in) :: oper_pert_orders(oper_num_pert)
                integer(kind=QINT), intent(in) :: num_dmat
                type(QcMat), intent(in) :: dens_mat(num_dmat)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_int
                type(QcMat), intent(inout) :: val_int(num_int)
            end subroutine get_two_oper_mat
            subroutine get_two_oper_exp(oper_num_pert,    &
                                        oper_pert_labels, &
                                        oper_pert_orders, &
                                        dmat_len_tuple,   &
                                        num_LHS_dmat,     &
                                        LHS_dens_mat,     &
                                        num_RHS_dmat,     &
                                        RHS_dens_mat,     &
#if defined(OPENRSP_F_USER_CONTEXT)
                                        len_ctx,          &
                                        user_ctx,         &
#endif
                                        num_exp,          &
                                        val_exp)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                use RSPPertBasicTypes_f, only: QcPertInt
                integer(kind=QINT), intent(in) :: oper_num_pert
                integer(kind=QcPertInt), intent(in) :: oper_pert_labels(oper_num_pert)
                integer(kind=QINT), intent(in) :: oper_pert_orders(oper_num_pert)
                integer(kind=QINT), intent(in) :: dmat_len_tuple
                integer(kind=QINT), intent(in) :: num_LHS_dmat(dmat_len_tuple)
                type(QcMat), intent(in) :: LHS_dens_mat(sum(num_LHS_dmat))
                integer(kind=QINT), intent(in) :: num_RHS_dmat(dmat_len_tuple)
                type(QcMat), intent(in) :: RHS_dens_mat(sum(num_RHS_dmat))
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_exp
                real(kind=QREAL), intent(inout) :: val_exp(2*num_exp)
            end subroutine get_two_oper_exp
            subroutine RSPTwoOperGetMat_f(oper_num_pert,    &
                                          oper_pert_labels, &
                                          oper_pert_orders, &
                                          num_dmat,         &
                                          dens_mat,         &
                                          user_ctx,         &
                                          num_int,          &
                                          val_int)          &
                bind(C, name="RSPTwoOperGetMat_f")
                use, intrinsic :: iso_c_binding
                use RSPPertBasicTypes_f, only: C_QCPERTINT
                integer(kind=C_QINT), value, intent(in) :: oper_num_pert
                integer(kind=C_QCPERTINT), intent(in) :: oper_pert_labels(oper_num_pert)
                integer(kind=C_QINT), intent(in) :: oper_pert_orders(oper_num_pert)
                integer(kind=C_QINT), value, intent(in) :: num_dmat
                type(C_PTR), intent(in) :: dens_mat(num_dmat)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_int
                type(C_PTR), intent(inout) :: val_int(num_int)
            end subroutine RSPTwoOperGetMat_f
            subroutine RSPTwoOperGetExp_f(oper_num_pert,    &
                                          oper_pert_labels, &
                                          oper_pert_orders, &
                                          dmat_len_tuple,   &
                                          num_LHS_dmat,     &
                                          LHS_dens_mat,     &
                                          num_RHS_dmat,     &
                                          RHS_dens_mat,     &
                                          user_ctx,         &
                                          num_exp,          &
                                          val_exp)          &
                bind(C, name="RSPTwoOperGetExp_f")
                use, intrinsic :: iso_c_binding
                use RSPPertBasicTypes_f, only: C_QCPERTINT
                integer(kind=C_QINT), value, intent(in) :: oper_num_pert
                integer(kind=C_QCPERTINT), intent(in) :: oper_pert_labels(oper_num_pert)
                integer(kind=C_QINT), intent(in) :: oper_pert_orders(oper_num_pert)
                integer(kind=C_QINT), value, intent(in) :: dmat_len_tuple
                integer(kind=C_QINT), intent(in) :: num_LHS_dmat(dmat_len_tuple)
                type(C_PTR), intent(in) :: LHS_dens_mat(sum(num_LHS_dmat))
                integer(kind=C_QINT), intent(in) :: num_RHS_dmat(dmat_len_tuple)
                type(C_PTR), intent(in) :: RHS_dens_mat(sum(num_RHS_dmat))
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_exp
                real(kind=C_QREAL), intent(inout) :: val_exp(2*num_exp)
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
                                 num_pert_lab,                     &
                                 pert_labels,                      &
                                 pert_max_orders,                  &
                                 c_loc(cur_two_oper%two_oper_fun), &
                                 c_funloc(RSPTwoOperGetMat_f),     &
                                 c_funloc(RSPTwoOperGetExp_f))
    end function OpenRSPAddTwoOper_f

    function OpenRSPAddXCFun_f(open_rsp,        &
                               num_pert_lab,    &
                               pert_labels,     &
                               pert_max_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                               user_ctx,        &
#endif
                               get_xc_fun_mat,  &
                               get_xc_fun_exp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: num_pert_lab
        integer(kind=QcPertInt), intent(in) :: pert_labels(num_pert_lab)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert_lab)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_xc_fun_mat(xc_len_tuple,       &
                                      xc_pert_tuple,      &
                                      num_freq_configs,   &
                                      pert_freq_category, &
                                      dmat_num_tuple,     &
                                      dmat_idx_tuple,     &
                                      num_dmat,           &
                                      dens_mat,           &
#if defined(OPENRSP_F_USER_CONTEXT)
                                      len_ctx,            &
                                      user_ctx,           &
#endif
                                      num_int,            &
                                      val_int)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                use RSPPertBasicTypes_f, only: QcPertInt
                integer(kind=QINT), intent(in) :: xc_len_tuple
                integer(kind=QcPertInt), intent(in) :: xc_pert_tuple(xc_len_tuple)
                integer(kind=QINT), intent(in) :: num_freq_configs
                integer(kind=QINT), intent(in) :: &
                    pert_freq_category(num_freq_configs*xc_len_tuple)
                integer(kind=QINT), intent(in) :: dmat_num_tuple
                integer(kind=QINT), intent(in) :: dmat_idx_tuple(dmat_num_tuple)
                integer(kind=QINT), intent(in) :: num_dmat
                type(QcMat), intent(in) :: dens_mat(num_dmat)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_int
                type(QcMat), intent(inout) :: val_int(num_int)
            end subroutine get_xc_fun_mat
            subroutine get_xc_fun_exp(xc_len_tuple,       &
                                      xc_pert_tuple,      &
                                      num_freq_configs,   &
                                      pert_freq_category, &
                                      dmat_num_tuple,     &
                                      dmat_idx_tuple,     &
                                      num_dmat,           &
                                      dens_mat,           &
#if defined(OPENRSP_F_USER_CONTEXT)
                                      len_ctx,            &
                                      user_ctx,           &
#endif
                                      num_exp,            &
                                      val_exp)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                use RSPPertBasicTypes_f, only: QcPertInt
                integer(kind=QINT), intent(in) :: xc_len_tuple
                integer(kind=QcPertInt), intent(in) :: xc_pert_tuple(xc_len_tuple)
                integer(kind=QINT), intent(in) :: num_freq_configs
                integer(kind=QINT), intent(in) :: &
                    pert_freq_category(num_freq_configs*xc_len_tuple)
                integer(kind=QINT), intent(in) :: dmat_num_tuple
                integer(kind=QINT), intent(in) :: dmat_idx_tuple(dmat_num_tuple)
                integer(kind=QINT), intent(in) :: num_dmat
                type(QcMat), intent(in) :: dens_mat(num_dmat)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_exp
                real(kind=QREAL), intent(inout) :: val_exp(2*num_exp)
            end subroutine get_xc_fun_exp
            subroutine RSPXCFunGetMat_f(xc_len_tuple,       &
                                        xc_pert_tuple,      &
                                        num_freq_configs,   &
                                        pert_freq_category, &
                                        dmat_num_tuple,     &
                                        dmat_idx_tuple,     &
                                        num_dmat,           &
                                        dens_mat,           &
                                        user_ctx,           &
                                        num_int,            &
                                        val_int)            &
                bind(C, name="RSPXCFunGetMat_f")
                use, intrinsic :: iso_c_binding
                use RSPPertBasicTypes_f, only: C_QCPERTINT
                integer(kind=C_QINT), value, intent(in) :: xc_len_tuple
                integer(kind=C_QCPERTINT), intent(in) :: xc_pert_tuple(xc_len_tuple)
                integer(kind=C_QINT), value, intent(in) :: num_freq_configs
                integer(kind=C_QINT), intent(in) :: &
                    pert_freq_category(num_freq_configs*xc_len_tuple)
                integer(kind=C_QINT), value, intent(in) :: dmat_num_tuple
                integer(kind=C_QINT), intent(in) :: dmat_idx_tuple(dmat_num_tuple)
                integer(kind=C_QINT), value, intent(in) :: num_dmat
                type(C_PTR), intent(in) :: dens_mat(num_dmat)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_int
                type(C_PTR), intent(inout) :: val_int(num_int)
            end subroutine RSPXCFunGetMat_f
            subroutine RSPXCFunGetExp_f(xc_len_tuple,       &
                                        xc_pert_tuple,      &
                                        num_freq_configs,   &
                                        pert_freq_category, &
                                        dmat_num_tuple,     &
                                        dmat_idx_tuple,     &
                                        num_dmat,           &
                                        dens_mat,           &
                                        user_ctx,           &
                                        num_exp,            &
                                        val_exp)            &
                bind(C, name="RSPXCFunGetExp_f")
                use, intrinsic :: iso_c_binding
                use RSPPertBasicTypes_f, only: C_QCPERTINT
                integer(kind=C_QINT), value, intent(in) :: xc_len_tuple
                integer(kind=C_QCPERTINT), intent(in) :: xc_pert_tuple(xc_len_tuple)
                integer(kind=C_QINT), value, intent(in) :: num_freq_configs
                integer(kind=C_QINT), intent(in) :: &
                    pert_freq_category(num_freq_configs*xc_len_tuple)
                integer(kind=C_QINT), value, intent(in) :: dmat_num_tuple
                integer(kind=C_QINT), intent(in) :: dmat_idx_tuple(dmat_num_tuple)
                integer(kind=C_QINT), value, intent(in) :: num_dmat
                type(C_PTR), intent(in) :: dens_mat(num_dmat)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: num_exp
                real(kind=C_QREAL), intent(inout) :: val_exp(2*num_exp)
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
        allocate(cur_xc_fun%xcfun_fun)
        nullify(cur_xc_fun%next_xc_fun)
        ! adds context of callback functions of the new XC functional
        call RSPXCFunCreate_f(cur_xc_fun%xcfun_fun, &
#if defined(OPENRSP_F_USER_CONTEXT)
                              user_ctx,             &
#endif
                              get_xc_fun_mat,       &
                              get_xc_fun_exp)
        ierr = OpenRSPAddXCFun(open_rsp%c_rsp,              &
                               num_pert_lab,                &
                               pert_labels,                 &
                               pert_max_orders,             &
                               c_loc(cur_xc_fun%xcfun_fun), &
                               c_funloc(RSPXCFunGetMat_f),  &
                               c_funloc(RSPXCFunGetExp_f))
    end function OpenRSPAddXCFun_f

    function OpenRSPSetNucHamilton_f(open_rsp,        &
                                     num_pert_lab,    &
                                     pert_labels,     &
                                     pert_max_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                     user_ctx,        &
#endif
                                     get_nuc_contrib, &
                                     num_atoms) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        integer(kind=QINT), intent(in) :: num_pert_lab
        integer(kind=QcPertInt), intent(in) :: pert_labels(num_pert_lab)
        integer(kind=QINT), intent(in) :: pert_max_orders(num_pert_lab)
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        integer(kind=QINT), intent(in) :: num_atoms
        interface
            subroutine get_nuc_contrib(nuc_num_pert,    &
                                       nuc_pert_labels, &
                                       nuc_pert_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                       len_ctx,         &
                                       user_ctx,        &
#endif
                                       size_pert,       &
                                       val_nuc)
                use qcmatrix_f, only: QINT,QREAL
                use RSPPertBasicTypes_f, only: QcPertInt
                integer(kind=QINT), intent(in) :: nuc_num_pert
                integer(kind=QcPertInt), intent(in) :: nuc_pert_labels(nuc_num_pert)
                integer(kind=QINT), intent(in) :: nuc_pert_orders(nuc_num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: size_pert
                real(kind=QREAL), intent(inout) :: val_nuc(2*size_pert)
            end subroutine get_nuc_contrib
            subroutine RSPNucHamiltonGetContributions_f(nuc_num_pert,    &
                                                        nuc_pert_labels, &
                                                        nuc_pert_orders, &
                                                        user_ctx,        &
                                                        size_pert,       &
                                                        val_nuc)         &
                bind(C, name="RSPNucHamiltonGetContributions_f")
                use, intrinsic :: iso_c_binding
                use RSPPertBasicTypes_f, only: C_QCPERTINT
                integer(kind=C_QINT), value, intent(in) :: nuc_num_pert
                integer(kind=C_QCPERTINT), intent(in) :: nuc_pert_labels(nuc_num_pert)
                integer(kind=C_QINT), intent(in) :: nuc_pert_orders(nuc_num_pert)
                type(C_PTR), value, intent(in) :: user_ctx
                integer(kind=C_QINT), value, intent(in) :: size_pert
                real(kind=C_QREAL), intent(inout) :: val_nuc(2*size_pert)
            end subroutine RSPNucHamiltonGetContributions_f
        end interface
        if (associated(open_rsp%nuc_hamilton_fun)) then
            call RSPNucHamiltonDestroy_f(open_rsp%nuc_hamilton_fun)
        else
            allocate(open_rsp%nuc_hamilton_fun)
        end if
        ! adds context of callback function of the nuclear Hamiltonian
        call RSPNucHamiltonCreate_f(open_rsp%nuc_hamilton_fun, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                    user_ctx,                  &
#endif
                                    get_nuc_contrib)
        ierr = OpenRSPSetNucHamilton(open_rsp%c_rsp,                             &
                                     num_pert_lab,                               &
                                     pert_labels,                                &
                                     pert_max_orders,                            &
                                     c_loc(open_rsp%nuc_hamilton_fun),           &
                                     c_funloc(RSPNucHamiltonGetContributions_f), &
                                     num_atoms)
    end function OpenRSPSetNucHamilton_f

    function OpenRSPAssemble_f(open_rsp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        ierr = OpenRSPAssemble(open_rsp%c_rsp)
    end function OpenRSPAssemble_f

    function OpenRSPWritebyFileName_f(open_rsp, file_name) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(in) :: open_rsp
        character*(*), intent(in) :: file_name
        ierr = OpenRSPWriteFortranAdapter(open_rsp%c_rsp, &
                                          file_name//C_NULL_CHAR)
    end function OpenRSPWritebyFileName_f

    function OpenRSPWritebyUnit_f(open_rsp, io_unit) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(in) :: open_rsp
        integer(kind=4), intent(in) :: io_unit
integer(kind=4), parameter :: MAX_LEN_QCCHAR = 256
integer(kind=4), parameter :: QCSTDOUT = 6
        logical(kind=4) is_opened
        logical(kind=4) is_named
        character(MAX_LEN_QCCHAR) file_name
        if (io_unit==QCSTDOUT) then
            ierr = OpenRSPWriteStdOutFortranAdapter(open_rsp%c_rsp)
        else
            ! gets the file name from its associated logical unit
            inquire(unit=io_unit,     &
                    opened=is_opened, &
                    named=is_named,   &
                    name=file_name)
            if (is_opened .and. is_named) then
                if (len_trim(file_name)==MAX_LEN_QCCHAR) then
                    write(QCSTDOUT,100) "file name too long, increase MAX_LEN_QCCHAR", &
                                        MAX_LEN_QCCHAR
                    call QErrorExit(QCSTDOUT, __LINE__, "OpenRSP.F90")
                else
                    ! writes by file name
                    ierr = OpenRSPWriteFortranAdapter(open_rsp%c_rsp, &
                                                      trim(file_name)//C_NULL_CHAR)
                end if
            else
                write(QCSTDOUT,100) "logical unit is not named and/or opened", &
                                    io_unit
                call QErrorExit(QCSTDOUT, __LINE__, "OpenRSP.F90")
            end if
        end if
100     format("OpenRSPWritebyUnit_f>> ",A,I6)
    end function OpenRSPWritebyUnit_f

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
        integer(kind=QcPertInt), intent(in) :: pert_tuple(sum(len_tuple))
        integer(kind=QINT), intent(in) :: num_freq_configs(num_props)
        real(kind=QREAL), intent(in) :: &
            pert_freqs(2*dot_product(len_tuple,num_freq_configs))
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

    function OpenRSPGetResidue_f(open_rsp,         &
                                 ref_ham,          &
                                 ref_state,        &
                                 ref_overlap,      &
                                 order_residue,    &
                                 num_excit,        &
                                 excit_energy,     &
                                 eigen_vector,     &
                                 num_props,        &
                                 len_tuple,        &
                                 pert_tuple,       &
                                 residue_num_pert, &
                                 residue_idx_pert, &
                                 num_freq_configs, &
                                 pert_freqs,       &
                                 kn_rules,         &
                                 size_residues,    &
                                 residues) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(in) :: open_rsp
        type(QcMat), target, intent(in) :: ref_ham
        type(QcMat), target, intent(in) :: ref_state
        type(QcMat), target, intent(in) :: ref_overlap
        integer(kind=QINT), intent(in) :: order_residue
        integer(kind=QINT), intent(in) :: num_excit
        real(kind=QREAL), intent(in) :: excit_energy(order_residue*num_excit)
        type(QcMat), intent(in) :: eigen_vector(order_residue*num_excit)
        integer(kind=QINT), intent(in) :: num_props
        integer(kind=QINT), intent(in) :: len_tuple(num_props)
        integer(kind=QcPertInt), intent(in) :: pert_tuple(sum(len_tuple))
        integer(kind=QINT), intent(in) :: residue_num_pert(order_residue*num_props)
        integer(kind=QINT), intent(in) :: residue_idx_pert(sum(residue_num_pert))
        integer(kind=QINT), intent(in) :: num_freq_configs(num_props)
        real(kind=QREAL), intent(in) :: &
            pert_freqs(2*dot_product(len_tuple,num_freq_configs)*num_excit)
        integer(kind=QINT), intent(in) :: kn_rules(num_props)
        integer(kind=QINT), intent(in) :: size_residues
        real(kind=QREAL), intent(out) :: residues(2*size_residues)
        type(C_PTR) c_ref_ham(1)
        type(C_PTR) c_ref_state(1)
        type(C_PTR) c_ref_overlap(1)
        type(C_PTR), allocatable :: c_eigen_vector(:)
        integer(kind=QINT) iext
        ierr = QcMat_C_LOC((/ref_ham/), c_ref_ham)
        if (ierr==QFAILURE) return
        ierr = QcMat_C_LOC((/ref_state/), c_ref_state)
        if (ierr==QFAILURE) return
        ierr = QcMat_C_LOC((/ref_overlap/), c_ref_overlap)
        if (ierr==QFAILURE) return
        allocate(c_eigen_vector(order_residue*num_excit), stat=ierr)
        if (ierr/=0) return
        ierr = QcMat_C_LOC(eigen_vector, c_eigen_vector)
        if (ierr==QFAILURE) return
        ierr = OpenRSPGetResidue(open_rsp%c_rsp,    &
                                 c_ref_ham(1),      &
                                 c_ref_state(1),    &
                                 c_ref_overlap(1),  &
                                 order_residue,     &
                                 num_excit,         &
                                 excit_energy,      &
                                 c_eigen_vector,    &
                                 num_props,         &
                                 len_tuple,         &
                                 pert_tuple,        &
                                 residue_num_pert,  &
                                 residue_idx_pert,  &
                                 num_freq_configs,  &
                                 pert_freqs,        &
                                 kn_rules,          &
                                 size_residues,     &
                                 residues)
        c_ref_ham(1) = C_NULL_PTR
        c_ref_state(1) = C_NULL_PTR
        c_ref_overlap(1) = C_NULL_PTR
        do iext = 1, order_residue*num_excit
            c_eigen_vector(iext) = C_NULL_PTR
        end do
        deallocate(c_eigen_vector)
    end function OpenRSPGetResidue_f

    function OpenRSPDestroy_f(open_rsp) result(ierr)
        integer(kind=4) :: ierr
        type(OpenRSP), intent(inout) :: open_rsp
        type(OneOperList_f), pointer :: cur_one_oper   !current one-electron operator
        type(OneOperList_f), pointer :: next_one_oper  !next one-electron operator
        type(TwoOperList_f), pointer :: cur_two_oper   !current two-electron operator
        type(TwoOperList_f), pointer :: next_two_oper  !next two-electron operator
        type(XCFunList_f), pointer :: cur_xc_fun       !current XC functional
        type(XCFunList_f), pointer :: next_xc_fun      !next XC functional
        ierr = OpenRSPDestroyFortranAdapter(open_rsp%c_rsp)
        ! cleans up callback subroutine of response equation solver
        if (associated(open_rsp%solver_fun)) then
            call RSPSolverDestroy_f(open_rsp%solver_fun)
            deallocate(open_rsp%solver_fun)
            nullify(open_rsp%solver_fun)
        end if
        ! cleans up callback subroutines of perturbations
        if (associated(open_rsp%pert_fun)) then
            call RSPPertDestroy_f(open_rsp%pert_fun)
            deallocate(open_rsp%pert_fun)
            nullify(open_rsp%pert_fun)
        end if
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
            if (associated(cur_xc_fun%xcfun_fun)) then
                call RSPXCFunDestroy_f(cur_xc_fun%xcfun_fun)
                deallocate(cur_xc_fun%xcfun_fun)
                nullify(cur_xc_fun%xcfun_fun)
            end if
            deallocate(cur_xc_fun)
            nullify(cur_xc_fun)
            cur_xc_fun => next_xc_fun
        end do
        ! cleans up callback subroutine of nuclear Hamiltonian
        if (associated(open_rsp%nuc_hamilton_fun)) then
            call RSPNucHamiltonDestroy_f(open_rsp%nuc_hamilton_fun)
            deallocate(open_rsp%nuc_hamilton_fun)
            nullify(open_rsp%nuc_hamilton_fun)
        end if
    end function OpenRSPDestroy_f

end module OpenRSP_f

