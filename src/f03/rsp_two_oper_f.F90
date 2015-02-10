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
!!  This file implements subroutines of two-electron operators used in the Fortran APIs.
!!
!!  2014-08-06, Bin Gao
!!  * first version

! basic data types
#include "api/qcmatrix_c_type.h"

#define OPENRSP_API_SRC "src/f03/rsp_two_oper_f.F90"

module rsp_two_oper_f

    use, intrinsic :: iso_c_binding
    use qcmatrix_f, only: QINT,QREAL,QcMat,QcMat_C_F_POINTER,QcMat_C_NULL_PTR

    implicit none

    integer(kind=4), private, parameter :: STDOUT = 6

    ! user specified callback subroutines
    abstract interface
        subroutine TwoOperGetMat_f(num_pert,     &
                                   pert_labels,  &
                                   pert_orders,  &
                                   num_var_dens, &
                                   var_ao_dens,  &
#if defined(OPENRSP_F_USER_CONTEXT)
                                   len_ctx,      &
                                   user_ctx,     &
#endif
                                   num_int,      &
                                   val_int)
            use qcmatrix_f, only: QINT,QREAL,QcMat
            integer(kind=QINT), intent(in) :: num_pert
            integer(kind=QINT), intent(in) :: pert_labels(num_pert)
            integer(kind=QINT), intent(in) :: pert_orders(num_pert)
            integer(kind=QINT), intent(in) :: num_var_dens
            type(QcMat), intent(in) :: var_ao_dens(num_var_dens)
#if defined(OPENRSP_F_USER_CONTEXT)
            integer(kind=QINT), intent(in) :: len_ctx
            character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
            integer(kind=QINT), intent(in) :: num_int
            type(QcMat), intent(inout) :: val_int(num_int)
        end subroutine TwoOperGetMat_f
        subroutine TwoOperGetExp_f(num_pert,       &
                                   pert_labels,    &
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
            use qcmatrix_f, only: QINT,QREAL,QcMat
            integer(kind=QINT), intent(in) :: num_pert
            integer(kind=QINT), intent(in) :: pert_labels(num_pert)
            integer(kind=QINT), intent(in) :: pert_orders(num_pert)
            integer(kind=QINT), intent(in) :: num_var_dens
            type(QcMat), intent(in) :: var_ao_dens(num_var_dens)
            integer(kind=QINT), intent(in) :: num_contr_dens
            type(QcMat), intent(in) :: contr_ao_dens(num_contr_dens)
#if defined(OPENRSP_F_USER_CONTEXT)
            integer(kind=QINT), intent(in) :: len_ctx
            character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
            integer(kind=QINT), intent(in) :: num_exp
            real(kind=QREAL), intent(inout) :: val_exp(num_exp)
        end subroutine TwoOperGetExp_f
    end interface

    ! context of callback subroutines of two-electron operator
    type, public :: TwoOperFun_f
        private
#if defined(OPENRSP_F_USER_CONTEXT)
        ! user-defined callback function context
        integer(kind=QINT) :: len_ctx = 0
        character(len=1), allocatable :: user_ctx(:)
#endif
        ! callback functions
        procedure(TwoOperGetMat_f), nopass, pointer :: get_two_oper_mat
        procedure(TwoOperGetExp_f), nopass, pointer :: get_two_oper_exp
    end type TwoOperFun_f

    public :: RSPTwoOperCreate_f
    public :: RSPTwoOperGetMat_f
    public :: RSPTwoOperGetExp_f
    public :: RSPTwoOperDestroy_f

    contains

    !% \brief creates the context of callback subroutines of two-electron operator
    !  \author Bin Gao
    !  \date 2014-08-06
    !  \param[TwoOperFun_f:type]{inout} two_oper_fun the context of callback subroutines
    !  \param[character]{in} user_ctx user-defined callback function context
    !  \param[subroutine]{in} get_two_oper_mat user specified function for
    !      getting integral matrices
    !  \param[subroutine]{in} get_two_oper_exp user specified function for
    !%     getting expectation values
    subroutine RSPTwoOperCreate_f(two_oper_fun,     &
#if defined(OPENRSP_F_USER_CONTEXT)
                                  user_ctx,         &
#endif
                                  get_two_oper_mat, &
                                  get_two_oper_exp)
        type(TwoOperFun_f), intent(inout) :: two_oper_fun
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_two_oper_mat(num_pert,     &
                                        pert_labels,  &
                                        pert_orders,  &
                                        num_var_dens, &
                                        var_ao_dens,  &
#if defined(OPENRSP_F_USER_CONTEXT)
                                        len_ctx,      &
                                        user_ctx,     &
#endif
                                        num_int,      &
                                        val_int)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                integer(kind=QINT), intent(in) :: num_pert
                integer(kind=QINT), intent(in) :: pert_labels(num_pert)
                integer(kind=QINT), intent(in) :: pert_orders(num_pert)
                integer(kind=QINT), intent(in) :: num_var_dens
                type(QcMat), intent(in) :: var_ao_dens(num_var_dens)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_int
                type(QcMat), intent(inout) :: val_int(num_int)
            end subroutine get_two_oper_mat
            subroutine get_two_oper_exp(num_pert,       &
                                        pert_labels,    &
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
                use qcmatrix_f, only: QINT,QREAL,QcMat
                integer(kind=QINT), intent(in) :: num_pert
                integer(kind=QINT), intent(in) :: pert_labels(num_pert)
                integer(kind=QINT), intent(in) :: pert_orders(num_pert)
                integer(kind=QINT), intent(in) :: num_var_dens
                type(QcMat), intent(in) :: var_ao_dens(num_var_dens)
                integer(kind=QINT), intent(in) :: num_contr_dens
                type(QcMat), intent(in) :: contr_ao_dens(num_contr_dens)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_exp
                real(kind=QREAL), intent(inout) :: val_exp(num_exp)
            end subroutine get_two_oper_exp
        end interface
#if defined(OPENRSP_F_USER_CONTEXT)
        integer(kind=4) ierr  !error information
        two_oper_fun%len_ctx = size(user_ctx)
        allocate(two_oper_fun%user_ctx(two_oper_fun%len_ctx), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPTwoOperCreate_f>> length", two_oper_fun%len_ctx
            stop "RSPTwoOperCreate_f>> failed to allocate memory for user_ctx"
        end if
        two_oper_fun%user_ctx = user_ctx
#endif
        two_oper_fun%get_two_oper_mat => get_two_oper_mat
        two_oper_fun%get_two_oper_exp => get_two_oper_exp
    end subroutine RSPTwoOperCreate_f

    !% \brief calls Fortran callback subroutine to get integral matrices of
    !      a two-electron operator
    !  \author Bin Gao
    !  \date 2014-08-06
    !  \param[integer]{in} num_pert number of perturbations
    !  \param[integer]{in} pert_labels labels of the perturbations
    !  \param[integer]{in} pert_orders orders of the perturbations
    !  \param[integer]{in} num_var_dens number of variable AO based density matrices
    !  \param[C_PTR:type]{in} var_ao_dens the variable AO based density matrices
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !  \param[integer]{in} num_int number of the integral matrices
    !% \param[C_PTR:type]{inout} val_int the integral matrices
    subroutine RSPTwoOperGetMat_f(num_pert,     &
                                  pert_labels,  &
                                  pert_orders,  &
                                  num_var_dens, &
                                  var_ao_dens,  &
                                  user_ctx,     &
                                  num_int,      &
                                  val_int)      &
        bind(C, name="RSPTwoOperGetMat_f")
        integer(kind=C_QINT), value, intent(in) :: num_pert
        integer(kind=C_QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
        integer(kind=C_QINT), value, intent(in) :: num_var_dens
        type(C_PTR), intent(in) :: var_ao_dens(num_var_dens)
        type(C_PTR), value, intent(in) :: user_ctx
        integer(kind=C_QINT), value, intent(in) :: num_int
        type(C_PTR), intent(inout) :: val_int(num_int)
        type(TwoOperFun_f), pointer :: two_oper_fun   !context of callback subroutines
        type(QcMat), allocatable :: f_var_ao_dens(:)  !variable AO based density matrices
        type(QcMat), allocatable :: f_val_int(:)      !integral matrices
        integer(kind=4) ierr                          !error information
        integer(kind=QINT) imat                       !incremental recorder over matrices
        ! converts C pointer to Fortran QcMat type
        allocate(f_var_ao_dens(num_var_dens), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPTwoOperGetMat_f>> num_var_dens", num_var_dens
            stop "RSPTwoOperGetMat_f>> failed to allocate memory for f_var_ao_dens"
        end if
        ierr = QcMat_C_F_POINTER(A=f_var_ao_dens, c_A=var_ao_dens)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        allocate(f_val_int(num_int), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPTwoOperGetMat_f>> num_int", num_int
            stop "RSPTwoOperGetMat_f>> failed to allocate memory for f_val_int"
        end if
        ierr = QcMat_C_F_POINTER(A=f_val_int, c_A=val_int)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, two_oper_fun)
        ! invokes Fortran callback subroutine to calculate the integral matrices
        call two_oper_fun%get_two_oper_mat(num_pert,              &
                                           pert_labels,           &
                                           pert_orders,           &
                                           num_var_dens,          &
                                           f_var_ao_dens,         &
#if defined(OPENRSP_F_USER_CONTEXT)
                                           two_oper_fun%len_ctx,  &
                                           two_oper_fun%user_ctx, &
#endif
                                           num_int,               &
                                           f_val_int)
        ! cleans up
        nullify(two_oper_fun)
        ierr = QcMat_C_NULL_PTR(A=f_val_int)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        ierr = QcMat_C_NULL_PTR(A=f_var_ao_dens)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        deallocate(f_val_int)
        deallocate(f_var_ao_dens)
    end subroutine RSPTwoOperGetMat_f

    !% \brief calls Fortran callback subroutine to get expectation values of
    !      a two-electron operator
    !  \author Bin Gao
    !  \date 2014-08-06
    !  \param[integer]{in} num_pert number of perturbations
    !  \param[integer]{in} pert_labels labels of the perturbations
    !  \param[integer]{in} pert_orders orders of the perturbations
    !  \param[integer]{in} num_var_dens number of variable AO based density matrices
    !  \param[C_PTR:type]{in} var_ao_dens the variable AO based density matrices
    !  \param[integer]{in} num_contr_dens number of contracted AO based density matrices
    !  \param[C_PTR:type]{in} contr_ao_dens the contracted AO based density matrices
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !  \param[integer]{in} num_exp number of expectation values
    !% \param[real]{out} val_exp the expectation values
    subroutine RSPTwoOperGetExp_f(num_pert,       &
                                  pert_labels,    &
                                  pert_orders,    &
                                  num_var_dens,   &
                                  var_ao_dens,    &
                                  num_contr_dens, &
                                  contr_ao_dens,  &
                                  user_ctx,       &
                                  num_exp,        &
                                  val_exp)        &
        bind(C, name="RSPTwoOperGetExp_f")
        integer(kind=C_QINT), value, intent(in) :: num_pert
        integer(kind=C_QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
        integer(kind=C_QINT), value, intent(in) :: num_var_dens
        type(C_PTR), intent(in) :: var_ao_dens(num_var_dens)
        integer(kind=C_QINT), value, intent(in) :: num_contr_dens
        type(C_PTR), intent(in) :: contr_ao_dens(num_contr_dens)
        type(C_PTR), value, intent(in) :: user_ctx
        integer(kind=C_QINT), value, intent(in) :: num_exp
        real(kind=C_QREAL), intent(inout) :: val_exp(num_exp)
        type(TwoOperFun_f), pointer :: two_oper_fun     !context of callback subroutines
        type(QcMat), allocatable :: f_var_ao_dens(:)    !variable AO based density matrices
        type(QcMat), allocatable :: f_contr_ao_dens(:)  !contracted AO based density matrices
        integer(kind=4) ierr                            !error information
        integer(kind=QINT) imat                         !incremental recorder over matrices
        ! converts C pointer to Fortran QcMat type
        allocate(f_var_ao_dens(num_var_dens), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPTwoOperGetExp_f>> num_var_dens", num_var_dens
            stop "RSPTwoOperGetExp_f>> failed to allocate memory for f_var_ao_dens"
        end if
        ierr = QcMat_C_F_POINTER(A=f_var_ao_dens, c_A=var_ao_dens)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        allocate(f_contr_ao_dens(num_contr_dens), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPTwoOperGetExp_f>> num_contr_dens", num_contr_dens
            stop "RSPTwoOperGetExp_f>> failed to allocate memory for f_contr_ao_dens"
        end if
        ierr = QcMat_C_F_POINTER(A=f_contr_ao_dens, c_A=contr_ao_dens)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, two_oper_fun)
        ! invokes Fortran callback subroutine to calculate the expectation values
        call two_oper_fun%get_two_oper_exp(num_pert,              &
                                           pert_labels,           &
                                           pert_orders,           &
                                           num_var_dens,          &
                                           f_var_ao_dens,         &
                                           num_contr_dens,        &
                                           f_contr_ao_dens,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                                           two_oper_fun%len_ctx,  &
                                           two_oper_fun%user_ctx, &
#endif
                                           num_exp,               &
                                           val_exp)
        ! cleans up
        nullify(two_oper_fun)
        ierr = QcMat_C_NULL_PTR(A=f_contr_ao_dens)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        ierr = QcMat_C_NULL_PTR(A=f_var_ao_dens)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        deallocate(f_contr_ao_dens)
        deallocate(f_var_ao_dens)
        return
    end subroutine RSPTwoOperGetExp_f

    !% \brief cleans the context of callback subroutines of two-electron operator
    !  \author Bin Gao
    !  \date 2014-08-06
    !% \param[TwoOperFun_f:type]{inout} two_oper_fun the context of callback subroutines
    subroutine RSPTwoOperDestroy_f(two_oper_fun)
        type(TwoOperFun_f), intent(inout) :: two_oper_fun
#if defined(OPENRSP_F_USER_CONTEXT)
        two_oper_fun%len_ctx = 0
        deallocate(two_oper_fun%user_ctx)
#endif
        nullify(two_oper_fun%get_two_oper_mat)
        nullify(two_oper_fun%get_two_oper_exp)
    end subroutine RSPTwoOperDestroy_f

end module rsp_two_oper_f

#undef OPENRSP_API_SRC
