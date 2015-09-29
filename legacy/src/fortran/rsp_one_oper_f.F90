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
!!  This file implements subroutines of one-electron operators used in the Fortran APIs.
!!
!!  2014-08-02, Bin Gao
!!  * first version

! basic data types
#include "api/qcmatrix_c_type.h"

#define OPENRSP_API_SRC "src/f03/rsp_one_oper_f.F90"

module rsp_one_oper_f

    use, intrinsic :: iso_c_binding
    use qcmatrix_f, only: QINT,QREAL,QcMat,QcMat_C_F_POINTER,QcMat_C_NULL_PTR

    implicit none

    integer(kind=4), private, parameter :: STDOUT = 6

    ! user specified callback subroutines
    abstract interface
        subroutine OneOperGetMat_f(len_tuple,  &
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
        end subroutine OneOperGetMat_f
        subroutine OneOperGetExp_f(len_tuple,  &
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
        end subroutine OneOperGetExp_f
    end interface

    ! context of callback subroutines of one-electron operator
    type, public :: OneOperFun_f
        private
#if defined(OPENRSP_F_USER_CONTEXT)
        ! user-defined callback function context
        integer(kind=QINT) :: len_ctx = 0
        character(len=1), allocatable :: user_ctx(:)
#endif
        ! callback functions
        procedure(OneOperGetMat_f), nopass, pointer :: get_one_oper_mat
        procedure(OneOperGetExp_f), nopass, pointer :: get_one_oper_exp
    end type OneOperFun_f

    public :: RSPOneOperCreate_f
    public :: RSPOneOperGetMat_f
    public :: RSPOneOperGetExp_f
    public :: RSPOneOperDestroy_f

    contains

    !% \brief creates the context of callback subroutines of one-electron operator
    !  \author Bin Gao
    !  \date 2014-08-03
    !  \param[OneOperFun_f:type]{inout} one_oper_fun the context of callback subroutines
    !  \param[character]{in} user_ctx user-defined callback function context
    !  \param[subroutine]{in} get_one_oper_mat user specified function for
    !      getting integral matrices
    !  \param[subroutine]{in} get_one_oper_exp user specified function for
    !%     getting expectation values
    subroutine RSPOneOperCreate_f(one_oper_fun,     &
#if defined(OPENRSP_F_USER_CONTEXT)
                                  user_ctx,         &
#endif
                                  get_one_oper_mat, &
                                  get_one_oper_exp)
        type(OneOperFun_f), intent(inout) :: one_oper_fun
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
        end interface
#if defined(OPENRSP_F_USER_CONTEXT)
        integer(kind=4) ierr  !error information
        one_oper_fun%len_ctx = size(user_ctx)
        allocate(one_oper_fun%user_ctx(one_oper_fun%len_ctx), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPOneOperCreate_f>> length", one_oper_fun%len_ctx
            stop "RSPOneOperCreate_f>> failed to allocate memory for user_ctx"
        end if
        one_oper_fun%user_ctx = user_ctx
#endif
        one_oper_fun%get_one_oper_mat => get_one_oper_mat
        one_oper_fun%get_one_oper_exp => get_one_oper_exp
    end subroutine RSPOneOperCreate_f

    !% \brief calls Fortran callback subroutine to get integral matrices of
    !      a one-electron operator
    !  \author Bin Gao
    !  \date 2014-08-02
    !  \param[integer]{in} len_tuple length of perturbation tuple on the one-electron operator
    !  \param[integer]{in} pert_tuple perturbation tuple on the one-electron operator
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !  \param[integer]{in} num_int number of the integral matrices
    !% \param[C_PTR:type]{inout} val_int the integral matrices
    subroutine RSPOneOperGetMat_f(len_tuple,  &
                                  pert_tuple, &
                                  user_ctx,   &
                                  num_int,    &
                                  val_int)    &
        bind(C, name="RSPOneOperGetMat_f")
        integer(kind=C_QINT), value, intent(in) :: len_tuple
        integer(kind=C_QINT), intent(in) :: pert_tuple(len_tuple)
        type(C_PTR), value, intent(in) :: user_ctx
        integer(kind=C_QINT), value, intent(in) :: num_int
        type(C_PTR), intent(inout) :: val_int(num_int)
        type(OneOperFun_f), pointer :: one_oper_fun  !context of callback subroutines
        type(QcMat), allocatable :: f_val_int(:)     !integral matrices
        integer(kind=4) ierr                         !error information
        ! converts C pointer to Fortran QcMat type
        allocate(f_val_int(num_int), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPOneOperGetMat_f>> num_int", num_int
            stop "RSPOneOperGetMat_f>> failed to allocate memory for f_val_int"
        end if
        ierr = QcMat_C_F_POINTER(A=f_val_int, c_A=val_int)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, one_oper_fun)
        ! invokes Fortran callback subroutine to calculate the integral matrices
        call one_oper_fun%get_one_oper_mat(len_tuple,             &
                                           pert_tuple,            &
#if defined(OPENRSP_F_USER_CONTEXT)
                                           one_oper_fun%len_ctx,  &
                                           one_oper_fun%user_ctx, &
#endif
                                           num_int,               &
                                           f_val_int)
        ! cleans up
        nullify(one_oper_fun)
        ierr = QcMat_C_NULL_PTR(A=f_val_int)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        deallocate(f_val_int)
    end subroutine RSPOneOperGetMat_f

    !% \brief calls Fortran callback subroutine to get expectation values of
    !      a one-electron operator
    !  \author Bin Gao
    !  \date 2014-08-02
    !  \param[integer]{in} len_tuple length of perturbation tuple on the one-electron operator
    !  \param[integer]{in} pert_tuple perturbation tuple on the one-electron operator
    !  \param[integer]{in} num_dmat number of atomic orbital (AO) based density matrices
    !  \param[C_PTR:type]{inout} dens_mat the AO based density matrices
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !  \param[integer]{in} num_exp number of expectation values
    !% \param[real]{out} val_exp the expectation values
    subroutine RSPOneOperGetExp_f(len_tuple,  &
                                  pert_tuple, &
                                  num_dmat,   &
                                  dens_mat,   &
                                  user_ctx,   &
                                  num_exp,    &
                                  val_exp)    &
        bind(C, name="RSPOneOperGetExp_f")
        integer(kind=C_QINT), value, intent(in) :: len_tuple
        integer(kind=C_QINT), intent(in) :: pert_tuple(len_tuple)
        integer(kind=C_QINT), value, intent(in) :: num_dmat
        type(C_PTR), intent(in) :: dens_mat(num_dmat)
        type(C_PTR), value, intent(in) :: user_ctx
        integer(kind=C_QINT), value, intent(in) :: num_exp
        real(kind=C_QREAL), intent(inout) :: val_exp(num_exp)
        type(OneOperFun_f), pointer :: one_oper_fun  !context of callback subroutines
        type(QcMat), allocatable :: f_dens_mat(:)     !AO based density matrices
        integer(kind=4) ierr                         !error information
        ! converts C pointer to Fortran QcMat type
        allocate(f_dens_mat(num_dmat), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPOneOperGetExp_f>> num_dmat", num_dmat
            stop "RSPOneOperGetExp_f>> failed to allocate memory for f_dens_mat"
        end if
        ierr = QcMat_C_F_POINTER(A=f_dens_mat, c_A=dens_mat)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, one_oper_fun)
        ! invokes Fortran callback subroutine to calculate the expectation values
        call one_oper_fun%get_one_oper_exp(len_tuple,             &
                                           pert_tuple,            &
                                           num_dmat,              &
                                           f_dens_mat,            &
#if defined(OPENRSP_F_USER_CONTEXT)
                                           one_oper_fun%len_ctx,  &
                                           one_oper_fun%user_ctx, &
#endif
                                           num_exp,               &
                                           val_exp)
        ! cleans up
        nullify(one_oper_fun)
        ierr = QcMat_C_NULL_PTR(A=f_dens_mat)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        deallocate(f_dens_mat)
        return
    end subroutine RSPOneOperGetExp_f

    !% \brief cleans the context of callback subroutines of one-electron operator
    !  \author Bin Gao
    !  \date 2014-08-03
    !% \param[OneOperFun_f:type]{inout} one_oper_fun the context of callback subroutines
    subroutine RSPOneOperDestroy_f(one_oper_fun)
        type(OneOperFun_f), intent(inout) :: one_oper_fun
#if defined(OPENRSP_F_USER_CONTEXT)
        one_oper_fun%len_ctx = 0
        deallocate(one_oper_fun%user_ctx)
#endif
        nullify(one_oper_fun%get_one_oper_mat)
        nullify(one_oper_fun%get_one_oper_exp)
    end subroutine RSPOneOperDestroy_f

end module rsp_one_oper_f

#undef OPENRSP_API_SRC
