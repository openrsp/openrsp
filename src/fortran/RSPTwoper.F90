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


!!  2014-08-06, Bin Gao
!!  * first version

! basic data types
#include "api/qcmatrix_c_type.h"

#define OPENRSP_API_SRC "src/fortran/RSPTwoOper.F90"

module RSPTwoOper_f

    use, intrinsic :: iso_c_binding
    use qcmatrix_f, only: QINT,QREAL,QcMat,QcMat_C_F_POINTER,QcMat_C_NULL_PTR
    use RSPPertBasicTypes_f, only: QcPertInt, &
                                   C_QCPERTINT

    implicit none

    integer(kind=4), private, parameter :: STDOUT = 6

    ! user specified callback subroutines
    abstract interface
        subroutine TwoOperGetMat_f(oper_num_pert,    &
                                   oper_pert_labels, &
                                   oper_pert_orders, &
                                   num_dmat,         &
                                   dens_mat,         &
#if defined(OPENRSP_F_USER_CONTEXT)
                                   !len_ctx,          &
                                   user_ctx,         &
#endif
                                   num_int,          &
                                   val_int)
#if defined(OPENRSP_F_USER_CONTEXT)
            use, intrinsic :: iso_c_binding
#endif
            use qcmatrix_f, only: QINT,QREAL,QcMat
            use RSPPertBasicTypes_f, only: QcPertInt
            integer(kind=QINT), intent(in) :: oper_num_pert
            integer(kind=QcPertInt), intent(in) :: oper_pert_labels(oper_num_pert)
            integer(kind=QINT), intent(in) :: oper_pert_orders(oper_num_pert)
            integer(kind=QINT), intent(in) :: num_dmat
            type(QcMat), intent(in) :: dens_mat(num_dmat)
#if defined(OPENRSP_F_USER_CONTEXT)
            !integer(kind=QINT), intent(in) :: len_ctx
            !character(len=1), intent(in) :: user_ctx(len_ctx)
            type(C_PTR), intent(in) :: user_ctx
#endif
            integer(kind=QINT), intent(in) :: num_int
            type(QcMat), intent(inout) :: val_int(num_int)
        end subroutine TwoOperGetMat_f
        subroutine TwoOperGetExp_f(oper_num_pert,    &
                                   oper_pert_labels, &
                                   oper_pert_orders, &
                                   dmat_len_tuple,   &
                                   num_LHS_dmat,     &
                                   LHS_dens_mat,     &
                                   num_RHS_dmat,     &
                                   RHS_dens_mat,     &
#if defined(OPENRSP_F_USER_CONTEXT)
                                   !len_ctx,          &
                                   user_ctx,         &
#endif
                                   num_exp,          &
                                   val_exp)
#if defined(OPENRSP_F_USER_CONTEXT)
            use, intrinsic :: iso_c_binding
#endif
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
            !integer(kind=QINT), intent(in) :: len_ctx
            !character(len=1), intent(in) :: user_ctx(len_ctx)
            type(C_PTR), intent(in) :: user_ctx
#endif
            integer(kind=QINT), intent(in) :: num_exp
            real(kind=QREAL), intent(inout) :: val_exp(2*num_exp)
        end subroutine TwoOperGetExp_f
    end interface

    ! context of callback subroutines of two-electron operator
    type, public :: TwoOperFun_f
        private
#if defined(OPENRSP_F_USER_CONTEXT)
        ! user-defined callback function context
        !integer(kind=QINT) :: len_ctx = 0
        !character(len=1), allocatable :: user_ctx(:)
        type(C_PTR) :: user_ctx
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
        !character(len=1), intent(in) :: user_ctx(:)
        type(C_PTR), intent(in) :: user_ctx
#endif
        interface
            subroutine get_two_oper_mat(oper_num_pert,    &
                                        oper_pert_labels, &
                                        oper_pert_orders, &
                                        num_dmat,         &
                                        dens_mat,         &
#if defined(OPENRSP_F_USER_CONTEXT)
                                        !len_ctx,          &
                                        user_ctx,         &
#endif
                                        num_int,          &
                                        val_int)
#if defined(OPENRSP_F_USER_CONTEXT)
                use, intrinsic :: iso_c_binding
#endif
                use qcmatrix_f, only: QINT,QREAL,QcMat
                use RSPPertBasicTypes_f, only: QcPertInt
                integer(kind=QINT), intent(in) :: oper_num_pert
                integer(kind=QcPertInt), intent(in) :: oper_pert_labels(oper_num_pert)
                integer(kind=QINT), intent(in) :: oper_pert_orders(oper_num_pert)
                integer(kind=QINT), intent(in) :: num_dmat
                type(QcMat), intent(in) :: dens_mat(num_dmat)
#if defined(OPENRSP_F_USER_CONTEXT)
                !integer(kind=QINT), intent(in) :: len_ctx
                !character(len=1), intent(in) :: user_ctx(len_ctx)
                type(C_PTR), intent(in) :: user_ctx
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
                                        !len_ctx,          &
                                        user_ctx,         &
#endif
                                        num_exp,          &
                                        val_exp)
#if defined(OPENRSP_F_USER_CONTEXT)
                use, intrinsic :: iso_c_binding
#endif
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
                !integer(kind=QINT), intent(in) :: len_ctx
                !character(len=1), intent(in) :: user_ctx(len_ctx)
                type(C_PTR), intent(in) :: user_ctx
#endif
                integer(kind=QINT), intent(in) :: num_exp
                real(kind=QREAL), intent(inout) :: val_exp(2*num_exp)
            end subroutine get_two_oper_exp
        end interface
#if defined(OPENRSP_F_USER_CONTEXT)
        !integer(kind=4) ierr  !error information
        !two_oper_fun%len_ctx = size(user_ctx)
        !allocate(two_oper_fun%user_ctx(two_oper_fun%len_ctx), stat=ierr)
        !if (ierr/=0) then
        !    write(STDOUT,"(A,I8)") "RSPTwoOperCreate_f>> length", two_oper_fun%len_ctx
        !    stop "RSPTwoOperCreate_f>> failed to allocate memory for user_ctx"
        !end if
        two_oper_fun%user_ctx = user_ctx
#endif
        two_oper_fun%get_two_oper_mat => get_two_oper_mat
        two_oper_fun%get_two_oper_exp => get_two_oper_exp
    end subroutine RSPTwoOperCreate_f

    !% \brief calls Fortran callback subroutine to get integral matrices of
    !      a two-electron operator
    !  \author Bin Gao
    !  \date 2014-08-06
    !  \param[integer]{in} len_tuple length of perturbation tuple on the two-electron operator
    !  \param[integer]{in} pert_tuple perturbation tuple on the two-electron operator
    !  \param[integer]{in} num_dmat number of AO based density matrices
    !  \param[C_PTR:type]{in} dens_mat the AO based density matrices
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !  \param[integer]{in} num_int number of the integral matrices
    !% \param[C_PTR:type]{inout} val_int the integral matrices
    subroutine RSPTwoOperGetMat_f(oper_num_pert,    &
                                  oper_pert_labels, &
                                  oper_pert_orders, &
                                  num_dmat,         &
                                  dens_mat,         &
                                  user_ctx,         &
                                  num_int,          &
                                  val_int)          &
        bind(C, name="RSPTwoOperGetMat_f")
        integer(kind=C_QINT), value, intent(in) :: oper_num_pert
        integer(kind=C_QCPERTINT), intent(in) :: oper_pert_labels(oper_num_pert)
        integer(kind=C_QINT), intent(in) :: oper_pert_orders(oper_num_pert)
        integer(kind=C_QINT), value, intent(in) :: num_dmat
        type(C_PTR), intent(in) :: dens_mat(num_dmat)
        type(C_PTR), value, intent(in) :: user_ctx
        integer(kind=C_QINT), value, intent(in) :: num_int
        type(C_PTR), intent(inout) :: val_int(num_int)
        type(TwoOperFun_f), pointer :: two_oper_fun  !context of callback subroutines
        type(QcMat), allocatable :: f_dens_mat(:)    !AO based density matrices
        type(QcMat), allocatable :: f_val_int(:)     !integral matrices
        integer(kind=4) ierr                         !error information
        ! converts C pointer to Fortran QcMat type
        allocate(f_dens_mat(num_dmat), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPTwoOperGetMat_f>> num_dmat", num_dmat
            stop "RSPTwoOperGetMat_f>> failed to allocate memory for f_dens_mat"
        end if
        ierr = QcMat_C_F_POINTER(A=f_dens_mat, c_A=dens_mat)
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
        call two_oper_fun%get_two_oper_mat(oper_num_pert,         &
                                           oper_pert_labels,      &
                                           oper_pert_orders,      &
                                           num_dmat,              &
                                           f_dens_mat,            &
#if defined(OPENRSP_F_USER_CONTEXT)
                                           !two_oper_fun%len_ctx,  &
                                           two_oper_fun%user_ctx, &
#endif
                                           num_int,               &
                                           f_val_int)
        ! cleans up
        nullify(two_oper_fun)
        ierr = QcMat_C_NULL_PTR(A=f_val_int)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        ierr = QcMat_C_NULL_PTR(A=f_dens_mat)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        deallocate(f_val_int)
        deallocate(f_dens_mat)
    end subroutine RSPTwoOperGetMat_f

    !% \brief calls Fortran callback subroutine to get expectation values of
    !      a two-electron operator
    !  \author Bin Gao
    !  \date 2014-08-06
    !  \param[integer]{in} len_tuple length of perturbation tuple on the two-electron operator
    !  \param[integer]{in} pert_tuple perturbation tuple on the two-electron operator
    !  \param[integer]{in} dmat_len_tuple length of different perturbation tuples
    !      of the left-hand-side (LHS) and right-hand-side (RHS) AO based density
    !      matrices passed
    !  \param[integer]{in} num_LHS_dmat number of LHS AO based density matrices
    !      passed for each LHS density matrix perturbation tuple
    !  \param[C_PTR:type]{in} LHS_dens_mat the LHS AO based density matrices
    !  \param[integer]{in} num_RHS_dmat number of RHS AO based density matrices
    !      passed for each RHS density matrix perturbation tuple
    !  \param[C_PTR:type]{in} RHS_dens_mat the RHS AO based density matrices
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !  \param[integer]{in} num_exp number of expectation values
    !% \param[real]{out} val_exp the expectation values
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
        type(TwoOperFun_f), pointer :: two_oper_fun    !context of callback subroutines
        type(QcMat), allocatable :: f_LHS_dens_mat(:)  !LHS AO based density matrices
        type(QcMat), allocatable :: f_RHS_dens_mat(:)  !RHS AO based density matrices
        integer(kind=4) ierr                           !error information
        ! converts C pointer to Fortran QcMat type
        allocate(f_LHS_dens_mat(sum(num_LHS_dmat)), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPTwoOperGetExp_f>> sum(num_LHS_dmat)", &
                                   sum(num_LHS_dmat)
            stop "RSPTwoOperGetExp_f>> failed to allocate memory for f_LHS_dens_mat"
        end if
        ierr = QcMat_C_F_POINTER(A=f_LHS_dens_mat, c_A=LHS_dens_mat)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        allocate(f_RHS_dens_mat(sum(num_RHS_dmat)), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPTwoOperGetExp_f>> sum(num_RHS_dmat)", &
                                   sum(num_RHS_dmat)
            stop "RSPTwoOperGetExp_f>> failed to allocate memory for f_RHS_dens_mat"
        end if
        ierr = QcMat_C_F_POINTER(A=f_RHS_dens_mat, c_A=RHS_dens_mat)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, two_oper_fun)
        ! invokes Fortran callback subroutine to calculate the expectation values
        call two_oper_fun%get_two_oper_exp(oper_num_pert,         &
                                           oper_pert_labels,      &
                                           oper_pert_orders,      &
                                           dmat_len_tuple,        &
                                           num_LHS_dmat,          &
                                           f_LHS_dens_mat,        &
                                           num_RHS_dmat,          &
                                           f_RHS_dens_mat,        &
#if defined(OPENRSP_F_USER_CONTEXT)
                                           !two_oper_fun%len_ctx,  &
                                           two_oper_fun%user_ctx, &
#endif
                                           num_exp,               &
                                           val_exp)
        ! cleans up
        nullify(two_oper_fun)
        ierr = QcMat_C_NULL_PTR(A=f_RHS_dens_mat)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        ierr = QcMat_C_NULL_PTR(A=f_LHS_dens_mat)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        deallocate(f_RHS_dens_mat)
        deallocate(f_LHS_dens_mat)
        return
    end subroutine RSPTwoOperGetExp_f

    !% \brief cleans the context of callback subroutines of two-electron operator
    !  \author Bin Gao
    !  \date 2014-08-06
    !% \param[TwoOperFun_f:type]{inout} two_oper_fun the context of callback subroutines
    subroutine RSPTwoOperDestroy_f(two_oper_fun)
        type(TwoOperFun_f), intent(inout) :: two_oper_fun
!#if defined(OPENRSP_F_USER_CONTEXT)
!        two_oper_fun%len_ctx = 0
!        deallocate(two_oper_fun%user_ctx)
!#endif
        nullify(two_oper_fun%get_two_oper_mat)
        nullify(two_oper_fun%get_two_oper_exp)
    end subroutine RSPTwoOperDestroy_f

end module RSPTwoOper_f

#undef OPENRSP_API_SRC

