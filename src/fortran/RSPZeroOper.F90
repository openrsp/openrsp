!! OpenRSP: open-ended library for response theory
!! Copyright 2015 Radovan Bast,
!!                Daniel H. Friese,
!!                Bin Gao,
!!                Dan J. Jonsson,
!!                Magnus Ringholm,
!!                Kenneth Ruud,
!!                Andreas Thorvaldsen
!!
!!  This source code form is subject to the terms of the
!!  GNU Lesser General Public License, version 2.1.
!!  If a copy of the GNU LGPL v2.1 was not distributed with this
!!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.


!!  2015-06-23, Bin Gao
!!  * first version

! basic data types
#include "api/qcmatrix_c_type.h"

#define OPENRSP_API_SRC "src/fortran/RSPZeroOper.F90"

module RSPZeroOper_f

    use, intrinsic :: iso_c_binding
    use qcmatrix_f, only: QINT,QREAL
    use RSPPertBasicTypes_f, only: QcPertInt, &
                                   C_QCPERTINT

    implicit none

    integer(kind=4), private, parameter :: STDOUT = 6

    ! user specified callback subroutine
    abstract interface
        subroutine ZeroOperGetContrib_f(oper_num_pert,    &
                                        oper_pert_labels, &
                                        oper_pert_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                        !len_ctx,         &
                                        user_ctx,        &
#endif
                                        size_pert,       &
                                        val_oper)
#if defined(OPENRSP_F_USER_CONTEXT)
            use, intrinsic :: iso_c_binding
#endif
            use qcmatrix_f, only: QINT,QREAL
            use RSPPertBasicTypes_f, only: QcPertInt
            integer(kind=QINT), intent(in) :: oper_num_pert
            integer(kind=QcPertInt), intent(in) :: oper_pert_labels(oper_num_pert)
            integer(kind=QINT), intent(in) :: oper_pert_orders(oper_num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
            !integer(kind=QINT), intent(in) :: len_ctx
            !character(len=1), intent(in) :: user_ctx(len_ctx)
            type(C_PTR), intent(in) :: user_ctx
#endif
            integer(kind=QINT), intent(in) :: size_pert
            real(kind=QREAL), intent(inout) :: val_oper(2*size_pert)
        end subroutine ZeroOperGetContrib_f
    end interface

    ! context of callback subroutine of zero-electron operator
    type, public :: ZeroOperFun_f
        private
#if defined(OPENRSP_F_USER_CONTEXT)
        ! user-defined callback function context
        !integer(kind=QINT) :: len_ctx = 0
        !character(len=1), allocatable :: user_ctx(:)
        type(C_PTR) :: user_ctx
#endif
        ! callback function
        procedure(ZeroOperGetContrib_f), nopass, pointer :: get_zero_oper_contrib
    end type ZeroOperFun_f

    public :: RSPZeroOperCreate_f
    public :: RSPZeroOperGetContribution_f
    public :: RSPZeroOperDestroy_f

    contains

    !% \brief creates the context of callback subroutine of zero-electron operator
    !  \author Bin Gao
    !  \date 2015-06-23
    !  \param[ZeroOperFun_f:type]{inout} zero_oper_fun the context of callback subroutine
    !  \param[character]{in} user_ctx user-defined callback function context
    !  \param[subroutine]{in} get_zero_oper_contrib user specified function for
    !%     getting zero-electron operator contributions
    subroutine RSPZeroOperCreate_f(zero_oper_fun, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                   user_ctx,      &
#endif
                                   get_zero_oper_contrib)
        type(ZeroOperFun_f), intent(inout) :: zero_oper_fun
#if defined(OPENRSP_F_USER_CONTEXT)
        !character(len=1), intent(in) :: user_ctx(:)
        type(C_PTR), intent(in) :: user_ctx
#endif
        interface
            subroutine get_zero_oper_contrib(oper_num_pert,    &
                                             oper_pert_labels, &
                                             oper_pert_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                             !len_ctx,          &
                                             user_ctx,         &
#endif
                                             size_pert,        &
                                             val_oper)
#if defined(OPENRSP_F_USER_CONTEXT)
                use, intrinsic :: iso_c_binding
#endif
                use qcmatrix_f, only: QINT,QREAL
                use RSPPertBasicTypes_f, only: QcPertInt
                integer(kind=QINT), intent(in) :: oper_num_pert
                integer(kind=QcPertInt), intent(in) :: oper_pert_labels(oper_num_pert)
                integer(kind=QINT), intent(in) :: oper_pert_orders(oper_num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
                !integer(kind=QINT), intent(in) :: len_ctx
                !character(len=1), intent(in) :: user_ctx(len_ctx)
                type(C_PTR), intent(in) :: user_ctx
#endif
                integer(kind=QINT), intent(in) :: size_pert
                real(kind=QREAL), intent(inout) :: val_oper(2*size_pert)
            end subroutine get_zero_oper_contrib
        end interface
#if defined(OPENRSP_F_USER_CONTEXT)
        !integer(kind=4) ierr  !error information
        !zero_oper_fun%len_ctx = size(user_ctx)
        !allocate(zero_oper_fun%user_ctx(zero_oper_fun%len_ctx), stat=ierr)
        !if (ierr/=0) then
        !    write(STDOUT,"(A,I8)") "RSPZeroOperCreate_f>> length", zero_oper_fun%len_ctx
        !    stop "RSPZeroOperCreate_f>> failed to allocate memory for user_ctx"
        !end if
        zero_oper_fun%user_ctx = user_ctx
#endif
        zero_oper_fun%get_zero_oper_contrib => get_zero_oper_contrib
    end subroutine RSPZeroOperCreate_f

    !% \brief calls Fortran callback subroutine to get zero-electron operator contributions
    !  \author Bin Gao
    !  \date 2015-06-23
    !  \param[integer]{in} len_tuple length of perturbation tuple on the zero-electron operator
    !  \param[integer]{in} pert_tuple perturbation tuple on the zero-electron operator
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !  \param[integer]{in} size_pert size of the perturbations on the zero-electron operator
    !% \param[real]{out} val_oper the zero-electron operator contributions
    subroutine RSPZeroOperGetContribution_f(oper_num_pert,    &
                                            oper_pert_labels, &
                                            oper_pert_orders, &
                                            user_ctx,         &
                                            size_pert,        &
                                            val_oper)         &
        bind(C, name="RSPZeroOperGetContribution_f")
        integer(kind=C_QINT), value, intent(in) :: oper_num_pert
        integer(kind=C_QCPERTINT), intent(in) :: oper_pert_labels(oper_num_pert)
        integer(kind=C_QINT), intent(in) :: oper_pert_orders(oper_num_pert)
        type(C_PTR), value, intent(in) :: user_ctx
        integer(kind=C_QINT), value, intent(in) :: size_pert
        real(C_QREAL), intent(inout) :: val_oper(2*size_pert)
        type(ZeroOperFun_f), pointer ::  zero_oper_fun !context of callback subroutine
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, zero_oper_fun)
        ! invokes Fortran callback subroutine to calculate the zero-electron operator contributions
        call zero_oper_fun%get_zero_oper_contrib(oper_num_pert,          &
                                                 oper_pert_labels,       &
                                                 oper_pert_orders,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                                                 !zero_oper_fun%len_ctx,  &
                                                 zero_oper_fun%user_ctx, &
#endif
                                                 size_pert,              &
                                                 val_oper)
        ! cleans up
        nullify(zero_oper_fun)
        return
    end subroutine RSPZeroOperGetContribution_f

    !% \brief cleans the context of callback subroutine of zero-electron operator
    !  \author Bin Gao
    !  \date 2015-06-23
    !% \param[ZeroOperFun_f:type]{inout} zero_oper_fun the context of callback subroutine
    subroutine RSPZeroOperDestroy_f(zero_oper_fun)
        type(ZeroOperFun_f), intent(inout) :: zero_oper_fun
!#if defined(OPENRSP_F_USER_CONTEXT)
!        zero_oper_fun%len_ctx = 0
!        deallocate(zero_oper_fun%user_ctx)
!#endif
        nullify(zero_oper_fun%get_zero_oper_contrib)
    end subroutine RSPZeroOperDestroy_f

end module RSPZeroOper_f

#undef OPENRSP_API_SRC

