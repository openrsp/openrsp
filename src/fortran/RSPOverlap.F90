!! OpenRSP: open-ended library for response theory
!! Copyright 2015 Radovan Bast,
!!                Daniel H. Friese,
!!                Bin Gao,
!!                Dan J. Jonsson,
!!                Magnus Ringholm,
!!                Kenneth Ruud
!!
!!  This source code form is subject to the terms of the
!!  GNU Lesser General Public License, version 2.1.
!!  If a copy of the GNU LGPL v2.1 was not distributed with this
!!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.


!!  2014-08-05, Bin Gao
!!  * first version

! basic data types
#include "api/qcmatrix_c_type.h"

#define OPENRSP_API_SRC "src/fortran/RSPOverlap.F90"

module RSPOverlap_f

    use, intrinsic :: iso_c_binding
    use qcmatrix_f, only: QINT,QREAL,QcMat,QcMat_C_F_POINTER,QcMat_C_NULL_PTR
    use RSPPertBasicTypes_f, only: QcPertInt, &
                                   C_QCPERTINT

    implicit none

    integer(kind=4), private, parameter :: STDOUT = 6

    ! user specified callback subroutines
    abstract interface
        subroutine OverlapGetMat_f(bra_num_pert,     &
                                   bra_pert_labels,  &
                                   bra_pert_orders,  &
                                   ket_num_pert,     &
                                   ket_pert_labels,  &
                                   ket_pert_orders,  &
                                   oper_num_pert,    &
                                   oper_pert_labels, &
                                   oper_pert_orders, &
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
            !integer(kind=QINT), intent(in) :: len_ctx
            !character(len=1), intent(in) :: user_ctx(len_ctx)
            type(C_PTR), intent(in) :: user_ctx
#endif
            integer(kind=QINT), intent(in) :: num_int
            type(QcMat), intent(inout) :: val_int(num_int)
        end subroutine OverlapGetMat_f
        subroutine OverlapGetExp_f(bra_num_pert,     &
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
            !integer(kind=QINT), intent(in) :: len_ctx
            !character(len=1), intent(in) :: user_ctx(len_ctx)
            type(C_PTR), intent(in) :: user_ctx
#endif
            integer(kind=QINT), intent(in) :: num_exp
            real(kind=QREAL), intent(inout) :: val_exp(2*num_exp)
        end subroutine OverlapGetExp_f
    end interface

    ! context of callback subroutines of overlap integrals
    type, public :: OverlapFun_f
        private
#if defined(OPENRSP_F_USER_CONTEXT)
        ! user-defined callback function context
        !integer(kind=QINT) :: len_ctx = 0
        !character(len=1), allocatable :: user_ctx(:)
        type(C_PTR) :: user_ctx
#endif
        ! callback functions
        procedure(OverlapGetMat_f), nopass, pointer :: get_overlap_mat
        procedure(OverlapGetExp_f), nopass, pointer :: get_overlap_exp
    end type OverlapFun_f

    public :: RSPOverlapCreate_f
    public :: RSPOverlapGetMat_f
    public :: RSPOverlapGetExp_f
    public :: RSPOverlapDestroy_f

    contains

    !% \brief creates the context of callback subroutines of overlap integrals
    !  \author Bin Gao
    !  \date 2014-08-05
    !  \param[OverlapFun_f:type]{inout} overlap_fun the context of callback subroutines
    !  \param[character]{in} user_ctx user-defined callback function context
    !  \param[subroutine]{in} get_overlap_mat user specified function for
    !      getting integral matrices
    !  \param[subroutine]{in} get_overlap_exp user specified function for
    !%     getting expectation values
    subroutine RSPOverlapCreate_f(overlap_fun,     &
#if defined(OPENRSP_F_USER_CONTEXT)
                                  user_ctx,        &
#endif
                                  get_overlap_mat, &
                                  get_overlap_exp)
        type(OverlapFun_f), intent(inout) :: overlap_fun
#if defined(OPENRSP_F_USER_CONTEXT)
        !character(len=1), intent(in) :: user_ctx(:)
        type(C_PTR), intent(in) :: user_ctx
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
                !integer(kind=QINT), intent(in) :: len_ctx
                !character(len=1), intent(in) :: user_ctx(len_ctx)
                type(C_PTR), intent(in) :: user_ctx
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
                !integer(kind=QINT), intent(in) :: len_ctx
                !character(len=1), intent(in) :: user_ctx(len_ctx)
                type(C_PTR), intent(in) :: user_ctx
#endif
                integer(kind=QINT), intent(in) :: num_exp
                real(kind=QREAL), intent(inout) :: val_exp(2*num_exp)
            end subroutine get_overlap_exp
        end interface
#if defined(OPENRSP_F_USER_CONTEXT)
        !integer(kind=4) ierr  !error information
        !overlap_fun%len_ctx = size(user_ctx)
        !allocate(overlap_fun%user_ctx(overlap_fun%len_ctx), stat=ierr)
        !if (ierr/=0) then
        !    write(STDOUT,"(A,I8)") "RSPOverlapCreate_f>> length", overlap_fun%len_ctx
        !    stop "RSPOverlapCreate_f>> failed to allocate memory for user_ctx"
        !end if
        overlap_fun%user_ctx = user_ctx
#endif
        overlap_fun%get_overlap_mat => get_overlap_mat
        overlap_fun%get_overlap_exp => get_overlap_exp
    end subroutine RSPOverlapCreate_f

    !% \brief calls Fortran callback subroutine to get integral matrices of overlap integrals
    !  \author Bin Gao
    !  \date 2014-08-05
    !  \param[integer]{in} bra_len_tuple length of the perturbation tuple on the bra
    !  \param[integer]{in} bra_pert_tuple perturbation tuple on the bra
    !  \param[integer]{in} ket_len_tuple length of the perturbation tuple on the ket
    !  \param[integer]{in} ket_pert_tuple perturbation tuple on the ket
    !  \param[integer]{in} len_tuple length of perturbation tuple on the overlap integrals
    !  \param[integer]{in} pert_tuple perturbation tuple on the overlap integrals
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !  \param[integer]{in} num_int number of the integral matrices
    !% \param[C_PTR:type]{inout} val_int the integral matrices
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
        type(OverlapFun_f), pointer :: overlap_fun  !context of callback subroutines
        type(QcMat), allocatable :: f_val_int(:)    !integral matrices
        integer(kind=4) ierr                        !error information
        ! converts C pointer to Fortran QcMat type
        allocate(f_val_int(num_int), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPOverlapGetMat_f>> num_int", num_int
            stop "RSPOverlapGetMat_f>> failed to allocate memory for f_val_int"
        end if
        ierr = QcMat_C_F_POINTER(A=f_val_int, c_A=val_int)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, overlap_fun)
        ! invokes Fortran callback subroutine to calculate the integral matrices
        call overlap_fun%get_overlap_mat(bra_num_pert,         &
                                         bra_pert_labels,      &
                                         bra_pert_orders,      &
                                         ket_num_pert,         &
                                         ket_pert_labels,      &
                                         ket_pert_orders,      &
                                         oper_num_pert,        &
                                         oper_pert_labels,     &
                                         oper_pert_orders,     &
#if defined(OPENRSP_F_USER_CONTEXT)
                                         !overlap_fun%len_ctx,  &
                                         overlap_fun%user_ctx, &
#endif
                                         num_int,              &
                                         f_val_int)
        ! cleans up
        nullify(overlap_fun)
        ierr = QcMat_C_NULL_PTR(A=f_val_int)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        deallocate(f_val_int)
    end subroutine RSPOverlapGetMat_f

    !% \brief calls Fortran callback subroutine to get expectation values of overlap integrals
    !  \author Bin Gao
    !  \date 2014-08-05
    !  \param[integer]{in} bra_len_tuple length of the perturbation tuple on the bra
    !  \param[integer]{in} bra_pert_tuple perturbation tuple on the bra
    !  \param[integer]{in} ket_len_tuple length of the perturbation tuple on the ket
    !  \param[integer]{in} ket_pert_tuple perturbation tuple on the ket
    !  \param[integer]{in} len_tuple length of perturbation tuple on the overlap integrals
    !  \param[integer]{in} pert_tuple perturbation tuple on the overlap integrals
    !  \param[integer]{in} num_dmat number of atomic orbital (AO) based density matrices
    !  \param[C_PTR:type]{inout} dens_mat the AO based density matrices
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !  \param[integer]{in} num_exp number of expectation values
    !% \param[real]{out} val_exp the expectation values
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
        real(C_QREAL), intent(inout) :: val_exp(2*num_exp)
        type(OverlapFun_f), pointer :: overlap_fun  !context of callback subroutines
        type(QcMat), allocatable :: f_dens_mat(:)    !AO based density matrices
        integer(kind=4) ierr                        !error information
        ! converts C pointer to Fortran QcMat type
        allocate(f_dens_mat(num_dmat), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPOverlapGetExp_f>> num_dmat", num_dmat
            stop "RSPOverlapGetExp_f>> failed to allocate memory for f_dens_mat"
        end if
        ierr = QcMat_C_F_POINTER(A=f_dens_mat, c_A=dens_mat)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, overlap_fun)
        ! invokes Fortran callback subroutine to calculate the expectation values
        call overlap_fun%get_overlap_exp(bra_num_pert,         &
                                         bra_pert_labels,      &
                                         bra_pert_orders,      &
                                         ket_num_pert,         &
                                         ket_pert_labels,      &
                                         ket_pert_orders,      &
                                         oper_num_pert,        &
                                         oper_pert_labels,     &
                                         oper_pert_orders,     &
                                         num_dmat,             &
                                         f_dens_mat,           &
#if defined(OPENRSP_F_USER_CONTEXT)
                                         !overlap_fun%len_ctx,  &
                                         overlap_fun%user_ctx, &
#endif
                                         num_exp,              &
                                         val_exp)
        ! cleans up
        nullify(overlap_fun)
        ierr = QcMat_C_NULL_PTR(A=f_dens_mat)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        deallocate(f_dens_mat)
        return
    end subroutine RSPOverlapGetExp_f

    !% \brief cleans the context of callback subroutines of overlap integrals
    !  \author Bin Gao
    !  \date 2014-08-05
    !% \param[OverlapFun_f:type]{inout} overlap_fun the context of callback subroutines
    subroutine RSPOverlapDestroy_f(overlap_fun)
        type(OverlapFun_f), intent(inout) :: overlap_fun
!#if defined(OPENRSP_F_USER_CONTEXT)
!        overlap_fun%len_ctx = 0
!        deallocate(overlap_fun%user_ctx)
!#endif
        nullify(overlap_fun%get_overlap_mat)
        nullify(overlap_fun%get_overlap_exp)
    end subroutine RSPOverlapDestroy_f

end module RSPOverlap_f

#undef OPENRSP_API_SRC

