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


!!  2015-06-23, Bin Gao
!!  * first version

! basic data types
#include "api/qcmatrix_c_type.h"

#define OPENRSP_API_SRC "src/fortran/RSPXCFun.F90"

module RSPXCFun_f

    use, intrinsic :: iso_c_binding
    use qcmatrix_f, only: QINT,QREAL,QcMat,QcMat_C_F_POINTER,QcMat_C_NULL_PTR
    use RSPPertBasicTypes_f, only: QcPertInt, &
                                   C_QCPERTINT

    implicit none

    integer(kind=4), private, parameter :: STDOUT = 6

    ! user specified callback subroutines
    abstract interface
        subroutine XCFunGetMat_f(xc_len_tuple,       &
                                 xc_pert_tuple,      &
                                 num_freq_configs,   &
                                 pert_freq_category, &
                                 dmat_num_tuple,     &
                                 dmat_idx_tuple,     &
                                 num_dmat,           &
                                 dens_mat,           &
#if defined(OPENRSP_F_USER_CONTEXT)
                                 !len_ctx,            &
                                 user_ctx,           &
#endif
                                 num_int,            &
                                 val_int)
#if defined(OPENRSP_F_USER_CONTEXT)
            use, intrinsic :: iso_c_binding
#endif
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
            !integer(kind=QINT), intent(in) :: len_ctx
            !character(len=1), intent(in) :: user_ctx(len_ctx)
            type(C_PTR), intent(in) :: user_ctx
#endif
            integer(kind=QINT), intent(in) :: num_int
            type(QcMat), intent(inout) :: val_int(num_int)
        end subroutine XCFunGetMat_f
        subroutine XCFunGetExp_f(xc_len_tuple,       &
                                 xc_pert_tuple,      &
                                 num_freq_configs,   &
                                 pert_freq_category, &
                                 dmat_num_tuple,     &
                                 dmat_idx_tuple,     &
                                 num_dmat,           &
                                 dens_mat,           &
#if defined(OPENRSP_F_USER_CONTEXT)
                                 !len_ctx,            &
                                 user_ctx,           &
#endif
                                 num_exp,            &
                                 val_exp)
#if defined(OPENRSP_F_USER_CONTEXT)
            use, intrinsic :: iso_c_binding
#endif
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
            !integer(kind=QINT), intent(in) :: len_ctx
            !character(len=1), intent(in) :: user_ctx(len_ctx)
            type(C_PTR), intent(in) :: user_ctx
#endif
            integer(kind=QINT), intent(in) :: num_exp
            real(kind=QREAL), intent(inout) :: val_exp(2*num_exp)
        end subroutine XCFunGetExp_f
    end interface

    ! context of callback subroutines of XC functional
    type, public :: XCFunFun_f
        private
#if defined(OPENRSP_F_USER_CONTEXT)
        ! user-defined callback function context
        !integer(kind=QINT) :: len_ctx = 0
        !character(len=1), allocatable :: user_ctx(:)
        type(C_PTR) :: user_ctx
#endif
        ! callback functions
        procedure(XCFunGetMat_f), nopass, pointer :: get_xc_fun_mat
        procedure(XCFunGetExp_f), nopass, pointer :: get_xc_fun_exp
    end type XCFunFun_f

    public :: RSPXCFunCreate_f
    public :: RSPXCFunGetMat_f
    public :: RSPXCFunGetExp_f
    public :: RSPXCFunDestroy_f

    contains

    !% \brief creates the context of callback subroutines of XC functional
    !  \author Bin Gao
    !  \date 2015-06-23
    !  \param[XCFunFun_f:type]{inout} xcfun_fun the context of callback subroutines
    !  \param[character]{in} user_ctx user-defined callback function context
    !  \param[subroutine]{in} get_xc_fun_mat user specified function for
    !      getting integral matrices
    !  \param[subroutine]{in} get_xc_fun_exp user specified function for
    !%     getting expectation values
    subroutine RSPXCFunCreate_f(xcfun_fun,      &
#if defined(OPENRSP_F_USER_CONTEXT)
                                user_ctx,       &
#endif
                                get_xc_fun_mat, &
                                get_xc_fun_exp)
        type(XCFunFun_f), intent(inout) :: xcfun_fun
#if defined(OPENRSP_F_USER_CONTEXT)
        !character(len=1), intent(in) :: user_ctx(:)
        type(C_PTR), intent(in) :: user_ctx
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
                                      !len_ctx,            &
                                      user_ctx,           &
#endif
                                      num_int,            &
                                      val_int)
#if defined(OPENRSP_F_USER_CONTEXT)
                use, intrinsic :: iso_c_binding
#endif
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
                !integer(kind=QINT), intent(in) :: len_ctx
                !character(len=1), intent(in) :: user_ctx(len_ctx)
                type(C_PTR), intent(in) :: user_ctx
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
                                      !len_ctx,            &
                                      user_ctx,           &
#endif
                                      num_exp,            &
                                      val_exp)
#if defined(OPENRSP_F_USER_CONTEXT)
                use, intrinsic :: iso_c_binding
#endif
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
                !integer(kind=QINT), intent(in) :: len_ctx
                !character(len=1), intent(in) :: user_ctx(len_ctx)
                type(C_PTR), intent(in) :: user_ctx
#endif
                integer(kind=QINT), intent(in) :: num_exp
                real(kind=QREAL), intent(inout) :: val_exp(2*num_exp)
            end subroutine get_xc_fun_exp
        end interface
#if defined(OPENRSP_F_USER_CONTEXT)
        !integer(kind=4) ierr  !error information
        !xcfun_fun%len_ctx = size(user_ctx)
        !allocate(xcfun_fun%user_ctx(xcfun_fun%len_ctx), stat=ierr)
        !if (ierr/=0) then
        !    write(STDOUT,"(A,I8)") "RSPXCFunCreate_f>> length", xcfun_fun%len_ctx
        !    stop "RSPXCFunCreate_f>> failed to allocate memory for user_ctx"
        !end if
        xcfun_fun%user_ctx = user_ctx
#endif
        xcfun_fun%get_xc_fun_mat => get_xc_fun_mat
        xcfun_fun%get_xc_fun_exp => get_xc_fun_exp
    end subroutine RSPXCFunCreate_f

    !% \brief calls Fortran callback subroutine to get integral matrices of
    !      an XC functional
    !  \author Bin Gao
    !  \date 2015-06-23
    !  \param[integer]{in} len_tuple length of perturbation tuple on the XC functional
    !  \param[integer]{in} pert_tuple perturbation tuple on the XC functional
    !  \param[integer]{in} num_freq_configs the number of different frequency
    !      configurations to be considered for the perturbation tuple
    !  \param[integer]{in} pert_freq_category category of perturbation frequencies
    !  \param[integer]{in} dmat_len_tuple the number of different perturbation
    !      tuples of the AO based density matrices passed
    !  \param[integer]{in} dmat_idx_tuple indices of the density matrix
    !      perturbation tuples passed (canonically ordered)
    !  \param[integer]{in} num_dmat number of collected AO based density matrices for
    !      the passed density matrix perturbation tuples and all frequency configurations
    !  \param[C_PTR:type]{in} dens_mat the collected AO based density matrices
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !  \param[integer]{in} num_int number of the integral matrices
    !% \param[C_PTR:type]{inout} val_int the integral matrices
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
        type(XCFunFun_f), pointer :: xcfun_fun     !context of callback subroutines
        type(QcMat), allocatable :: f_dens_mat(:)  !AO based density matrices
        type(QcMat), allocatable :: f_val_int(:)   !integral matrices
        integer(kind=4) ierr                       !error information
        ! converts C pointer to Fortran QcMat type
        allocate(f_dens_mat(num_dmat), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPXCFunGetMat_f>> num_dmat", num_dmat
            stop "RSPXCFunGetMat_f>> failed to allocate memory for f_dens_mat"
        end if
        ierr = QcMat_C_F_POINTER(A=f_dens_mat, c_A=dens_mat)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        allocate(f_val_int(num_int), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPXCFunGetMat_f>> num_int", num_int
            stop "RSPXCFunGetMat_f>> failed to allocate memory for f_val_int"
        end if
        ierr = QcMat_C_F_POINTER(A=f_val_int, c_A=val_int)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, xcfun_fun)
        ! invokes Fortran callback subroutine to calculate the integral matrices
        call xcfun_fun%get_xc_fun_mat(xc_len_tuple,       &
                                      xc_pert_tuple,      &
                                      num_freq_configs,   &
                                      pert_freq_category, &
                                      dmat_num_tuple,     &
                                      dmat_idx_tuple,     &
                                      num_dmat,           &
                                      f_dens_mat,         &
#if defined(OPENRSP_F_USER_CONTEXT)
                                      !xcfun_fun%len_ctx,  &
                                      xcfun_fun%user_ctx, &
#endif
                                      num_int,            &
                                      f_val_int)
        ! cleans up
        nullify(xcfun_fun)
        ierr = QcMat_C_NULL_PTR(A=f_val_int)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        ierr = QcMat_C_NULL_PTR(A=f_dens_mat)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        deallocate(f_val_int)
        deallocate(f_dens_mat)
    end subroutine RSPXCFunGetMat_f

    !% \brief calls Fortran callback subroutine to get expectation values of
    !      an XC functional
    !  \author Bin Gao
    !  \date 2015-06-23
    !  \param[integer]{in} len_tuple length of perturbation tuple on the XC functional
    !  \param[integer]{in} pert_tuple perturbation tuple on the XC functional
    !  \param[integer]{in} num_freq_configs the number of different frequency
    !      configurations to be considered for the perturbation tuple
    !  \param[integer]{in} pert_freq_category category of perturbation frequencies
    !  \param[integer]{in} dmat_len_tuple the number of different perturbation
    !      tuples of the AO based density matrices passed
    !  \param[integer]{in} dmat_idx_tuple indices of the density matrix
    !      perturbation tuples passed (canonically ordered)
    !  \param[integer]{in} num_dmat number of collected AO based density matrices for
    !      the passed density matrix perturbation tuples and all frequency configurations
    !  \param[C_PTR:type]{in} dens_mat the collected AO based density matrices
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !  \param[integer]{in} num_exp number of expectation values
    !% \param[real]{out} val_exp the expectation values
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
        type(XCFunFun_f), pointer :: xcfun_fun     !context of callback subroutines
        type(QcMat), allocatable :: f_dens_mat(:)  !AO based density matrices
        integer(kind=4) ierr                       !error information
        ! converts C pointer to Fortran QcMat type
        allocate(f_dens_mat(num_dmat), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPXCFunGetExp_f>> num_dmat", num_dmat
            stop "RSPXCFunGetExp_f>> failed to allocate memory for f_dens_mat"
        end if
        ierr = QcMat_C_F_POINTER(A=f_dens_mat, c_A=dens_mat)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, xcfun_fun)
        ! invokes Fortran callback subroutine to calculate the expectation values
        call xcfun_fun%get_xc_fun_exp(xc_len_tuple,       &
                                      xc_pert_tuple,      &
                                      num_freq_configs,   &
                                      pert_freq_category, &
                                      dmat_num_tuple,     &
                                      dmat_idx_tuple,     &
                                      num_dmat,           &
                                      f_dens_mat,         &
#if defined(OPENRSP_F_USER_CONTEXT)
                                      !xcfun_fun%len_ctx,  &
                                      xcfun_fun%user_ctx, &
#endif
                                      num_exp,            &
                                      val_exp)
        ! cleans up
        nullify(xcfun_fun)
        ierr = QcMat_C_NULL_PTR(A=f_dens_mat)
        call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_API_SRC)
        deallocate(f_dens_mat)
        return
    end subroutine RSPXCFunGetExp_f

    !% \brief cleans the context of callback subroutines of XC functional
    !  \author Bin Gao
    !  \date 2015-06-23
    !% \param[XCFunFun_f:type]{inout} xcfun_fun the context of callback subroutines
    subroutine RSPXCFunDestroy_f(xcfun_fun)
        type(XCFunFun_f), intent(inout) :: xcfun_fun
!#if defined(OPENRSP_F_USER_CONTEXT)
!        xcfun_fun%len_ctx = 0
!        deallocate(xcfun_fun%user_ctx)
!#endif
        nullify(xcfun_fun%get_xc_fun_mat)
        nullify(xcfun_fun%get_xc_fun_exp)
    end subroutine RSPXCFunDestroy_f

end module RSPXCFun_f

#undef OPENRSP_API_SRC

