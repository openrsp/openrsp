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
!!  This file calls C functions to get required quantities for the Fortran
!!  recursive routine of OpenRSP.
!!
!!  2014-12-10, Bin Gao
!!  * first version

#define OPENRSP_DEBUG

! basic data types
#include "api/qcmatrix_c_type.h"

#define OPENRSP_AO_DENS_CALLBACK "src/ao_dens/adapter/openrsp_callback_f.F90"

module openrsp_callback_f

    use, intrinsic :: iso_c_binding
    use qcmatrix_f, only: QINT,     &
                          QREAL,    &
                          QcMat,    &
                          QSUCCESS, &
                          QcMat_C_LOC

    implicit none

    integer(kind=4), private, parameter :: STDOUT = 6

    ! type saving C struct's for calling C functions
    type, private :: RSP_CTX
        private
        type(C_PTR) :: rsp_solver = C_NULL_PTR
        type(C_PTR) :: nuc_hamilton = C_NULL_PTR
        type(C_PTR) :: overlap = C_NULL_PTR
        type(C_PTR) :: one_oper = C_NULL_PTR
        type(C_PTR) :: two_oper = C_NULL_PTR
        type(C_PTR) :: xc_fun = C_NULL_PTR
    end type RSP_CTX

    type(RSP_CTX), save, private :: ctx_saved

    public :: RSP_CTX_Create
    public :: RSP_CTX_Destroy
    public :: f_callback_RSPSolverGetLinearRSPSolution
    public :: f_callback_RSPNucHamiltonGetContributions
    public :: f_callback_RSPOverlapGetMat
    public :: f_callback_RSPOverlapGetExp
    public :: f_callback_RSPOneOperGetMat
    public :: f_callback_RSPOneOperGetExp
    public :: f_callback_RSPTwoOperGetMat
    public :: f_callback_RSPTwoOperGetExp
    public :: f_callback_RSPXCFunGetMat
    public :: f_callback_RSPXCFunGetExp

    interface
        integer(C_INT) function RSPSolverGetLinearRSPSolution(rsp_solver,    &
                                                              num_freq_sums, &
                                                              freq_sums,     &
                                                              size_pert,     &
                                                              RHS_mat,       &
                                                              rsp_param)     &
            bind(C, name="RSPSolverGetLinearRSPSolution")
            use, intrinsic :: iso_c_binding
            implicit none
            type(C_PTR), value, intent(in) :: rsp_solver
            integer(kind=C_QINT), value, intent(in) :: num_freq_sums
            real(kind=C_QREAL), intent(in) :: freq_sums(num_freq_sums)
            integer(kind=C_QINT), value, intent(in) :: size_pert
            type(C_PTR), intent(in) :: RHS_mat(size_pert*num_freq_sums)
            type(C_PTR), intent(in) :: rsp_param(size_pert*num_freq_sums)
        end function RSPSolverGetLinearRSPSolution
        integer(C_INT) function RSPOverlapGetMat(overlap,         &
                                                 bra_num_pert,    &
                                                 bra_pert_labels, &
                                                 bra_pert_orders, &
                                                 ket_num_pert,    &
                                                 ket_pert_labels, &
                                                 ket_pert_orders, &
                                                 num_pert,        &
                                                 pert_labels,     &
                                                 pert_orders,     &
                                                 num_int,         &
                                                 val_int)         &
            bind(C, name="RSPOverlapGetMat")
            use, intrinsic :: iso_c_binding
            implicit none
            type(C_PTR), value, intent(in) :: overlap
            integer(kind=C_QINT), value, intent(in) :: bra_num_pert
            integer(kind=C_QINT), intent(in) :: bra_pert_labels(bra_num_pert)
            integer(kind=C_QINT), intent(in) :: bra_pert_orders(bra_num_pert)
            integer(kind=C_QINT), value, intent(in) :: ket_num_pert
            integer(kind=C_QINT), intent(in) :: ket_pert_labels(ket_num_pert)
            integer(kind=C_QINT), intent(in) :: ket_pert_orders(ket_num_pert)
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: pert_labels(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
            integer(kind=C_QINT), value, intent(in) :: num_int
            type(C_PTR), intent(in) :: val_int(num_int)
        end function RSPOverlapGetMat
        integer(C_INT) function RSPOverlapGetExp(overlap,         &
                                                 bra_num_pert,    &
                                                 bra_pert_labels, &
                                                 bra_pert_orders, &
                                                 ket_num_pert,    &
                                                 ket_pert_labels, &
                                                 ket_pert_orders, &
                                                 num_pert,        &
                                                 pert_labels,     &
                                                 pert_orders,     &
                                                 num_dens,        &
                                                 ao_dens,         &
                                                 num_exp,         &
                                                 val_exp)         &
            bind(C, name="RSPOverlapGetExp")
            use, intrinsic :: iso_c_binding
            implicit none
            type(C_PTR), value, intent(in) :: overlap
            integer(kind=C_QINT), value, intent(in) :: bra_num_pert
            integer(kind=C_QINT), intent(in) :: bra_pert_labels(bra_num_pert)
            integer(kind=C_QINT), intent(in) :: bra_pert_orders(bra_num_pert)
            integer(kind=C_QINT), value, intent(in) :: ket_num_pert
            integer(kind=C_QINT), intent(in) :: ket_pert_labels(ket_num_pert)
            integer(kind=C_QINT), intent(in) :: ket_pert_orders(ket_num_pert)
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: pert_labels(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
            integer(kind=C_QINT), value, intent(in) :: num_dens
            type(C_PTR), intent(in) :: ao_dens(num_dens)
            integer(kind=C_QINT), value, intent(in) :: num_exp
            real(kind=C_QREAL), intent(inout) :: val_exp(2*num_exp)
        end function RSPOverlapGetExp
        integer(C_INT) function RSPOneOperGetMat(one_oper,    &
                                                 num_pert,    &
                                                 pert_labels, &
                                                 pert_orders, &
                                                 num_int,     &
                                                 val_int)     &
            bind(C, name="RSPOneOperGetMat")
            use, intrinsic :: iso_c_binding
            implicit none
            type(C_PTR), value, intent(in) :: one_oper
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: pert_labels(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
            integer(kind=C_QINT), value, intent(in) :: num_int
            type(C_PTR), intent(in) :: val_int(num_int)
        end function RSPOneOperGetMat
        integer(C_INT) function RSPOneOperGetExp(one_oper,    &
                                                 num_pert,    &
                                                 pert_labels, &
                                                 pert_orders, &
                                                 num_dens,    &
                                                 ao_dens,     &
                                                 num_exp,     &
                                                 val_exp)     &
            bind(C, name="RSPOneOperGetExp")
            use, intrinsic :: iso_c_binding
            implicit none
            type(C_PTR), value, intent(in) :: one_oper
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: pert_labels(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
            integer(kind=C_QINT), value, intent(in) :: num_dens
            type(C_PTR), intent(in) :: ao_dens(num_dens)
            integer(kind=C_QINT), value, intent(in) :: num_exp
            real(kind=C_QREAL), intent(inout) :: val_exp(2*num_exp)
        end function RSPOneOperGetExp
        integer(C_INT) function RSPTwoOperGetMat(two_oper,     &
                                                 num_pert,     &
                                                 pert_labels,  &
                                                 pert_orders,  &
                                                 num_var_dens, &
                                                 var_ao_dens,  &
                                                 num_int,      &
                                                 val_int)      &
            bind(C, name="RSPTwoOperGetMat")
            use, intrinsic :: iso_c_binding
            implicit none
            type(C_PTR), value, intent(in) :: two_oper
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: pert_labels(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
            integer(kind=C_QINT), value, intent(in) :: num_var_dens
            type(C_PTR), intent(in) :: var_ao_dens(num_var_dens)
            integer(kind=C_QINT), value, intent(in) :: num_int
            type(C_PTR), intent(in) :: val_int(num_int)
        end function RSPTwoOperGetMat
        integer(C_INT) function RSPTwoOperGetExp(one_oper,       &
                                                 num_pert,       &
                                                 pert_labels,    &
                                                 pert_orders,    &
                                                 num_var_dens,   &
                                                 var_ao_dens,    &
                                                 num_contr_dens, &
                                                 contr_ao_dens,  &
                                                 num_exp,        &
                                                 val_exp)        &
            bind(C, name="RSPTwoOperGetExp")   
            use, intrinsic :: iso_c_binding    
            implicit none
            type(C_PTR), value, intent(in) :: one_oper
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: pert_labels(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
            integer(kind=C_QINT), value, intent(in) :: num_var_dens
            type(C_PTR), intent(in) :: var_ao_dens(num_var_dens)
            integer(kind=C_QINT), value, intent(in) :: num_contr_dens
            type(C_PTR), intent(in) :: contr_ao_dens(num_contr_dens)
            integer(kind=C_QINT), value, intent(in) :: num_exp
            real(kind=C_QREAL), intent(inout) :: val_exp(2*num_exp)
        end function RSPTwoOperGetExp
    end interface

    contains

    ! creates the context for calling C functions
    subroutine RSP_CTX_Create(rsp_solver,   &
                              nuc_hamilton, &
                              overlap,      &
                              one_oper,     &
                              two_oper,     &
                              xc_fun)
        type(C_PTR), value, intent(in) :: rsp_solver
        type(C_PTR), value, intent(in) :: nuc_hamilton
        type(C_PTR), value, intent(in) :: overlap
        type(C_PTR), value, intent(in) :: one_oper
        type(C_PTR), value, intent(in) :: two_oper
        type(C_PTR), value, intent(in) :: xc_fun
        ctx_saved%rsp_solver = rsp_solver
        ctx_saved%nuc_hamilton = nuc_hamilton
        ctx_saved%overlap = overlap
        ctx_saved%one_oper = one_oper
        ctx_saved%two_oper = two_oper
        ctx_saved%xc_fun = xc_fun
    end subroutine RSP_CTX_Create

    ! cleans up the context for calling C functions
    subroutine RSP_CTX_Destroy()
        ctx_saved%rsp_solver = C_NULL_PTR
        ctx_saved%nuc_hamilton = C_NULL_PTR
        ctx_saved%overlap = C_NULL_PTR
        ctx_saved%one_oper = C_NULL_PTR
        ctx_saved%two_oper = C_NULL_PTR
        ctx_saved%xc_fun = C_NULL_PTR
    end subroutine RSP_CTX_Destroy

    ! callback subroutine to get the solution of linear response equation
    subroutine f_callback_RSPSolverGetLinearRSPSolution(num_freq_sums, &
                                                        freq_sums,     &
                                                        size_pert,     &
                                                        RHS_mat,       &
                                                        rsp_param)
        integer(kind=QINT), intent(in) :: num_freq_sums
        real(kind=QREAL), intent(in) :: freq_sums(num_freq_sums)
        integer(kind=QINT), intent(in) :: size_pert
        type(QcMat), intent(in) :: RHS_mat(size_pert*num_freq_sums)
        type(QcMat), intent(inout) :: rsp_param(size_pert*num_freq_sums)
        type(C_PTR), allocatable :: c_RHS_mat(:)
        type(C_PTR), allocatable :: c_rsp_param(:)
        integer(kind=QINT) imat
        integer(kind=4) ierr
        if (c_associated(ctx_saved%rsp_solver)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", num_freq_sums, size_pert
#endif
            allocate(c_RHS_mat(size_pert*num_freq_sums), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_RHS_mat", &
                                  size_pert*num_freq_sums
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=RHS_mat, c_A=c_RHS_mat)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            allocate(c_rsp_param(size_pert*num_freq_sums), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_rsp_param", &
                                  size_pert*num_freq_sums
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=rsp_param, c_A=c_rsp_param)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            ierr = RSPSolverGetLinearRSPSolution(ctx_saved%rsp_solver, &
                                                 num_freq_sums,        &
                                                 freq_sums,            &
                                                 size_pert,            &
                                                 c_RHS_mat,            &
                                                 c_rsp_param)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            do imat = 1, size_pert*num_freq_sums
                c_RHS_mat(imat) = C_NULL_PTR
                c_rsp_param(imat) = C_NULL_PTR
            end do
            deallocate(c_RHS_mat)
            deallocate(c_rsp_param)
        else
            write(STDOUT,100) "null callback function for linear response equation solver"
            call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
        end if
100     format("f_callback_RSPSolverGetLinearRSPSolution>> ",A,2I12)
    end subroutine f_callback_RSPSolverGetLinearRSPSolution

    ! callback subroutine to get nuclear contributions
    subroutine f_callback_RSPNucHamiltonGetContributions(num_pert,     &
                                                         pert_labels,  &
                                                         pert_orders,  &
                                                         size_contrib, &
                                                         val_contrib)
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_orders(num_pert)
        integer(kind=QINT), intent(in) :: size_contrib
        real(kind=QREAL), intent(inout) :: val_contrib(size_contrib)
        if (c_associated(ctx_saved%nuc_hamilton)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", num_pert, size_contrib
#endif
            !write(STDOUT,100) "not implemented"
            !call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
        end if
100     format("f_callback_RSPNucHamiltonGetContributions>> ",A,2I12)
    end subroutine f_callback_RSPNucHamiltonGetContributions

    ! callback subroutine to get (perturbed) overlap integral matrices
    subroutine f_callback_RSPOverlapGetMat(bra_num_pert,    &
                                           bra_pert_labels, &
                                           bra_pert_orders, &
                                           ket_num_pert,    &
                                           ket_pert_labels, &
                                           ket_pert_orders, &
                                           num_pert,        &
                                           pert_labels,     &
                                           pert_orders,     &
                                           num_int,         &
                                           val_int)
        integer(kind=QINT), intent(in) :: bra_num_pert
        integer(kind=QINT), intent(in) :: bra_pert_labels(bra_num_pert)
        integer(kind=QINT), intent(in) :: bra_pert_orders(bra_num_pert)
        integer(kind=QINT), intent(in) :: ket_num_pert
        integer(kind=QINT), intent(in) :: ket_pert_labels(ket_num_pert)
        integer(kind=QINT), intent(in) :: ket_pert_orders(ket_num_pert)
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_orders(num_pert)
        integer(kind=QINT), intent(in) :: num_int
        type(QcMat), intent(inout) :: val_int(num_int)
        type(C_PTR), allocatable :: c_val_int(:)
        integer(kind=QINT) imat
        integer(kind=4) ierr
        if (c_associated(ctx_saved%overlap)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", bra_num_pert, ket_num_pert, num_pert, num_int
#endif
            allocate(c_val_int(num_int), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_val_int", num_int
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=val_int, c_A=c_val_int)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            ierr = RSPOverlapGetMat(ctx_saved%overlap, &
                                    bra_num_pert,      &
                                    bra_pert_labels,   &
                                    bra_pert_orders,   &
                                    ket_num_pert,      &
                                    ket_pert_labels,   &
                                    ket_pert_orders,   &
                                    num_pert,          &
                                    pert_labels,       &
                                    pert_orders,       &
                                    num_int,           &
                                    c_val_int)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            do imat = 1, num_int
                c_val_int(imat) = C_NULL_PTR
            end do
            deallocate(c_val_int)
        end if
100     format("f_callback_RSPOverlapGetMat>> ",A,4I12)
    end subroutine f_callback_RSPOverlapGetMat

    ! callback subroutine to get expectation values of (perturbed) overlap integrals
    subroutine f_callback_RSPOverlapGetExp(bra_num_pert,    &
                                           bra_pert_labels, &
                                           bra_pert_orders, &
                                           ket_num_pert,    &
                                           ket_pert_labels, &
                                           ket_pert_orders, &
                                           num_pert,        &
                                           pert_labels,     &
                                           pert_orders,     &
                                           num_dens,        &
                                           ao_dens,         &
                                           num_exp,         &
                                           val_exp)
        integer(kind=QINT), intent(in) :: bra_num_pert
        integer(kind=QINT), intent(in) :: bra_pert_labels(bra_num_pert)
        integer(kind=QINT), intent(in) :: bra_pert_orders(bra_num_pert)
        integer(kind=QINT), intent(in) :: ket_num_pert
        integer(kind=QINT), intent(in) :: ket_pert_labels(ket_num_pert)
        integer(kind=QINT), intent(in) :: ket_pert_orders(ket_num_pert)
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_orders(num_pert)
        integer(kind=QINT), intent(in) :: num_dens
        type(QcMat), intent(in) :: ao_dens(num_dens)
        integer(kind=QINT), intent(in) :: num_exp
        complex(kind=QREAL), intent(inout) :: val_exp(num_exp)
        type(C_PTR), allocatable :: c_ao_dens(:)
        real(kind=QREAL), allocatable :: c_val_exp(:)
        integer(kind=QINT) ival
        integer(kind=4) ierr
        if (c_associated(ctx_saved%overlap)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", bra_num_pert, ket_num_pert, num_pert, num_dens, num_exp
#endif
            allocate(c_ao_dens(num_dens), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_ao_dens", num_dens
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=ao_dens, c_A=c_ao_dens)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            allocate(c_val_exp(2*num_exp), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_val_exp", num_exp
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = RSPOverlapGetExp(ctx_saved%overlap, &
                                    bra_num_pert,      &
                                    bra_pert_labels,   &
                                    bra_pert_orders,   &
                                    ket_num_pert,      &
                                    ket_pert_labels,   &
                                    ket_pert_orders,   &
                                    num_pert,          &
                                    pert_labels,       &
                                    pert_orders,       &
                                    num_dens,          &
                                    c_ao_dens,         &
                                    num_exp,           &
                                    c_val_exp)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            do ival = 1, num_exp
                val_exp(ival) = cmplx(c_val_exp(2*ival-1), c_val_exp(2*ival), kind=QREAL)
            end do
            do ival = 1, num_dens
                c_ao_dens(ival) = C_NULL_PTR
            end do
            deallocate(c_ao_dens)
            deallocate(c_val_exp)
        end if
100     format("f_callback_RSPOverlapGetExp>> ",A,5I12)
    end subroutine f_callback_RSPOverlapGetExp

    ! callback subroutine to get (perturbed) one-electron integral matrices
    subroutine f_callback_RSPOneOperGetMat(num_pert,    &
                                           pert_labels, &
                                           pert_orders, &
                                           num_int,     &
                                           val_int)
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_orders(num_pert)
        integer(kind=QINT), intent(in) :: num_int
        type(QcMat), intent(inout) :: val_int(num_int)
        type(C_PTR), allocatable :: c_val_int(:)
        integer(kind=QINT) imat
        integer(kind=4) ierr
        if (c_associated(ctx_saved%one_oper)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", num_pert, num_int
#endif
            allocate(c_val_int(num_int), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_val_int", num_int
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=val_int, c_A=c_val_int)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            ierr = RSPOneOperGetMat(ctx_saved%one_oper, &
                                    num_pert,           &
                                    pert_labels,        &
                                    pert_orders,        &
                                    num_int,            &
                                    c_val_int)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            do imat = 1, num_int
                c_val_int(imat) = C_NULL_PTR
            end do
            deallocate(c_val_int)
        end if
100     format("f_callback_RSPOneOperGetMat>> ",A,2I12)
    end subroutine f_callback_RSPOneOperGetMat

    ! callback subroutine to get expectation values of (perturbed) one-electron integrals
    subroutine f_callback_RSPOneOperGetExp(num_pert,    &
                                           pert_labels, &
                                           pert_orders, &
                                           num_dens,    &
                                           ao_dens,     &
                                           num_exp,     &
                                           val_exp)
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_orders(num_pert)
        integer(kind=QINT), intent(in) :: num_dens
        type(QcMat), intent(in) :: ao_dens(num_dens)
        integer(kind=QINT), intent(in) :: num_exp
        complex(kind=QREAL), intent(inout) :: val_exp(num_exp)
        type(C_PTR), allocatable :: c_ao_dens(:)
        real(kind=QREAL), allocatable :: c_val_exp(:)
        integer(kind=QINT) ival
        integer(kind=4) ierr
        if (c_associated(ctx_saved%one_oper)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", num_pert, num_dens, num_exp
#endif
            allocate(c_ao_dens(num_dens), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_ao_dens", num_dens
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=ao_dens, c_A=c_ao_dens)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            allocate(c_val_exp(2*num_exp), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_val_exp", num_exp
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = RSPOneOperGetExp(ctx_saved%one_oper, &
                                    num_pert,           &
                                    pert_labels,        &
                                    pert_orders,        &
                                    num_dens,           &
                                    c_ao_dens,          &
                                    num_exp,            &
                                    c_val_exp)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            do ival = 1, num_exp
                val_exp(ival) = cmplx(c_val_exp(2*ival-1), c_val_exp(2*ival), kind=QREAL)
            end do
            do ival = 1, num_dens
                c_ao_dens(ival) = C_NULL_PTR
            end do
            deallocate(c_ao_dens)
            deallocate(c_val_exp)
        end if
100     format("f_callback_RSPOneOperGetExp>> ",A,3I12)
    end subroutine f_callback_RSPOneOperGetExp

    ! callback subroutine to get (perturbed) two-electron integral matrices
    subroutine f_callback_RSPTwoOperGetMat(num_pert,     &
                                           pert_labels,  &
                                           pert_orders,  &
                                           num_var_dens, &
                                           var_ao_dens,  &
                                           num_int,      &
                                           val_int)
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_orders(num_pert)
        integer(kind=QINT), intent(in) :: num_var_dens
        type(QcMat), intent(in) :: var_ao_dens(num_var_dens)
        integer(kind=QINT), intent(in) :: num_int
        type(QcMat), intent(inout) :: val_int(num_int)
        type(C_PTR), allocatable :: c_var_ao_dens(:)
        type(C_PTR), allocatable :: c_val_int(:)
        integer(kind=QINT) imat
        integer(kind=4) ierr
        if (c_associated(ctx_saved%two_oper)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", num_pert, num_var_dens, num_int
#endif
            allocate(c_var_ao_dens(num_var_dens), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_var_ao_dens", &
                                  num_var_dens
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=var_ao_dens, c_A=c_var_ao_dens)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            allocate(c_val_int(num_int), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_val_int", num_int
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=val_int, c_A=c_val_int)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            ierr = RSPTwoOperGetMat(ctx_saved%two_oper, &
                                    num_pert,           &
                                    pert_labels,        &
                                    pert_orders,        &
                                    num_var_dens,       &
                                    c_var_ao_dens,      &
                                    num_int,            &
                                    c_val_int)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            do imat = 1, num_var_dens
                c_var_ao_dens(imat) = C_NULL_PTR
            end do
            do imat = 1, num_int
                c_val_int(imat) = C_NULL_PTR
            end do
            deallocate(c_var_ao_dens)
            deallocate(c_val_int)
        end if
100     format("f_callback_RSPTwoOperGetMat>> ",A,3I12)
    end subroutine f_callback_RSPTwoOperGetMat

    ! callback subroutine to get expectation values of (perturbed) two-electron integrals
    subroutine f_callback_RSPTwoOperGetExp(num_pert,       &
                                           pert_labels,    &
                                           pert_orders,    &
                                           num_var_dens,   &
                                           var_ao_dens,    &
                                           num_contr_dens, &
                                           contr_ao_dens,  &
                                           num_exp,        &
                                           val_exp)
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_orders(num_pert)
        integer(kind=QINT), intent(in) :: num_var_dens
        type(QcMat), intent(in) :: var_ao_dens(num_var_dens)
        integer(kind=QINT), intent(in) :: num_contr_dens
        type(QcMat), intent(in) :: contr_ao_dens(num_contr_dens)
        integer(kind=QINT), intent(in) :: num_exp
        complex(kind=QREAL), intent(inout) :: val_exp(num_exp)
        type(C_PTR), allocatable :: c_var_ao_dens(:)
        type(C_PTR), allocatable :: c_contr_ao_dens(:)
        real(kind=QREAL), allocatable :: c_val_exp(:)
        integer(kind=QINT) ival
        integer(kind=4) ierr
        if (c_associated(ctx_saved%two_oper)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", num_pert, num_var_dens, num_contr_dens, num_exp
#endif
            allocate(c_var_ao_dens(num_var_dens), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_var_ao_dens", &
                                  num_var_dens
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=var_ao_dens, c_A=c_var_ao_dens)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            allocate(c_contr_ao_dens(num_contr_dens), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_contr_ao_dens", &
                                  num_contr_dens
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=contr_ao_dens, c_A=c_contr_ao_dens)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            allocate(c_val_exp(2*num_exp), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_val_exp", num_exp
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = RSPTwoOperGetExp(ctx_saved%two_oper, &
                                    num_pert,           &
                                    pert_labels,        &
                                    pert_orders,        &
                                    num_var_dens,       &
                                    c_var_ao_dens,      &
                                    num_contr_dens,     &
                                    c_contr_ao_dens,    &
                                    num_exp,            &
                                    c_val_exp)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            do ival = 1, num_exp
                val_exp(ival) = cmplx(c_val_exp(2*ival-1), c_val_exp(2*ival), kind=QREAL)
            end do
            do ival = 1, num_var_dens
                c_var_ao_dens(ival) = C_NULL_PTR
            end do
            do ival = 1, num_contr_dens
                c_contr_ao_dens(ival) = C_NULL_PTR
            end do
            deallocate(c_var_ao_dens)
            deallocate(c_contr_ao_dens)
            deallocate(c_val_exp)
        end if
100     format("f_callback_RSPTwoOperGetExp>> ",A,4I12)
    end subroutine f_callback_RSPTwoOperGetExp

    ! callback subroutine to get (perturbed) exchange-correlation functional matrices
    subroutine f_callback_RSPXCFunGetMat()
    end subroutine f_callback_RSPXCFunGetMat

    ! callback subroutine to get expectation values of (perturbed) exchange-correlation functional
    subroutine f_callback_RSPXCFunGetExp()
    end subroutine f_callback_RSPXCFunGetExp

end module openrsp_callback_f

#undef OPENRSP_AO_DENS_CALLBACK
