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

! basic data types
#include "api/qmatrix_c_type.h"

module openrsp_callback_f

    use, intrinsic :: iso_c_binding
    use qmatrix, only: QINT,  &
                       QREAL, &
                       QMat

    implicit none

    integer(kind=4), private, parameter :: STDOUT = 6

    ! type saving C struct's for calling C functions
    type, private :: RSP_CTX
        private
        type(C_PTR) :: rsp_solver = C_NULL_PTR
        type(C_PTR) :: nuc_contrib = C_NULL_PTR
        type(C_PTR) :: overlap = C_NULL_PTR
        type(C_PTR) :: one_oper = C_NULL_PTR
        type(C_PTR) :: two_oper = C_NULL_PTR
        type(C_PTR) :: xc_fun = C_NULL_PTR
    end type RSP_CTX

    type(RSP_CTX), save, private :: ct_saved

    public :: RSP_CTX_Create
    public :: RSP_CTX_Destroy
    public :: f_callback_RSPSolverGetSolution
    public :: f_callback_RSPNucContribGet
    public :: f_callback_RSPOverlapGetMat
    public :: f_callback_RSPOverlapGetExp
    public :: f_callback_RSPOneOperGetMat
    public :: f_callback_RSPOneOperGetExp
    public :: f_callback_RSPTwoOperGetMat
    public :: f_callback_RSPTwoOperGetExp
    public :: f_callback_RSPXCFunGetMat
    public :: f_callback_RSPXCFunGetExp

    interface
        integer(C_INT) function RSPOneOperGetMat(one_oper,      &
                                                 num_pert,      &
                                                 perturbations, &
                                                 pert_orders,   &
                                                 num_int,       &
                                                 val_int)       &
            bind(C, name="RSPOneOperGetMat")
            use, intrinsic :: iso_c_binding
            implicit none
            type(C_PTR), value, intent(in) :: one_oper
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: perturbations(num_pert)
            integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
            integer(kind=C_QINT), value, intent(in) :: num_int
            type(C_PTR), intent(in) :: val_int(num_int)
        end function RSPOneOperGetMat
    end interface

    contains

    ! creates the context for calling C functions
    subroutine RSP_CTX_Create(rsp_solver,      &
                              nuc_contrib,     &
                              overlap,         &
                              one_oper,        &
                              two_oper,        &
                              xc_fun)
        type(C_PTR), value, intent(in) :: rsp_solver
        type(C_PTR), value, intent(in) :: nuc_contrib
        type(C_PTR), value, intent(in) :: overlap
        type(C_PTR), value, intent(in) :: one_oper
        type(C_PTR), value, intent(in) :: two_oper
        type(C_PTR), value, intent(in) :: xc_fun
        ctx_saved%rsp_solver = rsp_solver
        ctx_saved%nuc_contrib = nuc_contrib
        ctx_saved%overlap = overlap
        ct_saved%one_oper = one_oper
        ctx_saved%two_oper = two_oper
        ctx_saved%xc_fun = xc_fun
    end subroutine RSP_CTX_Create

    ! cleans up the context for calling C functions
    subroutine RSP_CTX_Destroy()
        ctx_saved%rsp_solver = C_NULL_PTR
        ctx_saved%nuc_contrib = C_NULL_PTR
        ctx_saved%overlap = C_NULL_PTR
        ctx_saved%one_oper = C_NULL_PTR
        ctx_saved%two_oper = C_NULL_PTR
        ctx_saved%xc_fun = C_NULL_PTR
    end subroutine RSP_CTX_Destroy

    subroutine f_callback_RSPSolverGetSolution()
    end subroutine f_callback_RSPSolverGetSolution

    subroutine f_callback_RSPNucContribGet()
    end subroutine f_callback_RSPNucContribGet

    subroutine f_callback_RSPOverlapGetMat()
    end subroutine f_callback_RSPOverlapGetMat

    subroutine f_callback_RSPOverlapGetExp()
    end subroutine f_callback_RSPOverlapGetExp

    ! callback subroutine to get the integral matrices of one-electron operator
    subroutine f_callback_RSPOneOperGetMat(num_pert,      &
                                           perturbations, &
                                           pert_orders,   &
                                           num_int,       &
                                           val_int)
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: perturbations(num_pert)
        integer(kind=QINT), intent(in) :: pert_orders(num_pert)
        integer(kind=QINT), intent(in) :: num_int
        type(QMat), intent(inout) :: val_int(num_int)
        type(C_PTR), allocatable :: c_val_int(:)
        integer imat
        integer(kind=4) ierr
        allocate(c_val_int(num_int), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "f_callback_RSPOneOperGetMat>> num_int", num_int
            stop "f_callback_RSPOneOperGetMat>> failed to allocate memory for c_val_int"
        end if
        do imat = 1, num_int
            c_val_int(imat) = c_loc(val_int(imat))
        end do
        ierr = RSPOneOperGetMat(ct_saved%one_oper, &
                                num_pert,          &
                                perturbations,     &
                                pert_orders,       &
                                num_int,           &
                                c_val_int)
        deallocate(c_val_int)
    end subroutine f_callback_RSPOneOperGetMat

    subroutine f_callback_RSPOneOperGetExp()
    end subroutine f_callback_RSPOneOperGetExp

    subroutine f_callback_RSPTwoOperGetMat()
    end subroutine f_callback_RSPTwoOperGetMat

    subroutine f_callback_RSPTwoOperGetExp()
    end subroutine f_callback_RSPTwoOperGetExp

    subroutine f_callback_RSPXCFunGetMat()
    end subroutine f_callback_RSPXCFunGetMat

    subroutine f_callback_RSPXCFunGetExp()
    end subroutine f_callback_RSPXCFunGetExp

end module openrsp_callback_f
