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
!!  This file calls C functions to get required quantities for the Fortran
!!  recursive routine of OpenRSP.
!!
!!  2014-12-10, Bin Gao
!!  * first version

#define OPENRSP_DEBUG

! basic data types
#include "api/qcmatrix_c_type.h"

#define OPENRSP_AO_DENS_CALLBACK "src/dens_mat/adapter/openrsp_callback_f.F90"

module openrsp_callback_f

    use, intrinsic :: iso_c_binding
    use qcmatrix_f, only: QINT,     &
                          QREAL,    &
                          QcMat,    &
                          QSUCCESS, &
                          QcMat_C_LOC
    use RSPPertBasicTypes_f, only: QcPertInt

    implicit none

    integer(kind=4), private, parameter :: STDOUT = 6

    ! Temporary solution for printing
    integer(kind=4), private, save :: IO_USER_OUTPUT = 6

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

    public :: f_callback_SetUserOutput
    public :: f_callback_UserOutput

    interface
        integer(C_INT) function RSPSolverGetLinearRSPSolution(rsp_solver,    &
                                                              num_pert,      &
                                                              num_comps,     &
                                                              num_freq_sums, &
                                                              freq_sums,     &
                                                              RHS_mat,       &
                                                              rsp_param)     &
            bind(C, name="RSPSolverGetLinearRSPSolution")
            use, intrinsic :: iso_c_binding
            implicit none
            type(C_PTR), value, intent(in) :: rsp_solver
            integer(kind=C_QINT), value, intent(in) :: num_pert
            integer(kind=C_QINT), intent(in) :: num_comps(num_pert)
            integer(kind=C_QINT), intent(in) :: num_freq_sums(num_pert)
            real(kind=C_QREAL), intent(in) :: freq_sums(2*sum(num_freq_sums))
            type(C_PTR), intent(in) :: RHS_mat(dot_product(num_comps,num_freq_sums))
            type(C_PTR), intent(in) :: rsp_param(dot_product(num_comps,num_freq_sums))
        end function RSPSolverGetLinearRSPSolution
        integer(C_INT) function RSPNucHamiltonGetContributions(nuc_hamilton,   &
                                                               nuc_len_tuple,  &
                                                               nuc_pert_tuple, &
                                                               size_pert,      &
                                                               val_nuc)        &
            bind(C, name="RSPNucHamiltonGetContributions")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            implicit none
            type(C_PTR), value, intent(in) :: nuc_hamilton
            integer(kind=C_QINT), value, intent(in) :: nuc_len_tuple
            integer(kind=C_QCPERTINT), intent(in) :: nuc_pert_tuple(nuc_len_tuple)
            integer(kind=C_QINT), value, intent(in) :: size_pert
            real(kind=C_QREAL), intent(inout) :: val_nuc(2*size_pert)
        end function RSPNucHamiltonGetContributions
        integer(C_INT) function RSPOverlapGetMat(overlap,         &
                                                 bra_len_tuple,   &
                                                 bra_pert_tuple,  &
                                                 ket_len_tuple,   &
                                                 ket_pert_tuple,  &
                                                 oper_len_tuple,  &
                                                 oper_pert_tuple, &
                                                 num_int,         &
                                                 val_int)         &
            bind(C, name="RSPOverlapGetMat")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            implicit none
            type(C_PTR), value, intent(in) :: overlap
            integer(kind=C_QINT), value, intent(in) :: bra_len_tuple
            integer(kind=C_QCPERTINT), intent(in) :: bra_pert_tuple(bra_len_tuple)
            integer(kind=C_QINT), value, intent(in) :: ket_len_tuple
            integer(kind=C_QCPERTINT), intent(in) :: ket_pert_tuple(ket_len_tuple)
            integer(kind=C_QINT), value, intent(in) :: oper_len_tuple
            integer(kind=C_QCPERTINT), intent(in) :: oper_pert_tuple(oper_len_tuple)
            integer(kind=C_QINT), value, intent(in) :: num_int
            type(C_PTR), intent(in) :: val_int(num_int)
        end function RSPOverlapGetMat
        integer(C_INT) function RSPOverlapGetExp(overlap,         &
                                                 bra_len_tuple,   &
                                                 bra_pert_tuple,  &
                                                 ket_len_tuple,   &
                                                 ket_pert_tuple,  &
                                                 oper_len_tuple,  &
                                                 oper_pert_tuple, &
                                                 num_dmat,        &
                                                 dens_mat,        &
                                                 num_exp,         &
                                                 val_exp)         &
            bind(C, name="RSPOverlapGetExp")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            implicit none
            type(C_PTR), value, intent(in) :: overlap
            integer(kind=C_QINT), value, intent(in) :: bra_len_tuple
            integer(kind=C_QCPERTINT), intent(in) :: bra_pert_tuple(bra_len_tuple)
            integer(kind=C_QINT), value, intent(in) :: ket_len_tuple
            integer(kind=C_QCPERTINT), intent(in) :: ket_pert_tuple(ket_len_tuple)
            integer(kind=C_QINT), value, intent(in) :: oper_len_tuple
            integer(kind=C_QCPERTINT), intent(in) :: oper_pert_tuple(oper_len_tuple)
            integer(kind=C_QINT), value, intent(in) :: num_dmat
            type(C_PTR), intent(in) :: dens_mat(num_dmat)
            integer(kind=C_QINT), value, intent(in) :: num_exp
            real(kind=C_QREAL), intent(inout) :: val_exp(2*num_exp)
        end function RSPOverlapGetExp
        integer(C_INT) function RSPOneOperGetMat(one_oper,        &
                                                 oper_len_tuple,  &
                                                 oper_pert_tuple, &
                                                 num_int,         &
                                                 val_int)         &
            bind(C, name="RSPOneOperGetMat")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            implicit none
            type(C_PTR), value, intent(in) :: one_oper
            integer(kind=C_QINT), value, intent(in) :: oper_len_tuple
            integer(kind=C_QCPERTINT), intent(in) :: oper_pert_tuple(oper_len_tuple)
            integer(kind=C_QINT), value, intent(in) :: num_int
            type(C_PTR), intent(in) :: val_int(num_int)
        end function RSPOneOperGetMat
        integer(C_INT) function RSPOneOperGetExp(one_oper,        &
                                                 oper_len_tuple,  &
                                                 oper_pert_tuple, &
                                                 num_dmat,        &
                                                 dens_mat,        &
                                                 num_exp,         &
                                                 val_exp)         &
            bind(C, name="RSPOneOperGetExp")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            implicit none
            type(C_PTR), value, intent(in) :: one_oper
            integer(kind=C_QINT), value, intent(in) :: oper_len_tuple
            integer(kind=C_QCPERTINT), intent(in) :: oper_pert_tuple(oper_len_tuple)
            integer(kind=C_QINT), value, intent(in) :: num_dmat
            type(C_PTR), intent(in) :: dens_mat(num_dmat)
            integer(kind=C_QINT), value, intent(in) :: num_exp
            real(kind=C_QREAL), intent(inout) :: val_exp(2*num_exp)
        end function RSPOneOperGetExp
        integer(C_INT) function RSPTwoOperGetMat(two_oper,        &
                                                 oper_len_tuple,  &
                                                 oper_pert_tuple, &
                                                 num_dmat,        &
                                                 dens_mat,        &
                                                 num_int,         &
                                                 val_int)         &
            bind(C, name="RSPTwoOperGetMat")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            implicit none
            type(C_PTR), value, intent(in) :: two_oper
            integer(kind=C_QINT), value, intent(in) :: oper_len_tuple
            integer(kind=C_QCPERTINT), intent(in) :: oper_pert_tuple(oper_len_tuple)
            integer(kind=C_QINT), value, intent(in) :: num_dmat
            type(C_PTR), intent(in) :: dens_mat(num_dmat)
            integer(kind=C_QINT), value, intent(in) :: num_int
            type(C_PTR), intent(in) :: val_int(num_int)
        end function RSPTwoOperGetMat
        integer(C_INT) function RSPTwoOperGetExp(one_oper,        &
                                                 oper_len_tuple,  &
                                                 oper_pert_tuple, &
                                                 dmat_len_tuple,  &
                                                 num_LHS_dmat,    &
                                                 LHS_dens_mat,    &
                                                 num_RHS_dmat,    &
                                                 RHS_dens_mat,    &
                                                 num_exp,         &
                                                 val_exp)         &
            bind(C, name="RSPTwoOperGetExp")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            implicit none
            type(C_PTR), value, intent(in) :: one_oper
            integer(kind=C_QINT), value, intent(in) :: oper_len_tuple
            integer(kind=C_QCPERTINT), intent(in) :: oper_pert_tuple(oper_len_tuple)
            integer(kind=C_QINT), value, intent(in) :: dmat_len_tuple
            integer(kind=C_QINT), intent(in) :: num_LHS_dmat(dmat_len_tuple)
            type(C_PTR), intent(in) :: LHS_dens_mat(sum(num_LHS_dmat))
            integer(kind=C_QINT), intent(in) :: num_RHS_dmat(dmat_len_tuple)
            type(C_PTR), intent(in) :: RHS_dens_mat(sum(num_RHS_dmat))
            integer(kind=C_QINT), value, intent(in) :: num_exp
            real(kind=C_QREAL), intent(inout) :: val_exp(2*num_exp)
        end function RSPTwoOperGetExp
        integer(C_INT) function RSPXCFunGetMat(xc_fun,             &
                                               xc_len_tuple,       &
                                               xc_pert_tuple,      &
                                               num_freq_configs,   &
                                               pert_freq_category, &
                                               dmat_num_tuple,     &
                                               dmat_idx_tuple,     &
                                               num_dmat,           &
                                               dens_mat,           &
                                               num_int,            &
                                               val_int)            &
            bind(C, name="RSPXCFunGetMat")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            implicit none
            type(C_PTR), value, intent(in) :: xc_fun
            integer(kind=C_QINT), value, intent(in) :: xc_len_tuple
            integer(kind=C_QCPERTINT), intent(in) :: xc_pert_tuple(xc_len_tuple)
            integer(kind=C_QINT), value, intent(in) :: num_freq_configs
            integer(kind=C_QINT), intent(in) :: &
                pert_freq_category(num_freq_configs*xc_len_tuple)
            integer(kind=C_QINT), value, intent(in) :: dmat_num_tuple
            integer(kind=C_QINT), intent(in) :: dmat_idx_tuple(dmat_num_tuple)
            integer(kind=C_QINT), value, intent(in) :: num_dmat
            type(C_PTR), intent(in) :: dens_mat(num_dmat)
            integer(kind=C_QINT), value, intent(in) :: num_int
            type(C_PTR), intent(in) :: val_int(num_int)
        end function RSPXCFunGetMat
        integer(C_INT) function RSPXCFunGetExp(xc_fun,             &
                                               xc_len_tuple,       &
                                               xc_pert_tuple,      &
                                               num_freq_configs,   &
                                               pert_freq_category, &
                                               dmat_num_tuple,     &
                                               dmat_idx_tuple,     &
                                               num_dmat,           &
                                               dens_mat,           &
                                               num_exp,            &
                                               val_exp)            &
            bind(C, name="RSPXCFunGetExp")
            use, intrinsic :: iso_c_binding
            use RSPPertBasicTypes_f, only: C_QCPERTINT
            implicit none
            type(C_PTR), value, intent(in) :: xc_fun
            integer(kind=C_QINT), value, intent(in) :: xc_len_tuple
            integer(kind=C_QCPERTINT), intent(in) :: xc_pert_tuple(xc_len_tuple)
            integer(kind=C_QINT), value, intent(in) :: num_freq_configs
            integer(kind=C_QINT), intent(in) :: &
                pert_freq_category(num_freq_configs*xc_len_tuple)
            integer(kind=C_QINT), value, intent(in) :: dmat_num_tuple
            integer(kind=C_QINT), intent(in) :: dmat_idx_tuple(dmat_num_tuple)
            integer(kind=C_QINT), value, intent(in) :: num_dmat
            type(C_PTR), intent(in) :: dens_mat(num_dmat)
            integer(kind=C_QINT), value, intent(in) :: num_exp
            real(kind=C_QREAL), intent(inout) :: val_exp(2*num_exp)
        end function RSPXCFunGetExp
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
    subroutine f_callback_RSPSolverGetLinearRSPSolution(num_pert,      &
                                                        num_comps,     &
                                                        num_freq_sums, &
                                                        freq_sums,     &
                                                        RHS_mat,       &
                                                        rsp_param)
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: num_comps(num_pert)
        integer(kind=QINT), intent(in) :: num_freq_sums(num_pert)
        complex(kind=C_QREAL), intent(in) :: freq_sums(sum(num_freq_sums))
        type(QcMat), intent(in) :: RHS_mat(dot_product(num_comps,num_freq_sums))
        type(QcMat), intent(inout) :: rsp_param(dot_product(num_comps,num_freq_sums))
        real(kind=C_QREAL), allocatable :: c_freq_sums(:)
        type(C_PTR), allocatable :: c_RHS_mat(:)
        type(C_PTR), allocatable :: c_rsp_param(:)
        integer(kind=QINT) imat
        integer(kind=4) ierr
        if (c_associated(ctx_saved%rsp_solver)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", num_pert, num_comps, num_freq_sums
#endif
            allocate(c_freq_sums(2*sum(num_freq_sums)), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_freq_sums", &
                                  2*sum(num_freq_sums)
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            do imat = 1, sum(num_freq_sums)
                c_freq_sums(2*imat-1) = real(freq_sums(imat))
                c_freq_sums(2*imat) = aimag(freq_sums(imat))
            end do
            allocate(c_RHS_mat(dot_product(num_comps,num_freq_sums)), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_RHS_mat", &
                                  dot_product(num_comps,num_freq_sums)
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=RHS_mat, c_A=c_RHS_mat)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            allocate(c_rsp_param(dot_product(num_comps,num_freq_sums)), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_rsp_param", &
                                  dot_product(num_comps,num_freq_sums)
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=rsp_param, c_A=c_rsp_param)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            ierr = RSPSolverGetLinearRSPSolution(ctx_saved%rsp_solver, &
                                                 num_pert,             &
                                                 num_comps,            &
                                                 num_freq_sums,        &
                                                 c_freq_sums,          &
                                                 c_RHS_mat,            &
                                                 c_rsp_param)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            deallocate(c_freq_sums)
            do imat = 1, dot_product(num_comps,num_freq_sums)
                c_RHS_mat(imat) = C_NULL_PTR
                c_rsp_param(imat) = C_NULL_PTR
            end do
            deallocate(c_RHS_mat)
            deallocate(c_rsp_param)
        else
            write(STDOUT,100) "null callback function for linear response equation solver"
            call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
        end if
100     format("f_callback_RSPSolverGetLinearRSPSolution>> ",A,600I6)
    end subroutine f_callback_RSPSolverGetLinearRSPSolution

    ! callback subroutine to get nuclear contributions
    subroutine f_callback_RSPNucHamiltonGetContributions(nuc_len_tuple,  &
                                                         nuc_pert_tuple, &
                                                         size_pert,      &
                                                         val_nuc)
        integer(kind=QINT), intent(in) :: nuc_len_tuple
        integer(kind=QcPertInt), intent(in) :: nuc_pert_tuple(nuc_len_tuple)
        integer(kind=QINT), intent(in) :: size_pert
        complex(kind=QREAL), intent(inout) :: val_nuc(size_pert)
        real(kind=QREAL), allocatable :: c_val_nuc(:)
        integer(kind=QINT) ival
        integer(kind=4) ierr
        if (c_associated(ctx_saved%nuc_hamilton)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", nuc_len_tuple, size_pert
#endif
            allocate(c_val_nuc(2*size_pert), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_val_nuc", &
                                  size_pert
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            c_val_nuc = 0.0
            ierr = RSPNucHamiltonGetContributions(ctx_saved%nuc_hamilton, &
                                                  nuc_len_tuple,          &
                                                  nuc_pert_tuple,         &
                                                  size_pert,              &
                                                  c_val_nuc)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            do ival = 1, size_pert
                val_nuc(ival) = val_nuc(ival) &
                              + cmplx(c_val_nuc(2*ival-1), c_val_nuc(2*ival), kind=QREAL)
            end do
            deallocate(c_val_nuc)
        end if
100     format("f_callback_RSPNucHamiltonGetContributions>> ",A,2I12)
    end subroutine f_callback_RSPNucHamiltonGetContributions

    ! callback subroutine to get (perturbed) overlap integral matrices
    subroutine f_callback_RSPOverlapGetMat(bra_len_tuple,   &
                                           bra_pert_tuple,  &
                                           ket_len_tuple,   &
                                           ket_pert_tuple,  &
                                           oper_len_tuple,  &
                                           oper_pert_tuple, &
                                           num_int,         &
                                           val_int)
        integer(kind=QINT), intent(in) :: bra_len_tuple
        integer(kind=QcPertInt), intent(in) :: bra_pert_tuple(bra_len_tuple)
        integer(kind=QINT), intent(in) :: ket_len_tuple
        integer(kind=QcPertInt), intent(in) :: ket_pert_tuple(ket_len_tuple)
        integer(kind=QINT), intent(in) :: oper_len_tuple
        integer(kind=QcPertInt), intent(in) :: oper_pert_tuple(oper_len_tuple)
        integer(kind=QINT), intent(in) :: num_int
        type(QcMat), intent(inout) :: val_int(num_int)
        type(C_PTR), allocatable :: c_val_int(:)
        integer(kind=QINT) imat
        integer(kind=4) ierr
        if (c_associated(ctx_saved%overlap)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", bra_len_tuple, ket_len_tuple, oper_len_tuple, num_int
#endif
            allocate(c_val_int(num_int), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_val_int", num_int
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=val_int, c_A=c_val_int)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            ierr = RSPOverlapGetMat(ctx_saved%overlap, &
                                    bra_len_tuple,     &
                                    bra_pert_tuple,    &
                                    ket_len_tuple,     &
                                    ket_pert_tuple,    &
                                    oper_len_tuple,    &
                                    oper_pert_tuple,   &
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
    subroutine f_callback_RSPOverlapGetExp(bra_len_tuple,   &
                                           bra_pert_tuple,  &
                                           ket_len_tuple,   &
                                           ket_pert_tuple,  &
                                           oper_len_tuple,  &
                                           oper_pert_tuple, &
                                           num_dmat,        &
                                           dens_mat,        &
                                           num_exp,         &
                                           val_exp)
        integer(kind=QINT), intent(in) :: bra_len_tuple
        integer(kind=QcPertInt), intent(in) :: bra_pert_tuple(bra_len_tuple)
        integer(kind=QINT), intent(in) :: ket_len_tuple
        integer(kind=QcPertInt), intent(in) :: ket_pert_tuple(ket_len_tuple)
        integer(kind=QINT), intent(in) :: oper_len_tuple
        integer(kind=QcPertInt), intent(in) :: oper_pert_tuple(oper_len_tuple)
        integer(kind=QINT), intent(in) :: num_dmat
        type(QcMat), intent(in) :: dens_mat(num_dmat)
        integer(kind=QINT), intent(in) :: num_exp
        complex(kind=QREAL), intent(inout) :: val_exp(num_exp)
        type(C_PTR), allocatable :: c_dens_mat(:)
        real(kind=QREAL), allocatable :: c_val_exp(:)
        integer(kind=QINT) ival
        integer(kind=4) ierr
        if (c_associated(ctx_saved%overlap)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", bra_len_tuple, ket_len_tuple, oper_len_tuple, num_dmat, num_exp
#endif
            allocate(c_dens_mat(num_dmat), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_dens_mat", num_dmat
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=dens_mat, c_A=c_dens_mat)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            allocate(c_val_exp(2*num_exp), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_val_exp", num_exp
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            c_val_exp = 0.0
            ierr = RSPOverlapGetExp(ctx_saved%overlap, &
                                    bra_len_tuple,     &
                                    bra_pert_tuple,    &
                                    ket_len_tuple,     &
                                    ket_pert_tuple,    &
                                    oper_len_tuple,    &
                                    oper_pert_tuple,   &
                                    num_dmat,          &
                                    c_dens_mat,        &
                                    num_exp,           &
                                    c_val_exp)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            do ival = 1, num_exp
                val_exp(ival) = val_exp(ival) &
                              + cmplx(c_val_exp(2*ival-1), c_val_exp(2*ival), kind=QREAL)
            end do
            do ival = 1, num_dmat
                c_dens_mat(ival) = C_NULL_PTR
            end do
            deallocate(c_dens_mat)
            deallocate(c_val_exp)
        end if
100     format("f_callback_RSPOverlapGetExp>> ",A,5I12)
    end subroutine f_callback_RSPOverlapGetExp

    ! callback subroutine to get (perturbed) one-electron integral matrices
    subroutine f_callback_RSPOneOperGetMat(oper_len_tuple,  &
                                           oper_pert_tuple, &
                                           num_int,         &
                                           val_int)
        integer(kind=QINT), intent(in) :: oper_len_tuple
        integer(kind=QcPertInt), intent(in) :: oper_pert_tuple(oper_len_tuple)
        integer(kind=QINT), intent(in) :: num_int
        type(QcMat), intent(inout) :: val_int(num_int)
        type(C_PTR), allocatable :: c_val_int(:)
        integer(kind=QINT) imat
        integer(kind=4) ierr
        if (c_associated(ctx_saved%one_oper)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", oper_len_tuple, num_int
#endif
            allocate(c_val_int(num_int), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_val_int", num_int
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=val_int, c_A=c_val_int)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            ierr = RSPOneOperGetMat(ctx_saved%one_oper, &
                                    oper_len_tuple,     &
                                    oper_pert_tuple,    &
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
    subroutine f_callback_RSPOneOperGetExp(oper_len_tuple,  &
                                           oper_pert_tuple, &
                                           num_dmat,        &
                                           dens_mat,        &
                                           num_exp,         &
                                           val_exp)
        integer(kind=QINT), intent(in) :: oper_len_tuple
        integer(kind=QcPertInt), intent(in) :: oper_pert_tuple(oper_len_tuple)
        integer(kind=QINT), intent(in) :: num_dmat
        type(QcMat), intent(in) :: dens_mat(num_dmat)
        integer(kind=QINT), intent(in) :: num_exp
        complex(kind=QREAL), intent(inout) :: val_exp(num_exp)
        type(C_PTR), allocatable :: c_dens_mat(:)
        real(kind=QREAL), allocatable :: c_val_exp(:)
        integer(kind=QINT) ival
        integer(kind=4) ierr
        if (c_associated(ctx_saved%one_oper)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", oper_len_tuple, num_dmat, num_exp
#endif
            allocate(c_dens_mat(num_dmat), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_dens_mat", num_dmat
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=dens_mat, c_A=c_dens_mat)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            allocate(c_val_exp(2*num_exp), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_val_exp", num_exp
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            c_val_exp = 0.0
            ierr = RSPOneOperGetExp(ctx_saved%one_oper, &
                                    oper_len_tuple,     &
                                    oper_pert_tuple,    &
                                    num_dmat,           &
                                    c_dens_mat,         &
                                    num_exp,            &
                                    c_val_exp)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            do ival = 1, num_exp
                val_exp(ival) = val_exp(ival) &
                              + cmplx(c_val_exp(2*ival-1), c_val_exp(2*ival), kind=QREAL)
            end do
            do ival = 1, num_dmat
                c_dens_mat(ival) = C_NULL_PTR
            end do
            deallocate(c_dens_mat)
            deallocate(c_val_exp)
        end if
100     format("f_callback_RSPOneOperGetExp>> ",A,3I12)
    end subroutine f_callback_RSPOneOperGetExp

    ! callback subroutine to get (perturbed) two-electron integral matrices
    subroutine f_callback_RSPTwoOperGetMat(oper_len_tuple,  &
                                           oper_pert_tuple, &
                                           num_dmat,        &
                                           dens_mat,        &
                                           num_int,         &
                                           val_int)
        integer(kind=QINT), intent(in) :: oper_len_tuple
        integer(kind=QcPertInt), intent(in) :: oper_pert_tuple(oper_len_tuple)
        integer(kind=QINT), intent(in) :: num_dmat
        type(QcMat), intent(in) :: dens_mat(num_dmat)
        integer(kind=QINT), intent(in) :: num_int
        type(QcMat), intent(inout) :: val_int(num_int)
        type(C_PTR), allocatable :: c_dens_mat(:)
        type(C_PTR), allocatable :: c_val_int(:)
        integer(kind=QINT) imat
        integer(kind=4) ierr
        if (c_associated(ctx_saved%two_oper)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", oper_len_tuple, num_dmat, num_int
#endif
            allocate(c_dens_mat(num_dmat), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_dens_mat", &
                                  num_dmat
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=dens_mat, c_A=c_dens_mat)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            allocate(c_val_int(num_int), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_val_int", num_int
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=val_int, c_A=c_val_int)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            ierr = RSPTwoOperGetMat(ctx_saved%two_oper, &
                                    oper_len_tuple,     &
                                    oper_pert_tuple,    &
                                    num_dmat,           &
                                    c_dens_mat,         &
                                    num_int,            &
                                    c_val_int)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            do imat = 1, num_dmat
                c_dens_mat(imat) = C_NULL_PTR
            end do
            do imat = 1, num_int
                c_val_int(imat) = C_NULL_PTR
            end do
            deallocate(c_dens_mat)
            deallocate(c_val_int)
        end if
100     format("f_callback_RSPTwoOperGetMat>> ",A,3I12)
    end subroutine f_callback_RSPTwoOperGetMat

    ! callback subroutine to get expectation values of (perturbed) two-electron integrals
    subroutine f_callback_RSPTwoOperGetExp(oper_len_tuple,  &
                                           oper_pert_tuple, &
                                           dmat_len_tuple,  &
                                           num_LHS_dmat,    &
                                           LHS_dens_mat,    &
                                           num_RHS_dmat,    &
                                           RHS_dens_mat,    &
                                           num_exp,         &
                                           val_exp)
        integer(kind=QINT), intent(in) :: oper_len_tuple
        integer(kind=QcPertInt), intent(in) :: oper_pert_tuple(oper_len_tuple)
        integer(kind=QINT), intent(in) :: dmat_len_tuple
        integer(kind=QINT), intent(in) :: num_LHS_dmat(dmat_len_tuple)
        type(QcMat), intent(in) :: LHS_dens_mat(sum(num_LHS_dmat))
        integer(kind=QINT), intent(in) :: num_RHS_dmat(dmat_len_tuple)
        type(QcMat), intent(in) :: RHS_dens_mat(sum(num_RHS_dmat))
        integer(kind=QINT), intent(in) :: num_exp
        complex(kind=QREAL), intent(inout) :: val_exp(num_exp)
        type(C_PTR), allocatable :: c_LHS_dens_mat(:)
        type(C_PTR), allocatable :: c_RHS_dens_mat(:)
        real(kind=QREAL), allocatable :: c_val_exp(:)
        integer(kind=QINT) ival
        integer(kind=4) ierr
        if (c_associated(ctx_saved%two_oper)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", oper_len_tuple, sum(num_LHS_dmat), sum(num_RHS_dmat), num_exp
#endif
            allocate(c_LHS_dens_mat(sum(num_LHS_dmat)), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_LHS_dens_mat", &
                                  sum(num_LHS_dmat)
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=LHS_dens_mat, c_A=c_LHS_dens_mat)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            allocate(c_RHS_dens_mat(sum(num_RHS_dmat)), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_RHS_dens_mat", &
                                  sum(num_RHS_dmat)
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=RHS_dens_mat, c_A=c_RHS_dens_mat)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            allocate(c_val_exp(2*num_exp), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_val_exp", num_exp
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            c_val_exp = 0.0
            ierr = RSPTwoOperGetExp(ctx_saved%two_oper, &
                                    oper_len_tuple,     &
                                    oper_pert_tuple,    &
                                    dmat_len_tuple,     &
                                    num_LHS_dmat,       &
                                    c_LHS_dens_mat,     &
                                    num_RHS_dmat,       &
                                    c_RHS_dens_mat,     &
                                    num_exp,            &
                                    c_val_exp)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            do ival = 1, num_exp
                val_exp(ival) = val_exp(ival) &
                              + cmplx(c_val_exp(2*ival-1), c_val_exp(2*ival), kind=QREAL)
            end do
            do ival = 1, sum(num_LHS_dmat)
                c_LHS_dens_mat(ival) = C_NULL_PTR
            end do
            do ival = 1, sum(num_RHS_dmat)
                c_RHS_dens_mat(ival) = C_NULL_PTR
            end do
            deallocate(c_LHS_dens_mat)
            deallocate(c_RHS_dens_mat)
            deallocate(c_val_exp)
        end if
100     format("f_callback_RSPTwoOperGetExp>> ",A,4I12)
    end subroutine f_callback_RSPTwoOperGetExp

    ! callback subroutine to get (perturbed) exchange-correlation functional matrices
    subroutine f_callback_RSPXCFunGetMat(xc_len_tuple,       &
                                         xc_pert_tuple,      &
                                         num_freq_configs,   &
                                         pert_freq_category, &
                                         dmat_num_tuple,     &
                                         dmat_idx_tuple,     &
                                         num_dmat,           &
                                         dens_mat,           &
                                         num_int,            &
                                         val_int)
        integer(kind=QINT), intent(in) :: xc_len_tuple
        integer(kind=QcPertInt), intent(in) :: xc_pert_tuple(xc_len_tuple)
        integer(kind=QINT), intent(in) :: num_freq_configs
        integer(kind=QINT), intent(in) :: &
            pert_freq_category(num_freq_configs*xc_len_tuple)
        integer(kind=QINT), intent(in) :: dmat_num_tuple
        integer(kind=QINT), intent(in) :: dmat_idx_tuple(dmat_num_tuple)
        integer(kind=QINT), intent(in) :: num_dmat
        type(QcMat), intent(in) :: dens_mat(num_dmat)
        integer(kind=QINT), intent(in) :: num_int
        type(QcMat), intent(inout) :: val_int(num_int)
        type(C_PTR), allocatable :: c_dens_mat(:)
        type(C_PTR), allocatable :: c_val_int(:)
        integer(kind=QINT) imat
        integer(kind=4) ierr
        if (c_associated(ctx_saved%xc_fun)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", xc_len_tuple, num_dmat, num_int
#endif
            allocate(c_dens_mat(num_dmat), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_dens_mat", &
                                  num_dmat
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=dens_mat, c_A=c_dens_mat)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            allocate(c_val_int(num_int), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_val_int", num_int
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=val_int, c_A=c_val_int)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            ierr = RSPXCFunGetMat(ctx_saved%xc_fun,   &
                                  xc_len_tuple,       &
                                  xc_pert_tuple,      &
                                  num_freq_configs,   &
                                  pert_freq_category, &
                                  dmat_num_tuple,     &
                                  dmat_idx_tuple,     &
                                  num_dmat,           &
                                  c_dens_mat,         &
                                  num_int,            &
                                  c_val_int)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            do imat = 1, num_dmat
                c_dens_mat(imat) = C_NULL_PTR
            end do
            do imat = 1, num_int
                c_val_int(imat) = C_NULL_PTR
            end do
            deallocate(c_dens_mat)
            deallocate(c_val_int)
        end if
100     format("f_callback_RSPXCFunGetMat>> ",A,3I12)
    end subroutine f_callback_RSPXCFunGetMat

    ! callback subroutine to get expectation values of (perturbed) exchange-correlation functional
    subroutine f_callback_RSPXCFunGetExp(xc_len_tuple,       &
                                         xc_pert_tuple,      &
                                         num_freq_configs,   &
                                         pert_freq_category, &
                                         dmat_num_tuple,     &
                                         dmat_idx_tuple,     &
                                         num_dmat,           &
                                         dens_mat,           &
                                         num_exp,            &
                                         val_exp)
        integer(kind=QINT), intent(in) :: xc_len_tuple
        integer(kind=QcPertInt), intent(in) :: xc_pert_tuple(xc_len_tuple)
        integer(kind=QINT), intent(in) :: num_freq_configs
        integer(kind=QINT), intent(in) :: &
            pert_freq_category(num_freq_configs*xc_len_tuple)
        integer(kind=QINT), intent(in) :: dmat_num_tuple
        integer(kind=QINT), intent(in) :: dmat_idx_tuple(dmat_num_tuple)
        integer(kind=QINT), intent(in) :: num_dmat
        type(QcMat), intent(in) :: dens_mat(num_dmat)
        integer(kind=QINT), intent(in) :: num_exp
        complex(kind=QREAL), intent(inout) :: val_exp(num_exp)
        type(C_PTR), allocatable :: c_dens_mat(:)
        real(kind=QREAL), allocatable :: c_val_exp(:)
        integer(kind=QINT) ival
        integer(kind=4) ierr
        if (c_associated(ctx_saved%xc_fun)) then
#if defined(OPENRSP_DEBUG)
            write(STDOUT,100) "size", xc_len_tuple, num_dmat, num_exp
#endif
            allocate(c_dens_mat(num_dmat), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_dens_mat", &
                                  num_dmat
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            ierr = QcMat_C_LOC(A=dens_mat, c_A=c_dens_mat)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            allocate(c_val_exp(2*num_exp), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,100) "failed to allocate memory for c_val_exp", num_exp
                call QErrorExit(STDOUT, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            end if
            c_val_exp = 0.0
            ierr = RSPXCFunGetExp(ctx_saved%xc_fun,   &
                                  xc_len_tuple,       &
                                  xc_pert_tuple,      &
                                  num_freq_configs,   &
                                  pert_freq_category, &
                                  dmat_num_tuple,     &
                                  dmat_idx_tuple,     &
                                  num_dmat,           &
                                  c_dens_mat,         &
                                  num_exp,            &
                                  c_val_exp)
            call QErrorCheckCode(STDOUT, ierr, __LINE__, OPENRSP_AO_DENS_CALLBACK)
            do ival = 1, num_exp
                val_exp(ival) = val_exp(ival) &
                              + cmplx(c_val_exp(2*ival-1), c_val_exp(2*ival), kind=QREAL)
            end do
            do ival = 1, num_dmat
                c_dens_mat(ival) = C_NULL_PTR
            end do
            deallocate(c_dens_mat)
            deallocate(c_val_exp)
        end if
100     format("f_callback_RSPXCFunGetExp>> ",A,3I12)
    end subroutine f_callback_RSPXCFunGetExp

    ! Temporary solution for printing
    subroutine f_callback_SetUserOutput(io_output)
        integer(kind=4), intent(in) :: io_output
        IO_USER_OUTPUT = io_output
    end subroutine f_callback_SetUserOutput

    ! Temporary solution for printing
    subroutine f_callback_UserOutput(out_str, out_level)
        character(*), intent(in) :: out_str
        integer(kind=4), intent(in) :: out_level
        integer(kind=4), parameter :: OUT_DEBUG = 2
        integer(kind=4), parameter :: OUT_ERROR = 0
        select case(out_level)
        case(OUT_ERROR)
            call lsquit(trim(out_str), IO_USER_OUTPUT)
        case default
            if (out_level<OUT_DEBUG) then
                write(IO_USER_OUTPUT, "(2A)") "[OUT]-> ", trim(out_str)
#if defined(OPENRSP_DEBUG)
            else
                write(IO_USER_OUTPUT, "(2A)") "[DEBUG]-> ", trim(out_str)
#endif
            end if
        end select
    end subroutine f_callback_UserOutput

end module openrsp_callback_f

#undef OPENRSP_AO_DENS_CALLBACK
