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
!!  This file implements the adapter between C APIs and Fortran recursive
!!  routine of OpenRSP.
!!
!!  2014-12-10, Bin Gao
!!  * first version

! data types between C/Fortran
#include "api/qcmatrix_c_type.h"

    subroutine OpenRSPGetRSPFun_f(num_atoms,        &
                                  num_props,        &
                                  len_tuple,        &
                                  pert_tuple,       &
                                  num_freq_configs, &
                                  pert_freqs,       &
                                  kn_rules,         &
                                  F_unpert,         &
                                  S_unpert,         &
                                  D_unpert,         &
                                  rsp_solver,       &
                                  zero_oper,        &
                                  overlap,          &
                                  one_oper,         &
                                  two_oper,         &
                                  xc_fun,           &
                                  size_rsp_funs,    &
                                  rsp_funs)         &
        bind(C, name="OpenRSPGetRSPFun_f")
        use, intrinsic :: iso_c_binding
        use qcmatrix_f!, only: QINT,              &
                      !        QREAL,             &
                      !        QcMat,             &
                      !        QSUCCESS,          &
                      !        QcMat_C_F_POINTER, &
                      !        QcMat_C_NULL_PTR,  &
                      !        QcMatWrite_f
        use openrsp_callback_f
        use RSPPertBasicTypes_f, only: QcPertInt, &
                                       C_QCPERTINT
        use rsp_pert_table
        use rsp_general, only: openrsp_get_property
        implicit none
        logical :: mem_calibrate
        integer :: max_mat, mem_result
        integer(kind=C_QINT), value, intent(in) :: num_atoms
        integer(kind=C_QINT), value, intent(in) :: num_props
        integer(kind=C_QINT), intent(in) :: len_tuple(num_props)
        integer(kind=C_QCPERTINT), intent(in) :: pert_tuple(sum(len_tuple))
        integer(kind=C_QINT), intent(in) :: num_freq_configs(num_props)
        real(kind=C_QREAL), intent(in) :: pert_freqs(2*dot_product(len_tuple,num_freq_configs))
        integer(kind=C_QINT), intent(in) :: kn_rules(num_props)
        type(C_PTR), value, intent(in) :: F_unpert
        type(C_PTR), value, intent(in) :: S_unpert
        type(C_PTR), value, intent(in) :: D_unpert
        type(C_PTR), value, intent(in) :: rsp_solver
        type(C_PTR), value, intent(in) :: zero_oper
        type(C_PTR), value, intent(in) :: overlap
        type(C_PTR), value, intent(in) :: one_oper
        type(C_PTR), value, intent(in) :: two_oper
        type(C_PTR), value, intent(in) :: xc_fun
        integer(kind=C_QINT), value, intent(in) :: size_rsp_funs
        real(kind=C_QREAL), intent(out) :: rsp_funs(2*size_rsp_funs)
        ! local variables for converting C arguments to Fortran ones
        integer(kind=QINT) num_coord
        integer(kind=QINT) num_all_pert
        integer(kind=QINT), allocatable :: f_pert_dims(:)
        integer(kind=QINT), allocatable :: f_pert_first_comp(:)
        character(4), allocatable :: f_pert_tuple(:)
        complex(kind=QREAL), allocatable :: f_pert_freqs(:)
        type(QcMat) f_F_unpert(1)
        type(QcMat) f_S_unpert(1)
        type(QcMat) f_D_unpert(1)
        complex(kind=QREAL), allocatable :: f_rsp_funs(:)
        integer(kind=QINT) ipert, jpert
        integer(kind=4) ierr

        ! gets the number of coordinates
        num_coord = 3*num_atoms
        ! gets the number of all perturbations
        num_all_pert = sum(len_tuple)
        ! gets the dimensions and labels of perturbations
        allocate(f_pert_dims(num_all_pert), stat=ierr)
        if (ierr/=0) then
            call f_callback_UserOutput("OpenRSPGetRSPFun_f>> failed to allocate memory for f_pert_dims", OUT_ERROR)
        end if
        allocate(f_pert_first_comp(num_all_pert), stat=ierr)
        if (ierr/=0) then
            call f_callback_UserOutput("OpenRSPGetRSPFun_f>> failed to allocate memory for f_pert_first_comp", OUT_ERROR)
        end if
        allocate(f_pert_tuple(num_all_pert), stat=ierr)
        if (ierr/=0) then
            call f_callback_UserOutput("OpenRSPGetRSPFun_f>> failed to allocate memory for f_pert_tuple", OUT_ERROR)
        end if
        do ipert = 1, num_all_pert
            select case (pert_tuple(ipert))
            case (RSP_GEO_PERT)
                f_pert_dims(ipert) = num_coord
            case (RSP_ELGR_PERT)
                f_pert_dims(ipert) = 6
            case default
                f_pert_dims(ipert) = 3
            end select
            f_pert_first_comp(ipert) = 1  !always starting from 1 for the time being
            f_pert_tuple(ipert) = CHAR_PERT_TABLE(pert_tuple(ipert))
        end do
        ! gets the frequencies of perturbations
        allocate(f_pert_freqs(size(pert_freqs)/2), stat=ierr)
        if (ierr/=0) then
            call f_callback_UserOutput("OpenRSPGetRSPFun_f>> failed to allocate memory for f_pert_freqs", OUT_ERROR)
        end if
        do ipert = 1, size(f_pert_freqs)
            f_pert_freqs(ipert) = cmplx(pert_freqs(2*ipert-1), pert_freqs(2*ipert), kind=QREAL)
        end do
        ! gets the matrices
        ierr = QcMat_C_F_POINTER(f_F_unpert, (/F_unpert/))
        if (ierr/=QSUCCESS) then
            call f_callback_UserOutput("OpenRSPGetRSPFun_f>> failed to call QcMat_C_F_POINTER(F)", OUT_ERROR)
        end if
        ierr = QcMat_C_F_POINTER(f_S_unpert, (/S_unpert/))
        if (ierr/=QSUCCESS) then
            call f_callback_UserOutput("OpenRSPGetRSPFun_f>> failed to call QcMat_C_F_POINTER(S)", OUT_ERROR)
        end if
        ierr = QcMat_C_F_POINTER(f_D_unpert, (/D_unpert/))
        if (ierr/=QSUCCESS) then
            call f_callback_UserOutput("OpenRSPGetRSPFun_f>> failed to call QcMat_C_F_POINTER(D)", OUT_ERROR)
        end if
        ! sets the context of callback functions
        call RSP_CTX_Create(rsp_solver, &
                            zero_oper,  &
                            overlap,    &
                            one_oper,   &
                            two_oper,   &
                            xc_fun)
        ! allocates memory for the results
        allocate(f_rsp_funs(size_rsp_funs), stat=ierr)
        if (ierr/=0) then
            call f_callback_UserOutput("OpenRSPGetRSPFun_f>> failed to allocate memory for f_rsp_funs", OUT_ERROR)
        end if
        f_rsp_funs = 0.0

        mem_calibrate = .FALSE.
        ! MaR: max_mat set to very high number to take matrix limitations out of use
        ! during development of other features
        max_mat = 999999999

        call openrsp_get_property(num_props,                                 &
                                  len_tuple,                                 &
                                  f_pert_dims,                               &
                                  f_pert_first_comp,                         &
                                  f_pert_tuple,                              &
                                  num_freq_configs,                          &
                                  f_pert_freqs,                              &
                                  kn_rules,                                  &
                                  f_F_unpert(1),                             &
                                  f_S_unpert(1),                             &
                                  f_D_unpert(1),                             &
                                  f_callback_RSPSolverGetLinearRSPSolution,  &
                                  f_callback_RSPZeroOperGetContribution,     &
                                  f_callback_RSPOverlapGetMat,               &
                                  f_callback_RSPOverlapGetExp,               &
                                  f_callback_RSPOneOperGetMat,               &
                                  f_callback_RSPOneOperGetExp,               &
                                  f_callback_RSPTwoOperGetMat,               &
                                  f_callback_RSPTwoOperGetExp,               &
                                  f_callback_RSPXCFunGetMat,                 &
                                  f_callback_RSPXCFunGetExp,                 &
                                  f_callback_UserOutput,                     &
                                  size(f_rsp_funs),                          &
                                  f_rsp_funs,                                &
                                  0,                                         &
                                  mem_calibrate=mem_calibrate,               &
                                  max_mat=max_mat,                           &
                                  mem_result=mem_result)

        ! assigns the results
        jpert = 0
        do ipert = 1, size_rsp_funs
            jpert = jpert+1
            rsp_funs(jpert) = real(f_rsp_funs(ipert))
            jpert = jpert+1
            rsp_funs(jpert) = aimag(f_rsp_funs(ipert))
        end do
        ! cleans up
        ierr = QcMat_C_NULL_PTR(A=f_F_unpert)
        if (ierr/=QSUCCESS) then
            call f_callback_UserOutput("OpenRSPGetRSPFun_f>> failed to call QcMat_C_NULL_PTR(F)", OUT_ERROR)
        end if
        ierr = QcMat_C_NULL_PTR(A=f_S_unpert)
        if (ierr/=QSUCCESS) then
            call f_callback_UserOutput("OpenRSPGetRSPFun_f>> failed to call QcMat_C_NULL_PTR(S)", OUT_ERROR)
        end if
        ierr = QcMat_C_NULL_PTR(A=f_D_unpert)
        if (ierr/=QSUCCESS) then
            call f_callback_UserOutput("OpenRSPGetRSPFun_f>> failed to call QcMat_C_NULL_PTR(D)", OUT_ERROR)
        end if
        deallocate(f_pert_dims)
        deallocate(f_pert_first_comp)
        deallocate(f_pert_tuple)
        deallocate(f_pert_freqs)
        call RSP_CTX_Destroy()
        deallocate(f_rsp_funs)
        return
100     format(A,50I8)
    end subroutine OpenRSPGetRSPFun_f
