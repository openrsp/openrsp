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
!!  This file implements the adapter between C APIs and Fortran recursive
!!  routine of OpenRSP.
!!
!!  2014-12-10, Bin Gao
!!  * first version

! data types between C/Fortran
#include "api/qmatrix_c_type.h"

    subroutine OpenRSPGetRSPFun_f(num_pert,        &
                                  perturbations,   &
                                  pert_orders,     &
                                  pert_freqs,      &
                                  kn,              &
                                  F_unpert,        &
                                  S_unpert,        &
                                  D_unpert,        &
                                  rsp_solver,      &
                                  nuc_contrib,     &
                                  overlap,         &
                                  one_oper,        &
                                  two_oper,        &
                                  xc_fun,          &
                                  !id_outp,         &
                                  property_size,   &
                                  rsp_tensor)      &
                                  !len_file_tensor, &
                                  !file_rsp_tensor) &
        bind(C, name="OpenRSPGetRSPFun_f")
        use, intrinsic :: iso_c_binding
        use qmatrix, only: QINT,QREAL,QMat,QSUCCESS
        use openrsp_callback_f
        use rsp_pert_table
        use rsp_general, only: openrsp_get_property_2014
        implicit none
        integer(kind=C_QINT), value, intent(in) :: num_pert
        integer(kind=C_QINT), intent(in) :: perturbations(num_pert)
        integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
        real(kind=C_QREAL), intent(in) :: pert_freqs(2*num_pert)
        integer(kind=C_QINT), intent(in) :: kn(2)
        type(C_PTR), intent(in) :: F_unpert
        type(C_PTR), intent(in) :: S_unpert
        type(C_PTR), intent(in) :: D_unpert
        type(C_PTR), value, intent(in) :: rsp_solver
        type(C_PTR), value, intent(in) :: nuc_contrib
        type(C_PTR), value, intent(in) :: overlap
        type(C_PTR), value, intent(in) :: one_oper
        type(C_PTR), value, intent(in) :: two_oper
        type(C_PTR), value, intent(in) :: xc_fun
        interface
            integer(C_INT) function RSPNucContribGetNumAtoms(nuc_contrib, &
                                                             num_atoms)   &
                bind(C, name="RSPNucContribGetNumAtoms")
                use, intrinsic :: iso_c_binding
                implicit none
                type(C_PTR), value, intent(in) :: nuc_contrib
                integer(kind=C_QINT), intent(out) :: num_atoms
            end function RSPNucContribGetNumAtoms
        end interface
        !integer, intent(in) :: id_outp
        integer(kind=C_QINT), value, intent(in) :: property_size
        real(kind=C_QREAL), intent(out) :: rsp_tensor(2*property_size)
        !integer(kind=C_QINT), value, intent(in) :: len_file_tensor
        !type(C_PTR), value, intent(in) :: file_rsp_tensor
        ! local variables for converting C arguments to Fortran ones
        integer(kind=4), parameter :: STDOUT = 6
        integer(kind=QINT) num_coord
        integer(kind=QINT) f_num_pert
        integer(kind=QINT), allocatable :: f_pert_dims(:)
        integer(kind=QINT), allocatable :: f_pert_first_comp(:)
        character(4), allocatable :: f_pert_labels(:)
        complex(kind=QREAL), allocatable :: f_pert_freqs(:)
        type(QMat), pointer :: f_F_unpert
        type(QMat), pointer :: f_S_unpert
        type(QMat), pointer :: f_D_unpert
        complex(kind=QREAL), allocatable :: f_rsp_tensor(:)
        !character(kind=C_CHAR), pointer :: ptr_file_tensor(:)
        !character, allocatable :: f_file_tensor(:)
        integer(kind=QINT) ipert, jpert, iorder
        integer(kind=4) ierr
        ! gets the number of coordinates
        ierr = RSPNucContribGetNumAtoms(nuc_contrib, num_coord)
        if (ierr/=QSUCCESS) then
            stop "OpenRSPGetRSPFun_f>> failed to call RSPNucContribGetNumAtoms"
        end if
        num_coord = 3*num_coord
        ! gets the perturbations
        f_num_pert = sum(pert_orders)
        allocate(f_pert_dims(f_num_pert), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "OpenRSPGetRSPFun_f>> f_num_pert", f_num_pert
            stop "OpenRSPGetRSPFun_f>> failed to allocate memory for f_pert_dims"
        end if
        allocate(f_pert_first_comp(f_num_pert), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "OpenRSPGetRSPFun_f>> f_num_pert", f_num_pert
            stop "OpenRSPGetRSPFun_f>> failed to allocate memory for f_pert_first_comp"
        end if
        allocate(f_pert_labels(f_num_pert), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "OpenRSPGetRSPFun_f>> f_num_pert", f_num_pert
            stop "OpenRSPGetRSPFun_f>> failed to allocate memory for f_pert_labels"
        end if
        allocate(f_pert_freqs(f_num_pert), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "OpenRSPGetRSPFun_f>> f_num_pert", f_num_pert
            stop "OpenRSPGetRSPFun_f>> failed to allocate memory for f_pert_freqs"
        end if
        jpert = 0
        do ipert = 1, num_pert
            do iorder = 1, pert_orders(ipert)
                jpert = jpert+1
                select case (perturbations(ipert))
                case (RSP_GEO_PERT)
                    f_pert_dims(jpert) = num_coord
                case (RSP_ELGR_PERT)
                    f_pert_dims(jpert) = 6
                case default
                    f_pert_dims(jpert) = 3
                end select
                f_pert_first_comp(jpert) = 1  !always starting from 1 for the time being
                f_pert_labels(jpert) = CHAR_PERT_TABLE(perturbations(ipert))
                f_pert_freqs(jpert) = cmplx(pert_freqs(2*ipert-1), pert_freqs(2*ipert))
            end do
        end do
        ! gets the matrices
        call c_f_pointer(F_unpert, f_F_unpert)
        call c_f_pointer(S_unpert, f_S_unpert)
        call c_f_pointer(D_unpert, f_D_unpert)
        ! sets the context of callback functions
        call RSP_CTX_Create(rsp_solver,  &
                            nuc_contrib, &
                            overlap,     &
                            one_oper,    &
                            two_oper,    &
                            xc_fun)
        ! allocates memory for the results
        allocate(f_rsp_tensor(property_size), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "OpenRSPGetRSPFun_f>> property_size", property_size
            stop "OpenRSPGetRSPFun_f>> failed to allocate memory for f_rsp_tensor"
        end if
        !! gets the file name of results
        !if (c_associated(file_rsp_tensor)) then
        !    call c_f_pointer(file_rsp_tensor, ptr_file_tensor, [len_file_tensor])
        !    allocate(f_file_tensor(len_file_tensor), stat=ierr)
        !    if (ierr/=0) then
        !        write(6,"(A,I8)") "OpenRSPGetRSPFun_f>> len_file_tensor", len_file_tensor
        !        stop "OpenRSPGetRSPFun_f>> failed to allocate memory for f_file_tensor"
        !    end if
        !    do ipert = 1_QINT, len_file_tensor
        !        f_file_tensor(ipert) = ptr_rsp_tensor(ipert)(1:1)
        !    end do
        !    ! gets the properties
        !    call openrsp_get_property_2014(f_num_pert,                      &
        !                                   f_pert_dims,                     &
        !                                   f_pert_first_comp,               &
        !                                   f_pert_labels,                   &
        !                                   f_pert_freqs,                    &
        !                                   kn,                              &
        !                                   f_F_unpert,                      &
        !                                   f_S_unpert,                      &
        !                                   f_D_unpert,                      &
        !                                   f_callback_RSPSolverGetSolution, &
        !                                   f_callback_RSPNucContribGet,     &
        !                                   f_callback_RSPOverlapGetMat,     &
        !                                   f_callback_RSPOverlapGetExp,     &
        !                                   f_callback_RSPOneOperGetMat,     &
        !                                   f_callback_RSPOneOperGetExp,     &
        !                                   f_callback_RSPTwoOperGetMat,     &
        !                                   f_callback_RSPTwoOperGetExp,     &
        !                                   f_callback_RSPXCFunGetMat,       &
        !                                   f_callback_RSPXCFunGetExp,       &
        !                                   STDOUT,                          &
        !                                   property_size,                   &
        !                                   f_rsp_tensor,                    &
        !                                   f_file_tensor)
        !    ! cleans up
        !    deallocate(f_file_tensor)
        !    nullify(ptr_file_tensor)
        !else
            call openrsp_get_property_2014(f_num_pert,                      &
                                           f_pert_dims,                     &
                                           f_pert_first_comp,               &
                                           f_pert_labels,                   &
                                           f_pert_freqs,                    &
                                           kn,                              &
                                           f_F_unpert,                      &
                                           f_S_unpert,                      &
                                           f_D_unpert,                      &
                                           f_callback_RSPSolverGetSolution, &
                                           f_callback_RSPNucContribGet,     &
                                           f_callback_RSPOverlapGetMat,     &
                                           f_callback_RSPOverlapGetExp,     &
                                           f_callback_RSPOneOperGetMat,     &
                                           f_callback_RSPOneOperGetExp,     &
                                           f_callback_RSPTwoOperGetMat,     &
                                           f_callback_RSPTwoOperGetExp,     &
                                           f_callback_RSPXCFunGetMat,       &
                                           f_callback_RSPXCFunGetExp,       &
                                           STDOUT,                          &
                                           property_size,                   &
                                           f_rsp_tensor)
        !end if
        ! assigns the results
        jpert = 0
        do ipert = 1, property_size
            jpert = jpert+1
            rsp_tensor(jpert) = real(f_rsp_tensor(ipert))
            jpert = jpert+1
            rsp_tensor(jpert) = aimag(f_rsp_tensor(ipert))
        end do
        ! cleans up
        deallocate(f_pert_dims)
        deallocate(f_pert_first_comp)
        deallocate(f_pert_labels)
        deallocate(f_pert_freqs)
        nullify(f_F_unpert)
        nullify(f_S_unpert)
        nullify(f_D_unpert)
        call RSP_CTX_Destroy()
        deallocate(f_rsp_tensor)
        return
    end subroutine OpenRSPGetRSPFun_f