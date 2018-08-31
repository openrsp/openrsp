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

    subroutine OpenRSPGetResidue_f(num_atoms,        &
                                   num_props,        &
                                   len_tuple,        &
                                   pert_tuple,       &
                                   residue_num_pert, &
                                   residue_idx_pert, &
                                   num_freq_configs, &
                                   pert_freqs,       &
                                   kn_rules,         &
                                   F_unpert,         &
                                   S_unpert,         &
                                   D_unpert,         &
                                   order_residue,    &
                                   num_excit,        &
                                   excit_energy,     &
                                   eigen_vector,     &
                                   rsp_solver,       &
                                   !zero_oper,        &
                                   overlap,          &
                                   one_oper,         &
                                   two_oper,         &
                                   xc_fun,           &
                                   r_flag,           &
                                   size_residues,    &
                                   residues)         &
        bind(C, name="OpenRSPGetResidue_f")
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
        integer(kind=C_QINT), value, intent(in) :: order_residue
        integer(kind=C_QINT), intent(in) :: residue_num_pert(order_residue*num_props)
        integer(kind=C_QINT), intent(in) :: residue_idx_pert(sum(residue_num_pert))
        integer(kind=C_QINT), intent(in) :: num_freq_configs(num_props)
        real(kind=C_QREAL), intent(in) :: &
            pert_freqs(2*dot_product(len_tuple,num_freq_configs)*num_excit)
        integer(kind=C_QINT), intent(in) :: kn_rules(num_props)
        type(C_PTR), value, intent(in) :: F_unpert
        type(C_PTR), value, intent(in) :: S_unpert
        type(C_PTR), value, intent(in) :: D_unpert
        integer(kind=C_QINT), value, intent(in) :: num_excit
        real(kind=C_QREAL), intent(in) :: excit_energy(order_residue*num_excit)
        type(C_PTR), intent(in) :: eigen_vector(order_residue*num_excit)
        type(C_PTR), value, intent(in) :: rsp_solver
        !type(C_PTR), value, intent(in) :: zero_oper
        type(C_PTR), value, intent(in) :: overlap
        type(C_PTR), value, intent(in) :: one_oper
        type(C_PTR), value, intent(in) :: two_oper
        type(C_PTR), value, intent(in) :: xc_fun
        integer(kind=C_QINT) r_flag
        integer(kind=C_QINT), value, intent(in) :: size_residues
        real(kind=C_QREAL), intent(out) :: residues(2*size_residues)
        ! local variables for converting C arguments to Fortran ones
        integer(kind=QINT) num_coord
        integer(kind=QINT) num_all_pert
        integer(kind=QINT), allocatable :: f_pert_dims(:)
        integer(kind=QINT), allocatable :: f_pert_first_comp(:)
        character(4), allocatable :: f_pert_tuple(:)
        integer(kind=QINT), allocatable :: residue_spec_index(:,:)
        integer(kind=QINT) offset_residue  !offset of residues per excited state
        complex(kind=QREAL) exenerg(2)
        complex(kind=QREAL), allocatable :: f_pert_freqs(:)
        type(QcMat) f_F_unpert(1)
        type(QcMat) f_S_unpert(1)
        type(QcMat) f_D_unpert(1)
        type(QcMat) X_unpert(2)
        integer(kind=QINT) resize_per_excit  !size of residues per excited state
        complex(kind=QREAL), allocatable :: f_residues(:)
        !character(kind=C_CHAR), pointer :: ptr_file_tensor(:)
        !character, allocatable :: f_file_tensor(:)
        integer(kind=QINT) ipert, jpert
        integer(kind=QINT) iext
        integer(kind=4) ierr

        ! gets the number of coordinates
        num_coord = 3*num_atoms
        ! gets the number of all perturbations
        num_all_pert = sum(len_tuple)
        ! gets the dimensions and labels of perturbations
        allocate(f_pert_dims(num_all_pert), stat=ierr)
        if (ierr/=0) then
            call f_callback_UserOutput("OpenRSPGetResidue_f>> failed to allocate memory for f_pert_dims", OUT_ERROR)
        end if
        allocate(f_pert_first_comp(num_all_pert), stat=ierr)
        if (ierr/=0) then
            call f_callback_UserOutput("OpenRSPGetResidue_f>> failed to allocate memory for f_pert_first_comp", OUT_ERROR)
        end if
        allocate(f_pert_tuple(num_all_pert), stat=ierr)
        if (ierr/=0) then
            call f_callback_UserOutput("OpenRSPGetResidue_f>> failed to allocate memory for f_pert_tuple", OUT_ERROR)
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
        ! gets the indices of perturbations w.r.t. residue is calculated
        if (num_props/=1) then
            call f_callback_UserOutput("FIXME>> openrsp_get_residue() does not support multiple properties", OUT_ERROR)
        end if
        ! C memory: [num_props][num_excit][num_freq_configs][pert_tuple][2]
        ! and num_props==1
        resize_per_excit = size_residues/num_excit
        if (order_residue>2) then
            call f_callback_UserOutput("OpenRSPGetResidue_f>> only supports single and double residues", OUT_ERROR)
        end if
        allocate(residue_spec_index(max(residue_num_pert(1),              &
                                        residue_num_pert(order_residue)), &
                                    order_residue),                       &
                 stat=ierr)
        if (ierr/=0) then
            call f_callback_UserOutput("OpenRSPGetResidue_f>> failed to allocate memory for residue_spec_index", OUT_ERROR)
        end if
        jpert = 0
        do iext = 1, order_residue
            do ipert = 1, residue_num_pert(iext)
                jpert = jpert+1
                residue_spec_index(ipert,iext) = residue_idx_pert(jpert)
            end do
        end do
        ! allocates memory for the frequencies of perturbations
        allocate(f_pert_freqs(size(pert_freqs)/num_excit/2), stat=ierr)
        if (ierr/=0) then
            call f_callback_UserOutput("OpenRSPGetResidue_f>> failed to allocate memory for f_pert_freqs", OUT_ERROR)
        end if
        ! gets the matrices
        ierr = QcMat_C_F_POINTER(f_F_unpert, (/F_unpert/))
        if (ierr/=QSUCCESS) then
            call f_callback_UserOutput("OpenRSPGetResidue_f>> failed to call QcMat_C_F_POINTER(F)", OUT_ERROR)
        end if
        ierr = QcMat_C_F_POINTER(f_S_unpert, (/S_unpert/))
        if (ierr/=QSUCCESS) then
            call f_callback_UserOutput("OpenRSPGetResidue_f>> failed to call QcMat_C_F_POINTER(S)", OUT_ERROR)
        end if
        ierr = QcMat_C_F_POINTER(f_D_unpert, (/D_unpert/))
        if (ierr/=QSUCCESS) then
            call f_callback_UserOutput("OpenRSPGetResidue_f>> failed to call QcMat_C_F_POINTER(D)", OUT_ERROR)
        end if
        ! sets the context of callback functions
        call RSP_CTX_Create(rsp_solver,   &
                            c_null_ptr,   &
                            overlap,      &
                            one_oper,     &
                            two_oper,     &
                            xc_fun)
        ! allocates memory for the results
        allocate(f_residues(size_residues), stat=ierr)
        if (ierr/=0) then
            call f_callback_UserOutput("OpenRSPGetResidue_f>> failed to allocate memory for f_residues", OUT_ERROR)
        end if
        f_residues = 0.0

        mem_calibrate = .FALSE.
        ! MaR: max_mat set to very high number to take matrix limitations out of use
        ! during development of other features
        max_mat = 999999999

        jpert = 0
        do iext = 1, num_excit
            ! gets the frequencies of perturbations
            do ipert = 1, size(f_pert_freqs)
                jpert = jpert+1
                f_pert_freqs(ipert) = cmplx(pert_freqs(2*jpert-1), &
                                            pert_freqs(2*jpert),   &
                                            kind=QREAL)
            end do
            ! gets the excitation energies and eigenvectors
            offset_residue = (iext-1)*order_residue
            do ipert = 1, order_residue
                exenerg(ipert) = cmplx(excit_energy(offset_residue+ipert), &
                                       kind=QREAL)
            end do
            ierr = QcMat_C_F_POINTER(X_unpert(1:order_residue),     &
                                     eigen_vector(offset_residue+1: &
                                                  offset_residue+order_residue))
            if (ierr/=QSUCCESS) then
                call f_callback_UserOutput("OpenRSPGetResidue_f>> failed to call QcMat_C_F_POINTER(X)", OUT_ERROR)
            end if
            ! calculates residues for the current excited state
            ipert = (iext-1)*resize_per_excit

            call openrsp_get_property(num_props,                                  &
                                      len_tuple,                                  &
                                      f_pert_dims,                                &
                                      f_pert_first_comp,                          &
                                      f_pert_tuple,                               &
                                      num_freq_configs,                           &
                                      f_pert_freqs,                               &
                                      kn_rules,                                   &
                                      f_F_unpert(1),                              &
                                      f_S_unpert(1),                              &
                                      f_D_unpert(1),                              &
                                      f_callback_RSPSolverGetLinearRSPSolution,   &
                                      f_callback_RSPZeroOperGetContribution,      &
                                      f_callback_RSPOverlapGetMat,                &
                                      f_callback_RSPOverlapGetExp,                &
                                      f_callback_RSPOneOperGetMat,                &
                                      f_callback_RSPOneOperGetExp,                &
                                      f_callback_RSPTwoOperGetMat,                &
                                      f_callback_RSPTwoOperGetExp,                &
                                      f_callback_RSPXCFunGetMat,                  &
                                      f_callback_RSPXCFunGetExp,                  &
                                      f_callback_UserOutput,                      &
                                      r_flag,                                     &
                                      resize_per_excit,                           &
                                      f_residues(ipert+1:ipert+resize_per_excit), &
                                      order_residue,                              &
                                      mem_calibrate=mem_calibrate,                &
                                      max_mat=max_mat,                            &
                                      mem_result=mem_result,                      &
                                      residue_spec_pert=residue_num_pert,         &
                                      size_rsi_1=size(residue_spec_index, 1),     &
                                      residue_spec_index=residue_spec_index,      &
                                      exenerg=exenerg,                            &
                                      Xf_unpert=X_unpert)
        end do

        ! assigns the results
        jpert = 0
        do ipert = 1, size_residues
            jpert = jpert+1
            residues(jpert) = real(f_residues(ipert))
            jpert = jpert+1
            residues(jpert) = aimag(f_residues(ipert))
        end do
        ! cleans up
        ierr = QcMat_C_NULL_PTR(A=f_F_unpert)
        if (ierr/=QSUCCESS) then
            call f_callback_UserOutput("OpenRSPGetResidue_f>> failed to call QcMat_C_NULL_PTR(F)", OUT_ERROR)
        end if
        ierr = QcMat_C_NULL_PTR(A=f_S_unpert)
        if (ierr/=QSUCCESS) then
            call f_callback_UserOutput("OpenRSPGetResidue_f>> failed to call QcMat_C_NULL_PTR(S)", OUT_ERROR)
        end if
        ierr = QcMat_C_NULL_PTR(A=f_D_unpert)
        if (ierr/=QSUCCESS) then
            call f_callback_UserOutput("OpenRSPGetResidue_f>> failed to call QcMat_C_NULL_PTR(D)", OUT_ERROR)
        end if
        ierr = QcMat_C_NULL_PTR(A=X_unpert(1:order_residue))
        if (ierr/=QSUCCESS) then
            call f_callback_UserOutput("OpenRSPGetResidue_f>> failed to call QcMat_C_NULL_PTR(X)", OUT_ERROR)
        end if
        deallocate(f_pert_dims)
        deallocate(f_pert_first_comp)
        deallocate(f_pert_tuple)
        deallocate(residue_spec_index)
        deallocate(f_pert_freqs)
        call RSP_CTX_Destroy()
        deallocate(f_residues)
        return
100     format(A,50I8)
    end subroutine OpenRSPGetResidue_f
