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
!!  This file implements subroutines of perturbations used in the Fortran APIs.
!!
!!  2014-08-18, Bin Gao
!!  * first version

! basic data types
#include "api/qcmatrix_c_type.h"

module rsp_pert_f

    use, intrinsic :: iso_c_binding
    use qcmatrix_f, only: QINT

    implicit none

    integer(kind=4), private, parameter :: STDOUT = 6

    ! user specified callback subroutines
    abstract interface
        subroutine GetPertComp_f(pert_label,      &
                                 pert_order,      &
                                 pert_rank,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                                 len_ctx,         &
                                 user_ctx,        &
#endif
                                 pert_num_comp,   &
                                 pert_components, &
                                 pert_comp_orders)
            use qcmatrix_f, only: QINT
            integer(kind=QINT), intent(in) :: pert_label
            integer(kind=QINT), intent(in) :: pert_order
            integer(kind=QINT), intent(in) :: pert_rank
#if defined(OPENRSP_F_USER_CONTEXT)
            integer(kind=QINT), intent(in) :: len_ctx
            character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
            integer(kind=QINT), intent(out) :: pert_num_comp
            integer(kind=QINT), intent(out) :: pert_components(pert_order)
            integer(kind=QINT), intent(out) :: pert_comp_orders(pert_order)
        end subroutine GetPertComp_f
        subroutine GetPertRank_f(pert_label,       &
                                 pert_num_comp,    &
                                 pert_components,  &
                                 pert_comp_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                 len_ctx,          &
                                 user_ctx,         &
#endif
                                 pert_rank)
            use qcmatrix_f, only: QINT
            integer(kind=QINT), intent(in) :: pert_label
            integer(kind=QINT), intent(in) :: pert_num_comp
            integer(kind=QINT), intent(in) :: pert_components(pert_num_comp)
            integer(kind=QINT), intent(in) :: pert_comp_orders(pert_num_comp)
#if defined(OPENRSP_F_USER_CONTEXT)
            integer(kind=QINT), intent(in) :: len_ctx
            character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
            integer(kind=QINT), intent(out) :: pert_rank
        end subroutine GetPertRank_f
    end interface

    ! context of callback subroutine of response equation solver
    type, public :: PertFun_f
        private
#if defined(OPENRSP_F_USER_CONTEXT)
        ! user-defined callback function context
        integer(kind=QINT) :: len_ctx = 0
        character(len=1), allocatable :: user_ctx(:)
#endif
        ! callback functions
        procedure(GetPertComp_f), nopass, pointer :: get_pert_comp
        procedure(GetPertRank_f), nopass, pointer :: get_pert_rank
    end type PertFun_f

    public :: RSPPertCreate_f
    public :: RSPPertGetComp_f
    public :: RSPPertGetRank_f
    public :: RSPPertDestroy_f

    contains

    !% \brief creates the context of callback subroutines of perturbations
    !  \author Bin Gao
    !  \date 2014-08-18
    !  \param[PertFun_f:type]{inout} pert_fun the context of callback subroutines
    !  \param[character]{in} user_ctx user-defined callback function context
    !  \param[subroutine]{in} get_pert_comp user specified function for
    !      getting components of a perturbation
    !  \param[subroutine]{in} get_pert_rank user specified function for
    !%     getting rank of a perturbation
    subroutine RSPPertCreate_f(pert_fun,      &
#if defined(OPENRSP_F_USER_CONTEXT)
                               user_ctx,      &
#endif
                               get_pert_comp, &
                               get_pert_rank)
        type(PertFun_f), intent(inout) :: pert_fun
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_pert_comp(pert_label,      &
                                     pert_order,      &
                                     pert_rank,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                                     len_ctx,         &
                                     user_ctx,        &
#endif
                                     pert_num_comp,   &
                                     pert_components, &
                                     pert_comp_orders)
                use qcmatrix_f, only: QINT
                integer(kind=QINT), intent(in) :: pert_label
                integer(kind=QINT), intent(in) :: pert_order
                integer(kind=QINT), intent(in) :: pert_rank
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(out) :: pert_num_comp
                integer(kind=QINT), intent(out) :: pert_components(pert_order)
                integer(kind=QINT), intent(out) :: pert_comp_orders(pert_order)
            end subroutine get_pert_comp
            subroutine get_pert_rank(pert_label,       &
                                     pert_num_comp,    &
                                     pert_components,  &
                                     pert_comp_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                     len_ctx,          &
                                     user_ctx,         &
#endif
                                     pert_rank)
                use qcmatrix_f, only: QINT
                integer(kind=QINT), intent(in) :: pert_label
                integer(kind=QINT), intent(in) :: pert_num_comp
                integer(kind=QINT), intent(in) :: pert_components(pert_num_comp)
                integer(kind=QINT), intent(in) :: pert_comp_orders(pert_num_comp)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(out) :: pert_rank
            end subroutine get_pert_rank
        end interface
#if defined(OPENRSP_F_USER_CONTEXT)
        integer(kind=4) ierr  !error information
        pert_fun%len_ctx = size(user_ctx)
        allocate(pert_fun%user_ctx(pert_fun%len_ctx), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPPertCreate_f>> length", pert_fun%len_ctx
            stop "RSPPertCreate_f>> failed to allocate memory for user_ctx"
        end if
        pert_fun%user_ctx = user_ctx
#endif
        pert_fun%get_pert_comp => get_pert_comp
        pert_fun%get_pert_rank => get_pert_rank
    end subroutine RSPPertCreate_f

    !% \brief calls Fortran callback subroutine to get components of a perturbation
    !  \author Bin Gao
    !  \date 2014-08-18
    !  \param[integer]{in} pert_label lable of the perturbation
    !  \param[integer]{in} pert_order the order of the perturbation
    !  \param[integer]{in} pert_rank the rank of the perturbation
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !  \param[integer]{out} pert_num_comp number of components of the perturbation
    !  \param[integer]{out} pert_components components of the perturbation
    !% \param[integer]{out} pert_comp_orders orders of the components
    subroutine RSPPertGetComp_f(pert_label,       &
                                pert_order,       &
                                pert_rank,        &
                                user_ctx,         &
                                pert_num_comp,    &
                                pert_components,  &
                                pert_comp_orders) &
        bind(C, name="RSPPertGetComp_f")
        integer(kind=C_QINT), value, intent(in) :: pert_label
        integer(kind=C_QINT), value, intent(in) :: pert_order
        integer(kind=C_QINT), value, intent(in) :: pert_rank
        type(C_PTR), value, intent(in) :: user_ctx
        integer(kind=C_QINT), intent(out) :: pert_num_comp
        integer(kind=C_QINT), intent(out) :: pert_components(pert_order)
        integer(kind=C_QINT), intent(out) :: pert_comp_orders(pert_order)
        type(PertFun_f), pointer :: pert_fun   !context of callback subroutines
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, pert_fun)
        ! invokes Fortran callback subroutine to get components of the perturbation
        call pert_fun%get_pert_comp(pert_label,        &
                                    pert_order,        &
                                    pert_rank,         &
#if defined(OPENRSP_F_USER_CONTEXT)
                                    pert_fun%len_ctx,  &
                                    pert_fun%user_ctx, &
#endif
                                    pert_num_comp,     &
                                    pert_components,   &
                                    pert_comp_orders)
        ! cleans up
        nullify(pert_fun)
        return
    end subroutine RSPPertGetComp_f

    !% \brief calls Fortran callback subroutine to get the rank of
    !      a perturbation with its components
    !  \author Bin Gao
    !  \date 2014-08-18
    !  \param[integer]{in} pert_label lable of the perturbation
    !  \param[integer]{in} pert_num_comp number of components of the perturbation
    !  \param[integer]{in} pert_components components of the perturbation
    !  \param[integer]{in} pert_comp_orders orders of the components
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !% \param[integer]{out} pert_rank the rank of the perturbation
    subroutine RSPPertGetRank_f(pert_label,       &
                                pert_num_comp,    &
                                pert_components,  &
                                pert_comp_orders, &
                                user_ctx,         &
                                pert_rank)        &
        bind(C, name="RSPPertGetRank_f")
        integer(kind=C_QINT), value, intent(in) :: pert_label
        integer(kind=C_QINT), value, intent(in) :: pert_num_comp
        integer(kind=C_QINT), intent(in) :: pert_components(pert_num_comp)
        integer(kind=C_QINT), intent(in) :: pert_comp_orders(pert_num_comp)
        type(C_PTR), value, intent(in) :: user_ctx
        integer(kind=C_QINT), intent(out) :: pert_rank
        type(PertFun_f), pointer :: pert_fun   !context of callback subroutines
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, pert_fun)
        ! invokes Fortran callback subroutine to get the rank of the perturbation
        call pert_fun%get_pert_rank(pert_label,        &
                                    pert_num_comp,     &
                                    pert_components,   &
                                    pert_comp_orders,  &
#if defined(OPENRSP_F_USER_CONTEXT)
                                    pert_fun%len_ctx,  &
                                    pert_fun%user_ctx, &
#endif
                                    pert_rank)
        ! cleans up
        nullify(pert_fun)
        return
    end subroutine RSPPertGetRank_f

    !% \brief cleans the context of callback subroutines of perturbations
    !  \author Bin Gao
    !  \date 2014-08-18
    !% \param[PertFun_f:type]{inout} pert_fun the context of callback subroutines
    subroutine RSPPertDestroy_f(pert_fun)
        type(PertFun_f), intent(inout) :: pert_fun
#if defined(OPENRSP_F_USER_CONTEXT)
        pert_fun%len_ctx = 0
        deallocate(pert_fun%user_ctx)
#endif
        nullify(pert_fun%get_pert_comp)
        nullify(pert_fun%get_pert_rank)
    end subroutine RSPPertDestroy_f

end module rsp_pert_f
