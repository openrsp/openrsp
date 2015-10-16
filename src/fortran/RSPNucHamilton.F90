!! OpenRSP: open-ended library for response theory
!! Copyright 2015 Radovan Bast,
!!                Daniel H. Friese,
!!                Bin Gao,
!!                Dan J. Jonsson,
!!                Magnus Ringholm,
!!                Kenneth Ruud,
!!                Andreas Thorvaldsen
!!
!! OpenRSP is free software: you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as
!! published by the Free Software Foundation, either version 3 of
!! the License, or (at your option) any later version.
!!
!! OpenRSP is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!! GNU Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public
!! License along with OpenRSP. If not, see <http://www.gnu.org/licenses/>.


!!  2015-06-23, Bin Gao
!!  * first version

! basic data types
#include "api/qcmatrix_c_type.h"

#define OPENRSP_API_SRC "src/fortran/RSPNucHamilton.F90"

module RSPNucHamilton_f

    use, intrinsic :: iso_c_binding
    use qcmatrix_f, only: QINT,QREAL
    use RSPPertBasicTypes_f, only: QcPertInt, &
                                   C_QCPERTINT

    implicit none

    integer(kind=4), private, parameter :: STDOUT = 6

    ! user specified callback subroutine
    abstract interface
        subroutine NucHamiltonGetContrib_f(nuc_num_pert,    &
                                           nuc_pert_labels, &
                                           nuc_pert_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                           len_ctx,         &
                                           user_ctx,        &
#endif
                                           size_pert,       &
                                           val_nuc)
            use qcmatrix_f, only: QINT,QREAL
            use RSPPertBasicTypes_f, only: QcPertInt
            integer(kind=QINT), intent(in) :: nuc_num_pert
            integer(kind=QcPertInt), intent(in) :: nuc_pert_labels(nuc_num_pert)
            integer(kind=QINT), intent(in) :: nuc_pert_orders(nuc_num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
            integer(kind=QINT), intent(in) :: len_ctx
            character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
            integer(kind=QINT), intent(in) :: size_pert
            real(kind=QREAL), intent(inout) :: val_nuc(size_pert)
        end subroutine NucHamiltonGetContrib_f
    end interface

    ! context of callback subroutine of nuclear Hamiltonian
    type, public :: NucHamiltonFun_f
        private
#if defined(OPENRSP_F_USER_CONTEXT)
        ! user-defined callback function context
        integer(kind=QINT) :: len_ctx = 0
        character(len=1), allocatable :: user_ctx(:)
#endif
        ! callback function
        procedure(NucHamiltonGetContrib_f), nopass, pointer :: get_nuc_contrib
    end type NucHamiltonFun_f

    public :: RSPNucHamiltonCreate_f
    public :: RSPNucHamiltonGetContributions_f
    public :: RSPNucHamiltonDestroy_f

    contains

    !% \brief creates the context of callback subroutine of nuclear Hamiltonian
    !  \author Bin Gao
    !  \date 2015-06-23
    !  \param[NucHamiltonFun_f:type]{inout} nuc_hamilton_fun the context of callback subroutine
    !  \param[character]{in} user_ctx user-defined callback function context
    !  \param[subroutine]{in} get_nuc_contrib user specified function for
    !%     getting nuclear contributions
    subroutine RSPNucHamiltonCreate_f(nuc_hamilton_fun, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                      user_ctx,         &
#endif
                                      get_nuc_contrib)
        type(NucHamiltonFun_f), intent(inout) :: nuc_hamilton_fun
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_nuc_contrib(nuc_num_pert,    &
                                       nuc_pert_labels, &
                                       nuc_pert_orders, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                       len_ctx,         &
                                       user_ctx,        &
#endif
                                       size_pert,       &
                                       val_nuc)
                use qcmatrix_f, only: QINT,QREAL
                use RSPPertBasicTypes_f, only: QcPertInt
                integer(kind=QINT), intent(in) :: nuc_num_pert
                integer(kind=QcPertInt), intent(in) :: nuc_pert_labels(nuc_num_pert)
                integer(kind=QINT), intent(in) :: nuc_pert_orders(nuc_num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: size_pert
                real(kind=QREAL), intent(inout) :: val_nuc(size_pert)
            end subroutine get_nuc_contrib
        end interface
#if defined(OPENRSP_F_USER_CONTEXT)
        integer(kind=4) ierr  !error information
        nuc_hamilton_fun%len_ctx = size(user_ctx)
        allocate(nuc_hamilton_fun%user_ctx(nuc_hamilton_fun%len_ctx), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPNucHamiltonCreate_f>> length", nuc_hamilton_fun%len_ctx
            stop "RSPNucHamiltonCreate_f>> failed to allocate memory for user_ctx"
        end if
        nuc_hamilton_fun%user_ctx = user_ctx
#endif
        nuc_hamilton_fun%get_nuc_contrib => get_nuc_contrib
    end subroutine RSPNucHamiltonCreate_f

    !% \brief calls Fortran callback subroutine to get nuclear contributions
    !  \author Bin Gao
    !  \date 2015-06-23
    !  \param[integer]{in} len_tuple length of perturbation tuple on the nuclear Hamiltonian
    !  \param[integer]{in} pert_tuple perturbation tuple on the nuclear Hamiltonian
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !  \param[integer]{in} size_pert size of the perturbations on the nuclear Hamiltonian
    !% \param[real]{out} val_nuc the nuclear contributions
    subroutine RSPNucHamiltonGetContributions_f(nuc_num_pert,    &
                                                nuc_pert_labels, &
                                                nuc_pert_orders, &
                                                user_ctx,        &
                                                size_pert,       &
                                                val_nuc)         &
        bind(C, name="RSPNucHamiltonGetContributions_f")
        integer(kind=C_QINT), value, intent(in) :: nuc_num_pert
        integer(kind=C_QCPERTINT), intent(in) :: nuc_pert_labels(nuc_num_pert)
        integer(kind=C_QINT), intent(in) :: nuc_pert_orders(nuc_num_pert)
        type(C_PTR), value, intent(in) :: user_ctx
        integer(kind=C_QINT), value, intent(in) :: size_pert
        real(C_QREAL), intent(inout) :: val_nuc(size_pert)
        type(NucHamiltonFun_f), pointer :: nuc_hamilton_fun  !context of callback subroutine
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, nuc_hamilton_fun)
        ! invokes Fortran callback subroutine to calculate the nuclear contributions
        call nuc_hamilton_fun%get_nuc_contrib(nuc_num_pert,              &
                                              nuc_pert_labels,           &
                                              nuc_pert_orders,           &
#if defined(OPENRSP_F_USER_CONTEXT)
                                              nuc_hamilton_fun%len_ctx,  &
                                              nuc_hamilton_fun%user_ctx, &
#endif
                                              size_pert,                 &
                                              val_nuc)
        ! cleans up
        nullify(nuc_hamilton_fun)
        return
    end subroutine RSPNucHamiltonGetContributions_f

    !% \brief cleans the context of callback subroutine of nuclear Hamiltonian
    !  \author Bin Gao
    !  \date 2015-06-23
    !% \param[NucHamiltonFun_f:type]{inout} nuc_hamilton_fun the context of callback subroutine
    subroutine RSPNucHamiltonDestroy_f(nuc_hamilton_fun)
        type(NucHamiltonFun_f), intent(inout) :: nuc_hamilton_fun
#if defined(OPENRSP_F_USER_CONTEXT)
        nuc_hamilton_fun%len_ctx = 0
        deallocate(nuc_hamilton_fun%user_ctx)
#endif
        nullify(nuc_hamilton_fun%get_nuc_contrib)
    end subroutine RSPNucHamiltonDestroy_f

end module RSPNucHamilton_f

#undef OPENRSP_API_SRC

