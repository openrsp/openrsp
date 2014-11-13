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
!!  This file implements subroutines of overlap integrals used in the Fortran APIs.
!!
!!  2014-08-05, Bin Gao
!!  * first version

! basic data types
#include "api/qmatrix_c_type.h"

module rsp_overlap_f

    use, intrinsic :: iso_c_binding
    use qmatrix, only: QINT,QREAL,QMat

    implicit none

    integer(kind=4), private, parameter :: STDOUT = 6

    ! user specified callback subroutines
    abstract interface
        subroutine OverlapGetMat_f(bra_num_pert,      &
                                   bra_perturbations, &
                                   bra_pert_orders,   &
                                   ket_num_pert,      &
                                   ket_perturbations, &
                                   ket_pert_orders,   &
                                   num_pert,          &
                                   perturbations,     &
                                   pert_orders,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                                   len_ctx,           &
                                   user_ctx,          &
#endif
                                   num_int,           &
                                   val_int)
            use qmatrix, only: QINT,QREAL,QMat
            integer(kind=QINT), intent(in) :: bra_num_pert
            integer(kind=QINT), intent(in) :: bra_perturbations(bra_num_pert)
            integer(kind=QINT), intent(in) :: bra_pert_orders(bra_num_pert)
            integer(kind=QINT), intent(in) :: ket_num_pert
            integer(kind=QINT), intent(in) :: ket_perturbations(ket_num_pert)
            integer(kind=QINT), intent(in) :: ket_pert_orders(ket_num_pert)
            integer(kind=QINT), intent(in) :: num_pert
            integer(kind=QINT), intent(in) :: perturbations(num_pert)
            integer(kind=QINT), intent(in) :: pert_orders(num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
            integer(kind=QINT), intent(in) :: len_ctx
            character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
            integer(kind=QINT), intent(in) :: num_int
            type(QMat), intent(inout) :: val_int(num_int)
        end subroutine OverlapGetMat_f
        subroutine OverlapGetExp_f(bra_num_pert,      &
                                   bra_perturbations, &
                                   bra_pert_orders,   &
                                   ket_num_pert,      &
                                   ket_perturbations, &
                                   ket_pert_orders,   &
                                   num_pert,          &
                                   perturbations,     &
                                   pert_orders,       &
                                   num_dens,          &
                                   ao_dens,           &
#if defined(OPENRSP_F_USER_CONTEXT)
                                   len_ctx,           &
                                   user_ctx,          &
#endif
                                   num_exp,           &
                                   val_exp)
            use qmatrix, only: QINT,QREAL,QMat
            integer(kind=QINT), intent(in) :: bra_num_pert
            integer(kind=QINT), intent(in) :: bra_perturbations(bra_num_pert)
            integer(kind=QINT), intent(in) :: bra_pert_orders(bra_num_pert)
            integer(kind=QINT), intent(in) :: ket_num_pert
            integer(kind=QINT), intent(in) :: ket_perturbations(ket_num_pert)
            integer(kind=QINT), intent(in) :: ket_pert_orders(ket_num_pert)
            integer(kind=QINT), intent(in) :: num_pert
            integer(kind=QINT), intent(in) :: perturbations(num_pert)
            integer(kind=QINT), intent(in) :: pert_orders(num_pert)
            integer(kind=QINT), intent(in) :: num_dens
            type(QMat), intent(in) :: ao_dens(num_dens)
#if defined(OPENRSP_F_USER_CONTEXT)
            integer(kind=QINT), intent(in) :: len_ctx
            character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
            integer(kind=QINT), intent(in) :: num_exp
            real(kind=QREAL), intent(inout) :: val_exp(num_exp)
        end subroutine OverlapGetExp_f
    end interface

    ! context of callback subroutines of overlap integrals
    type, public :: OverlapFun_f
        private
#if defined(OPENRSP_F_USER_CONTEXT)
        ! user-defined callback function context
        integer(kind=QINT) :: len_ctx = 0
        character(len=1), allocatable :: user_ctx(:)
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
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_overlap_mat(bra_num_pert,      &
                                       bra_perturbations, &
                                       bra_pert_orders,   &
                                       ket_num_pert,      &
                                       ket_perturbations, &
                                       ket_pert_orders,   &
                                       num_pert,          &
                                       perturbations,     &
                                       pert_orders,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                                       len_ctx,           &
                                       user_ctx,          &
#endif
                                       num_int,           &
                                       val_int)
                use qmatrix, only: QINT,QREAL,QMat
                integer(kind=QINT), intent(in) :: bra_num_pert
                integer(kind=QINT), intent(in) :: bra_perturbations(bra_num_pert)
                integer(kind=QINT), intent(in) :: bra_pert_orders(bra_num_pert)
                integer(kind=QINT), intent(in) :: ket_num_pert
                integer(kind=QINT), intent(in) :: ket_perturbations(ket_num_pert)
                integer(kind=QINT), intent(in) :: ket_pert_orders(ket_num_pert)
                integer(kind=QINT), intent(in) :: num_pert
                integer(kind=QINT), intent(in) :: perturbations(num_pert)
                integer(kind=QINT), intent(in) :: pert_orders(num_pert)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_int
                type(QMat), intent(inout) :: val_int(num_int)
            end subroutine get_overlap_mat
            subroutine get_overlap_exp(bra_num_pert,      &
                                       bra_perturbations, &
                                       bra_pert_orders,   &
                                       ket_num_pert,      &
                                       ket_perturbations, &
                                       ket_pert_orders,   &
                                       num_pert,          &
                                       perturbations,     &
                                       pert_orders,       &
                                       num_dens,          &
                                       ao_dens,           &
#if defined(OPENRSP_F_USER_CONTEXT)
                                       len_ctx,           &
                                       user_ctx,          &
#endif
                                       num_exp,           &
                                       val_exp)
                use qmatrix, only: QINT,QREAL,QMat
                integer(kind=QINT), intent(in) :: bra_num_pert
                integer(kind=QINT), intent(in) :: bra_perturbations(bra_num_pert)
                integer(kind=QINT), intent(in) :: bra_pert_orders(bra_num_pert)
                integer(kind=QINT), intent(in) :: ket_num_pert
                integer(kind=QINT), intent(in) :: ket_perturbations(ket_num_pert)
                integer(kind=QINT), intent(in) :: ket_pert_orders(ket_num_pert)
                integer(kind=QINT), intent(in) :: num_pert
                integer(kind=QINT), intent(in) :: perturbations(num_pert)
                integer(kind=QINT), intent(in) :: pert_orders(num_pert)
                integer(kind=QINT), intent(in) :: num_dens
                type(QMat), intent(in) :: ao_dens(num_dens)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                integer(kind=QINT), intent(in) :: num_exp
                real(kind=QREAL), intent(inout) :: val_exp(num_exp)
            end subroutine get_overlap_exp
        end interface
#if defined(OPENRSP_F_USER_CONTEXT)
        integer(kind=4) ierr  !error information
        overlap_fun%len_ctx = size(user_ctx)
        allocate(overlap_fun%user_ctx(overlap_fun%len_ctx), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPOverlapCreate_f>> length", overlap_fun%len_ctx
            stop "RSPOverlapCreate_f>> failed to allocate memory for user_ctx"
        end if
        overlap_fun%user_ctx = user_ctx
#endif
        overlap_fun%get_overlap_mat => get_overlap_mat
        overlap_fun%get_overlap_exp => get_overlap_exp
    end subroutine RSPOverlapCreate_f

    !% \brief calls Fortran callback subroutine to get integral matrices of overlap integrals
    !  \author Bin Gao
    !  \date 2014-08-05
    !  \param[integer]{in} bra_num_pert number of perturbations on the bra
    !  \param[integer]{in} bra_perturbations the perturbations on the bra
    !  \param[integer]{in} bra_pert_orders orders of the perturbations on the bra
    !  \param[integer]{in} ket_num_pert number of perturbations on the ket
    !  \param[integer]{in} ket_perturbations the perturbations on the ket
    !  \param[integer]{in} ket_pert_orders orders of the perturbations on the ket
    !  \param[integer]{in} num_pert number of perturbations
    !  \param[integer]{in} perturbations the perturbations
    !  \param[integer]{in} pert_orders orders of the perturbations
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !  \param[integer]{in} num_int number of the integral matrices
    !% \param[C_PTR:type]{inout} val_int the integral matrices
    subroutine RSPOverlapGetMat_f(bra_num_pert,      &
                                  bra_perturbations, &
                                  bra_pert_orders,   &
                                  ket_num_pert,      &
                                  ket_perturbations, &
                                  ket_pert_orders,   &
                                  num_pert,          &
                                  perturbations,     &
                                  pert_orders,       &
                                  user_ctx,          &
                                  num_int,           &
                                  val_int)           &
        bind(C, name="RSPOverlapGetMat_f")
        integer(kind=C_QINT), value, intent(in) :: bra_num_pert
        integer(kind=C_QINT), intent(in) :: bra_perturbations(bra_num_pert)
        integer(kind=C_QINT), intent(in) :: bra_pert_orders(bra_num_pert)
        integer(kind=C_QINT), value, intent(in) :: ket_num_pert
        integer(kind=C_QINT), intent(in) :: ket_perturbations(ket_num_pert)
        integer(kind=C_QINT), intent(in) :: ket_pert_orders(ket_num_pert)
        integer(kind=C_QINT), value, intent(in) :: num_pert
        integer(kind=C_QINT), intent(in) :: perturbations(num_pert)
        integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
        type(C_PTR), value, intent(in) :: user_ctx
        integer(kind=C_QINT), value, intent(in) :: num_int
        type(C_PTR), intent(inout) :: val_int(num_int)
        type(OverlapFun_f), pointer :: overlap_fun  !context of callback subroutines
        type(QMat), allocatable :: f_val_int(:)     !integral matrices
        character(len=1), allocatable :: enc(:)     !encoded data as an array of characters
        integer(kind=QINT) len_enc                             !length of encoded data
        integer(kind=4) ierr                                !error information
        integer(kind=QINT) imat                                !incremental recorder over matrices
        ! converts C pointer to Fortran QMat type, inspired by
        ! http://stackoverflow.com/questions/6998995/fortran-array-of-pointer-arrays
        ! and
        ! http://jblevins.org/log/transfer
        allocate(f_val_int(num_int), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPOverlapGetMat_f>> num_int", num_int
            stop "RSPOverlapGetMat_f>> failed to allocate memory for f_val_int"
        end if
        do imat = 1, num_int
            ! encodes the C pointer in a character array
            len_enc = size(transfer(val_int(imat), enc))
            allocate(enc(len_enc), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,"(A,I8)") "RSPOverlapGetMat_f>> length", len_enc
                stop "RSPOverlapGetMat_f>> failed to allocate memory for enc"
            end if
            enc = transfer(val_int(imat), enc)
            ! decodes as QMat type
            f_val_int(imat) = transfer(enc, f_val_int(imat))
            ! cleans up
            deallocate(enc)
        end do
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, overlap_fun)
        ! invokes Fortran callback subroutine to calculate the integral matrices
        call overlap_fun%get_overlap_mat(bra_num_pert,         &
                                         bra_perturbations,    &
                                         bra_pert_orders,      &
                                         ket_num_pert,         &
                                         ket_perturbations,    &
                                         ket_pert_orders,      &
                                         num_pert,             &
                                         perturbations,        &
                                         pert_orders,          &
#if defined(OPENRSP_F_USER_CONTEXT)
                                         overlap_fun%len_ctx,  &
                                         overlap_fun%user_ctx, &
#endif
                                         num_int,              &
                                         f_val_int)
        ! cleans up
        nullify(overlap_fun)
        deallocate(f_val_int)
    end subroutine RSPOverlapGetMat_f

    !% \brief calls Fortran callback subroutine to get expectation values of overlap integrals
    !  \author Bin Gao
    !  \date 2014-08-05
    !  \param[integer]{in} bra_num_pert number of perturbations on the bra
    !  \param[integer]{in} bra_perturbations the perturbations on the bra
    !  \param[integer]{in} bra_pert_orders orders of the perturbations on the bra
    !  \param[integer]{in} ket_num_pert number of perturbations on the ket
    !  \param[integer]{in} ket_perturbations the perturbations on the ket
    !  \param[integer]{in} ket_pert_orders orders of the perturbations on the ket
    !  \param[integer]{in} num_pert number of perturbations
    !  \param[integer]{in} perturbations the perturbations
    !  \param[integer]{in} pert_orders orders of the perturbations
    !  \param[integer]{in} num_dens number of atomic orbital (AO) based density matrices
    !  \param[C_PTR:type]{inout} ao_dens the AO based density matrices
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !  \param[integer]{in} num_exp number of expectation values
    !% \param[real]{out} val_exp the expectation values
    subroutine RSPOverlapGetExp_f(bra_num_pert,      &
                                  bra_perturbations, &
                                  bra_pert_orders,   &
                                  ket_num_pert,      &    
                                  ket_perturbations, &
                                  ket_pert_orders,   &
                                  num_pert,          &
                                  perturbations,     &
                                  pert_orders,       &
                                  num_dens,          &
                                  ao_dens,           &
                                  user_ctx,          &
                                  num_exp,           &
                                  val_exp)           &
        bind(C, name="RSPOverlapGetExp_f")
        integer(kind=C_QINT), value, intent(in) :: bra_num_pert
        integer(kind=C_QINT), intent(in) :: bra_perturbations(bra_num_pert)
        integer(kind=C_QINT), intent(in) :: bra_pert_orders(bra_num_pert)
        integer(kind=C_QINT), value, intent(in) :: ket_num_pert
        integer(kind=C_QINT), intent(in) :: ket_perturbations(ket_num_pert)
        integer(kind=C_QINT), intent(in) :: ket_pert_orders(ket_num_pert)
        integer(kind=C_QINT), value, intent(in) :: num_pert
        integer(kind=C_QINT), intent(in) :: perturbations(num_pert)
        integer(kind=C_QINT), intent(in) :: pert_orders(num_pert)
        integer(kind=C_QINT), value, intent(in) :: num_dens
        type(C_PTR), intent(in) :: ao_dens(num_dens)
        type(C_PTR), value, intent(in) :: user_ctx
        integer(kind=C_QINT), value, intent(in) :: num_exp
        real(C_QREAL), intent(inout) :: val_exp(num_exp)
        type(OverlapFun_f), pointer :: overlap_fun  !context of callback subroutines
        type(QMat), allocatable :: f_ao_dens(:)     !AO based density matrices
        character(len=1), allocatable :: enc(:)     !encoded data as an array of characters
        integer(kind=QINT) len_enc                             !length of encoded data
        integer(kind=4) ierr                                !error information
        integer(kind=QINT) imat                                !incremental recorder over matrices
        ! converts C pointer to Fortran QMat type, inspired by
        ! http://stackoverflow.com/questions/6998995/fortran-array-of-pointer-arrays
        ! and
        ! http://jblevins.org/log/transfer
        allocate(f_ao_dens(num_dens), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPOverlapGetExp_f>> num_dens", num_dens
            stop "RSPOverlapGetExp_f>> failed to allocate memory for f_ao_dens"
        end if
        do imat = 1, num_dens
            ! encodes the C pointer in a character array
            len_enc = size(transfer(ao_dens(imat), enc))
            allocate(enc(len_enc), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,"(A,I8)") "RSPOverlapGetExp_f>> length", len_enc
                stop "RSPOverlapGetExp_f>> failed to allocate memory for enc"
            end if
            enc = transfer(ao_dens(imat), enc)
            ! decodes as QMat type
            f_ao_dens(imat) = transfer(enc, f_ao_dens(imat))
            ! cleans up
            deallocate(enc)
        end do
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, overlap_fun)
        ! invokes Fortran callback subroutine to calculate the expectation values
        call overlap_fun%get_overlap_exp(bra_num_pert,         &
                                         bra_perturbations,    &
                                         bra_pert_orders,      &
                                         ket_num_pert,         &
                                         ket_perturbations,    &
                                         ket_pert_orders,      &
                                         num_pert,             &
                                         perturbations,        &
                                         pert_orders,          &
                                         num_dens,             &
                                         f_ao_dens,            &
#if defined(OPENRSP_F_USER_CONTEXT)
                                         overlap_fun%len_ctx,  &
                                         overlap_fun%user_ctx, &
#endif
                                         num_exp,              &
                                         val_exp)
        ! cleans up
        nullify(overlap_fun)
        deallocate(f_ao_dens)
        return
    end subroutine RSPOverlapGetExp_f

    !% \brief cleans the context of callback subroutines of overlap integrals
    !  \author Bin Gao
    !  \date 2014-08-05
    !% \param[OverlapFun_f:type]{inout} overlap_fun the context of callback subroutines
    subroutine RSPOverlapDestroy_f(overlap_fun)
        type(OverlapFun_f), intent(inout) :: overlap_fun
#if defined(OPENRSP_F_USER_CONTEXT)
        overlap_fun%len_ctx = 0
        deallocate(overlap_fun%user_ctx)
#endif
        nullify(overlap_fun%get_overlap_mat)
        nullify(overlap_fun%get_overlap_exp)
    end subroutine RSPOverlapDestroy_f

end module rsp_overlap_f
