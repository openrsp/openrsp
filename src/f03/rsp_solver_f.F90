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
!!  This file implements subroutines of response equation solver used in the Fortran APIs.
!!
!!  2014-08-06, Bin Gao
!!  * first version

! basic data types
#include "api/qcmatrix_c_type.h"

module rsp_solver_f

    use, intrinsic :: iso_c_binding
    use qcmatrix_f, only: QINT,QREAL,QcMat

    implicit none

    integer(kind=4), private, parameter :: STDOUT = 6

    ! user specified callback subroutine
    abstract interface
        subroutine SolverRun_f(ref_ham,       &
                               ref_state,     &
                               ref_overlap,   &
                               num_freq_sums, &
                               freq_sums,     &
                               size_pert,     &
                               RHS_mat,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                               len_ctx,       &
                               user_ctx,      &
#endif
                               rsp_param)
            use qcmatrix_f, only: QINT,QREAL,QcMat
            type(QcMat), intent(in) :: ref_ham
            type(QcMat), intent(in) :: ref_state
            type(QcMat), intent(in) :: ref_overlap
            integer(kind=QINT), intent(in) :: num_freq_sums
            real(kind=QREAL), intent(in) :: freq_sums(num_freq_sums)
            integer(kind=QINT), intent(in) :: size_pert
            type(QcMat), intent(in) :: RHS_mat(size_pert*num_freq_sums)
#if defined(OPENRSP_F_USER_CONTEXT)
            integer(kind=QINT), intent(in) :: len_ctx
            character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
            type(QcMat), intent(inout) :: rsp_param(size_pert*num_freq_sums)
        end subroutine SolverRun_f
    end interface

    ! context of callback subroutine of response equation solver
    type, public :: SolverFun_f
        private
#if defined(OPENRSP_F_USER_CONTEXT)
        ! user-defined callback function context
        integer(kind=QINT) :: len_ctx = 0
        character(len=1), allocatable :: user_ctx(:)
#endif
        ! callback function
        procedure(SolverRun_f), nopass, pointer :: get_rsp_solution
    end type SolverFun_f

    public :: RSPSolverCreate_f
    public :: RSPSolverGetSolution_f
    public :: RSPSolverDestroy_f

    contains

    !% \brief creates the context of callback subroutine of response equation solver
    !  \author Bin Gao
    !  \date 2014-08-06
    !  \param[SolverFun_f:type]{inout} solver_fun the context of callback subroutine
    !  \param[character]{in} user_ctx user-defined callback function context
    !  \param[subroutine]{in} get_rsp_solution user specified function of
    !%     response equation solver
    subroutine RSPSolverCreate_f(solver_fun, &
#if defined(OPENRSP_F_USER_CONTEXT)
                                 user_ctx,   &
#endif
                                 get_rsp_solution)
        type(SolverFun_f), intent(inout) :: solver_fun
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1), intent(in) :: user_ctx(:)
#endif
        interface
            subroutine get_rsp_solution(ref_ham,       &
                                        ref_state,     &
                                        ref_overlap,   &
                                        num_freq_sums, &
                                        freq_sums,     &
                                        size_pert,     &
                                        RHS_mat,       &
#if defined(OPENRSP_F_USER_CONTEXT)
                                        len_ctx,       &
                                        user_ctx,      &
#endif
                                        rsp_param)
                use qcmatrix_f, only: QINT,QREAL,QcMat
                type(QcMat), intent(in) :: ref_ham
                type(QcMat), intent(in) :: ref_state
                type(QcMat), intent(in) :: ref_overlap
                integer(kind=QINT), intent(in) :: num_freq_sums
                real(kind=QREAL), intent(in) :: freq_sums(num_freq_sums)
                integer(kind=QINT), intent(in) :: size_pert
                type(QcMat), intent(in) :: RHS_mat(size_pert*num_freq_sums)
#if defined(OPENRSP_F_USER_CONTEXT)
                integer(kind=QINT), intent(in) :: len_ctx
                character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
                type(QcMat), intent(inout) :: rsp_param(size_pert*num_freq_sums)
            end subroutine get_rsp_solution
        end interface
#if defined(OPENRSP_F_USER_CONTEXT)
        integer(kind=4) ierr  !error information
        solver_fun%len_ctx = size(user_ctx)
        allocate(solver_fun%user_ctx(solver_fun%len_ctx), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPSolverCreate_f>> length", solver_fun%len_ctx
            stop "RSPSolverCreate_f>> failed to allocate memory for user_ctx"
        end if
        solver_fun%user_ctx = user_ctx
#endif
        solver_fun%get_rsp_solution => get_rsp_solution
    end subroutine RSPSolverCreate_f

    !% \brief calls Fortran callback subroutine to get solution of response equation
    !  \author Bin Gao
    !  \date 2014-08-06
    !  \param[C_PTR:type]{in} ref_ham Hamiltonian of referenced state
    !  \param[C_PTR:type]{in} ref_state electronic state of referenced state
    !  \param[C_PTR:type]{in} ref_overlap overlap integral matrix of referenced state
    !  \param[integer]{in} num_freq_sums number of frequency sums on the left hand side
    !  \param[real]{in} freq_sums the frequency sums on the left hand side
    !  \param[integer]{in} size_pert size of perturbaed matrices
    !  \param[C_PTR:type]{in} RHS_mat RHS matrices, size is \var{size_pert}*\var{num_freq_sums}
    !  \param[C_PTR:type]{in} user_ctx user-defined callback function context
    !% \param[C_PTR:type]{out} rsp_param solved response parameters, size is \var{size_pert}*\var{num_freq_sums}
    subroutine RSPSolverGetSolution_f(ref_ham,       &
                                      ref_state,     &
                                      ref_overlap,   &
                                      num_freq_sums, &
                                      freq_sums,     &
                                      size_pert,     &
                                      RHS_mat,       &
                                      user_ctx,      &
                                      rsp_param)     &
        bind(C, name="RSPSolverGetSolution_f")
        type(C_PTR), intent(in) :: ref_ham
        type(C_PTR), intent(in) :: ref_state
        type(C_PTR), intent(in) :: ref_overlap
        integer(kind=C_QINT), value, intent(in) :: num_freq_sums
        real(kind=C_QREAL), intent(in) :: freq_sums(num_freq_sums)
        integer(kind=C_QINT), value, intent(in) :: size_pert
        type(C_PTR), intent(in) :: RHS_mat(size_pert*num_freq_sums)
        type(C_PTR), value, intent(in) :: user_ctx
        type(C_PTR), intent(inout) :: rsp_param(size_pert*num_freq_sums)
        type(SolverFun_f), pointer :: solver_fun   !context of callback subroutine
        type(QcMat), pointer :: f_ref_ham           !Hamiltonian of referenced state
        type(QcMat), pointer :: f_ref_state         !electronic state of referenced state
        type(QcMat), pointer :: f_ref_overlap       !overlap integral matrix of referenced state
        integer(kind=QINT) size_solution           !size of solution of response equation
        type(QcMat), allocatable :: f_RHS_mat(:)    !RHS matrices
        type(QcMat), allocatable :: f_rsp_param(:)  !response parameters
        character(len=1), allocatable :: enc(:)    !encoded data as an array of characters
        integer(kind=QINT) len_enc                 !length of encoded data
        integer(kind=4) ierr                       !error information
        integer(kind=QINT) imat                    !incremental recorder over matrices
        ! gets the matrices of referenced state
        call c_f_pointer(ref_ham, f_ref_ham)
        call c_f_pointer(ref_state, f_ref_state)
        call c_f_pointer(ref_overlap, f_ref_overlap)
        ! converts C pointer to Fortran QcMat type, inspired by
        ! http://stackoverflow.com/questions/6998995/fortran-array-of-pointer-arrays
        ! and
        ! http://jblevins.org/log/transfer
        size_solution = size_pert*num_freq_sums
        allocate(f_RHS_mat(size_solution), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPSolverGetSolution_f>> size_solution", &
                                   size_solution
            stop "RSPSolverGetSolution_f>> failed to allocate memory for f_RHS_mat"
        end if
        do imat = 1, size_solution
            ! encodes the C pointer in a character array
            len_enc = size(transfer(RHS_mat(imat), enc))
            allocate(enc(len_enc), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,"(A,I8)") "RSPSolverGetSolution_f>> length", len_enc
                stop "RSPSolverGetSolution_f>> failed to allocate memory for enc"
            end if
            enc = transfer(RHS_mat(imat), enc)
            ! decodes as QcMat type
            f_RHS_mat(imat) = transfer(enc, f_RHS_mat(imat))
            ! cleans up
            deallocate(enc)
        end do
        allocate(f_rsp_param(size_solution), stat=ierr)
        if (ierr/=0) then
            write(STDOUT,"(A,I8)") "RSPSolverGetSolution_f>> size_solution", &
                                   size_solution
            stop "RSPSolverGetSolution_f>> failed to allocate memory for f_rsp_param"
        end if
        do imat = 1, size_solution
            ! encodes the C pointer in a character array
            len_enc = size(transfer(rsp_param(imat), enc))
            allocate(enc(len_enc), stat=ierr)
            if (ierr/=0) then
                write(STDOUT,"(A,I8)") "RSPSolverGetSolution_f>> length", len_enc
                stop "RSPSolverGetSolution_f>> failed to allocate memory for enc"
            end if
            enc = transfer(rsp_param(imat), enc)
            ! decodes as QcMat type
            f_rsp_param(imat) = transfer(enc, f_rsp_param(imat))
            ! cleans up
            deallocate(enc)
        end do
        ! gets the Fortran callback subroutine
        call c_f_pointer(user_ctx, solver_fun)
        ! invokes Fortran callback subroutine to solve the response equation
        call solver_fun%get_rsp_solution(f_ref_ham,           &
                                         f_ref_state,         &
                                         f_ref_overlap,       &
                                         num_freq_sums,       &
                                         freq_sums,           &
                                         size_pert,           &
                                         f_RHS_mat,           &
#if defined(OPENRSP_F_USER_CONTEXT)
                                         solver_fun%len_ctx,  &
                                         solver_fun%user_ctx, &
#endif
                                         f_rsp_param)
        ! cleans up
        nullify(solver_fun)
        deallocate(f_rsp_param)
        deallocate(f_RHS_mat)
        nullify(f_ref_overlap)
        nullify(f_ref_state)
        nullify(f_ref_ham)
        return
    end subroutine RSPSolverGetSolution_f

    !% \brief cleans the context of callback subroutine of response equation solver
    !  \author Bin Gao
    !  \date 2014-08-06
    !% \param[SolverFun_f:type]{inout} solver_fun the context of callback subroutine
    subroutine RSPSolverDestroy_f(solver_fun)
        type(SolverFun_f), intent(inout) :: solver_fun
#if defined(OPENRSP_F_USER_CONTEXT)
        solver_fun%len_ctx = 0
        deallocate(solver_fun%user_ctx)
#endif
        nullify(solver_fun%get_rsp_solution)
    end subroutine RSPSolverDestroy_f

end module rsp_solver_f
