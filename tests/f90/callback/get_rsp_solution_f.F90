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
!!  This file implements the callback subroutine get_rsp_solution_f.
!!
!!  2014-08-03, Bin Gao
!!  * first version

! configuration file of QcMatrix library
#include "qcmatrix_config.h"

#define OPENRSP_F_TEST_SRC "tests/f90/callback/get_rsp_solution_f.F90"

    subroutine get_rsp_solution_f(ref_ham,       &
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
        implicit none
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
        return
    end subroutine get_rsp_solution_f

#undef OPENRSP_F_TEST_SRC
