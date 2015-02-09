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
!!  This file implements the callback subroutine get_two_oper_exp_f.
!!
!!  2014-08-03, Bin Gao
!!  * first version

! configuration file of QcMatrix library
#include "qcmatrix_config.h"

#define OPENRSP_F_TEST_SRC "tests/f90/callback/get_two_oper_exp_f.F90"

    subroutine get_two_oper_exp_f(num_pert,       &
                                  pert_labels,    &
                                  pert_orders,    &
                                  num_var_dens,   &
                                  var_ao_dens,    &
                                  num_contr_dens, &
                                  contr_ao_dens,  &
#if defined(OPENRSP_F_USER_CONTEXT)
                                  len_ctx,        &
                                  user_ctx,       &
#endif
                                  num_exp,        &
                                  val_exp)
        use qcmatrix_f, only: QINT,  &
                              QREAL, &
                              QcMat
        implicit none
        integer(kind=QINT), intent(in) :: num_pert
        integer(kind=QINT), intent(in) :: pert_labels(num_pert)
        integer(kind=QINT), intent(in) :: pert_orders(num_pert)
        integer(kind=QINT), intent(in) :: num_var_dens
        type(QcMat), intent(in) :: var_ao_dens(num_var_dens)
        integer(kind=QINT), intent(in) :: num_contr_dens
        type(QcMat), intent(in) :: contr_ao_dens(num_contr_dens)
#if defined(OPENRSP_F_USER_CONTEXT)
        integer(kind=QINT), intent(in) :: len_ctx
        character(len=1), intent(in) :: user_ctx(len_ctx)
#endif
        integer(kind=QINT), intent(in) :: num_exp
        real(kind=QREAL), intent(inout) :: val_exp(2*num_exp)
! defined perturbations and their maximum orders
#include "tests/openrsp_f_perturbations.h90"
#if defined(OPENRSP_F_USER_CONTEXT)
        character(len=1) :: twoel_context(6) = (/"N","O","N","L","A","O"/)
#endif
        integer(kind=4) ierr
#if defined(OPENRSP_F_USER_CONTEXT)
        if (all(user_ctx==twoel_contex)) then
            write(6,100) "not implemented"
            call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
        else
            ! electric fields
            if (num_pert==1 .and. pert_labels(1)==PERT_DIPOLE) then
                val_exp = 0
            else
                write(6,100) "not implemented"
                call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
            end if
        end if
#else
        ! electric fields
        if (num_pert==1 .and. pert_labels(1)==PERT_DIPOLE) then
            val_exp = 0
        else
            write(6,100) "not implemented"
            call QErrorExit(6, __LINE__, OPENRSP_F_TEST_SRC)
        end if
#endif
        return
100     format("get_two_oper_exp_f>> ",A)
    end subroutine get_two_oper_exp_f

#undef OPENRSP_F_TEST_SRC
