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
!!  This file tests the APIs of OpenRSP library.
!!
!!  2014-08-03, Bin Gao
!!  * first version

! parameters for the test suite
#include "tests/openrsp_test_param.h"

#define OPENRSP_F_TEST_SRC "tests/f90/test_f_OpenRSP.F90"

#if defined(OPENRSP_TEST_EXECUTABLE)
    program test_f_OpenRSP
#else
    subroutine test_f_OpenRSP(io_log)
#endif
        use qcmatrix_f, only: QINT,QREAL
        use openrsp_f, only: OpenRSP,                         &
                             OpenRSPCreate_f,                 &
#if defined(OPENRSP_PERTURBATION_FREE)
                             OpenRSPSetPerturbations_f,       &
#endif
                             OpenRSPSetNucGeoPerturbations_f, &
                             OpenRSPSetNucScalarPotential_f,  &
                             OpenRSPSetNucVectorPotential_f,  &
                             OpenRSPDestroy_f
        implicit none
        ! IO of standard output
#if defined(OPENRSP_TEST_EXECUTABLE)
        integer(kind=4), parameter :: io_log = 6
#else
        integer(kind=4), intent(in) :: io_log
#endif
! defined perturbations and their maximum orders
#include "tests/openrsp_f_perturbations.h90"
! atoms and origins
#include "tests/openrsp_f_molecule.h90"
        ! context of response theory calculations
        type(OpenRSP) open_rsp
#if defined(OPENRSP_PERTURBATION_FREE)
        ! labels of all perturbations
        integer(kind=QINT), parameter :: ALL_PERT_LABELS(NUM_ALL_PERT) = (/ &
            PERT_GEOMETRIC,PERT_DIPOLE,PERT_MAGNETIC/)
        ! maximum allowed orders of all perturbations
        integer(kind=QINT), parameter :: ALL_PERT_MAX_ORDERS(NUM_ALL_PERT) = (/ &
            MAX_ORDER_GEOMETRIC,MAX_ORDER_DIPOLE,MAX_ORDER_MAGNETIC/)
        ! sizes of all perturbations up to their maximum orders
        integer(kind=QINT), parameter :: ALL_PERT_SIZES(15) = (/ &
            12,78,364,1365,4368,12376,31824,                     &  !geometric derivatives (4 atoms)
            3,                                                   &  !electric dipole
            3,6,10,15,21,28,36/)                                    !magnetic derivatives
#if defined(OPENRSP_F_USER_CONTEXT)
        ! user defined context for perturbations
        character(len=1) :: pert_context(7) = (/"N","R","N","Z","G","E","O"/)
#endif
        ! callback subroutines getting the components and rank of a perturbation
        external get_pert_comp_f
        external get_pert_rank_f
#endif
        ! error information
        integer(kind=4) ierr

        ! creates the context of response theory calculations
        ierr = OpenRSPCreate_f(open_rsp)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPCreate_f() passed"

#if defined(OPENRSP_PERTURBATION_FREE)
        ! sets information of all perturbations
        ierr = OpenRSPSetPerturbations_f(open_rsp,            &
                                         NUM_ALL_PERT,        &
                                         ALL_PERT_LABELS,     &
                                         ALL_PERT_MAX_ORDERS, &
                                         ALL_PERT_SIZES,      &
#if defined(OPENRSP_F_USER_CONTEXT)
                                         pert_context,        &
#endif
                                         get_pert_comp_f,     &
                                         get_pert_rank_f)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPSetPerturbations_f() passed"
#endif

        ! sets the geometric perturbations for nuclear Hamiltonian
        ierr = OpenRSPSetNucGeoPerturbations_f(open_rsp,   &
                                               NUM_ATOMS,  &
                                               ATOM_COORD, &
                                               ATOM_CHARGE)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPSetNucGeoPerturbations_f() passed"

        ! sets the terms in nuclear Hamiltonian due to the scalar potential
        ierr = OpenRSPSetNucScalarPotential_f(open_rsp, &
                                              DIPOLE_ORIGIN)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPSetNucScalarPotential_f() passed"

        ! sets the terms in nuclear Hamiltonian due to the vector potential
        ierr = OpenRSPSetNucVectorPotential_f(open_rsp, &
                                              GAUGE_ORIGIN)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPSetNucVectorPotential_f() passed"

        ! tests the density matrix-based response theory
        call test_f_OpenRSP_AO(open_rsp, io_log)
        write(io_log,100) "density matrix-based response theory passed"

        ! destroys the context of response theory calculations
        ierr = OpenRSPDestroy_f(open_rsp)
        call QErrorCheckCode(io_log, ierr, __LINE__, OPENRSP_F_TEST_SRC)
        write(io_log,100) "OpenRSPDestroy_f() passed"

100     format("test_f_OpenRSP>> ",A)
#if defined(OPENRSP_TEST_EXECUTABLE)
    end program test_f_OpenRSP
#else
        return
    end subroutine test_f_OpenRSP
#endif

#undef OPENRSP_F_TEST_SRC
