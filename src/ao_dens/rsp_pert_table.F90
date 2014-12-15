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
!!  This file defines the table of perturbations (used for time being).
!!
!!  2014-12-11, Bin Gao
!!  * first version

module rsp_pert_table

    use qmatrix, only: QINT

    implicit none

    public

    integer(kind=QINT), parameter :: RSP_GEO_PERT = 1
    integer(kind=QINT), parameter :: RSP_EL_PERT = 2
    integer(kind=QINT), parameter :: RSP_ELGR_PERT = 3
    integer(kind=QINT), parameter :: RSP_MAG0_PERT = 4
    integer(kind=QINT), parameter :: RSP_MAG_PERT = 5
    character(4), parameter :: CHAR_PERT_TABLE(4) = (/'GEO ', 'EL  ', 'ELGR ', 'MAG0', 'MAG '/)

end module rsp_pert_table
