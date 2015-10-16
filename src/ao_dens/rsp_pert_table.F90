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

    use RSPPertBasicTypes_f, only: QcPertInt

    implicit none

    public

    integer(kind=QcPertInt), parameter :: RSP_GEO_PERT = 1_QcPertInt
    integer(kind=QcPertInt), parameter :: RSP_EL_PERT = 2_QcPertInt
    integer(kind=QcPertInt), parameter :: RSP_ELGR_PERT = 3_QcPertInt
    integer(kind=QcPertInt), parameter :: RSP_MAG0_PERT = 4_QcPertInt
    integer(kind=QcPertInt), parameter :: RSP_MAG_PERT = 5_QcPertInt
    character(4), parameter :: CHAR_PERT_TABLE(5) = (/'GEO ', 'EL  ', 'ELGR', 'MAG0', 'MAG '/)

end module rsp_pert_table
