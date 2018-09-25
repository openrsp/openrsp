!!  OpenRSP: open-ended library for response theory
!!  Copyright 2014
!!
!!  This source code form is subject to the terms of the
!!  GNU Lesser General Public License, version 2.1.
!!  If a copy of the GNU LGPL v2.1 was not distributed with this
!!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
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
