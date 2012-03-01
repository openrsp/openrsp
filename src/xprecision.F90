!!  gen1int: compute derivatives of one-electron integrals using Hermite Gaussians
!!  Copyright 2009 Bin Gao, Andreas Thorvaldsen, Radovan Bast, and Kenneth Ruud
!!
!!  This file is part of gen1int.
!!
!!  gen1int is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as
!!  published by the Free Software Foundation, either version 3 of
!!  the License, or (at your option) any later version.
!!
!!  gen1int is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public
!!  License along with gen1int. If not, see <http://www.gnu.org/licenses/>.
!!
!!  This file defines the precision of real numbers, and some tranformation functions.
!!
!!  2009-08-21, Bin Gao:
!!  * add function xfloat to transfer integer to real number
!!
!!  2009-04-08, Bin Gao:
!!  * first version

!> \brief module of the precision of real numbers
!> \details defines the precision of real numbers, and some tranformation functions
!> \author Bin Gao
!> \date 2009-04-08
module xprecision
  implicit none
  ! single precision
  integer, parameter :: sgl_t = kind(1.0)
  ! double precision
  integer, parameter :: dbl_t = kind(1.0D+00)
  ! used precision
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
  integer, parameter :: xp = kind(1.0)
#else
  integer, parameter :: xp = kind(1.0D+00)
#endif

  public :: xfloat

  contains

  !> \brief transfers an integer to real number
  !> \author Bin Gao
  !> \date 2009-08-21
  !> \param this_int is the integer to be transferred
  !> \return this_float is the transferred real number
  function xfloat( this_int ) result( this_float )
    integer, intent(in) :: this_int
    real(xp) this_float
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL) || defined (SYS_T90)
    !> single precision
    this_float = float( this_int )
#else
    !> double precision
    this_float = dfloat( this_int )
#endif
  end function xfloat

end module xprecision

