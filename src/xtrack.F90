!!  gen1int: compute derivatives of one-electron integrals using Hermite Gaussians
!!  Copyright 2009 Bin Gao, Andreas Thorvaldsen, Radovan Bast, and Kenneth Ruud
!!
!!  This file is part of gen1int.
!!
!!  gen1int is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  gen1int is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!!  
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with gen1int. If not, see <http://www.gnu.org/licenses/>.
!!
!!  This file tracks the errors
!!
!!  2009-10-24, Bin Gao:
!!  * modified from previous error tracking subroutines

!> \brief module of tracking errors
!> \details aborts due to critical errors and back tracks the subroutines
!> \author Bin Gao
!> \date 2009-10-24
module xtrack
  ! precision
  use xprecision
  implicit none
  !> standard input and output units
  integer, private, parameter :: x_stdinp = 5, x_stdout = 6
  !> maximum length of the name of subroutines
  integer, private, parameter :: x_len_sub = 50
  !> size of the stack of subroutines
  integer, private, parameter :: x_substack_size = 200
  !> stack of subroutines
  character*(x_len_sub), private :: x_sub_stack( x_substack_size )
  !> top of the stack of subroutines
  integer, private :: x_top_subs = 0

  interface xsub_error
    module procedure xsub_err_int
    module procedure xsub_err_real
  end interface xsub_error

  public :: xsub_enter
  public :: xsub_leave

  private :: xsub_track

  contains

  !> \brief error stop with some additional integer information
  !> \author Bin Gao
  !> \date 2009-10-24
  !> \param err_brief is the brief error message
  !> \param err_int is the additional integer information
  !> \param err_detail is the detailed error message
  subroutine xsub_err_int( err_brief, err_int, err_detail )
#include <config.h>
    character*(*), intent(in) :: err_brief
    integer, intent(in) :: err_int(:)
    character*(*), optional, intent(in) :: err_detail
    ! incremental recorder over integer information
    integer ierr
    ! brief error message
    do ierr = 1, size(err_int)
      write(x_stdout,100) err_brief, err_int(ierr)
    end do
    ! dumps the subroutines in stack
    call xsub_track
    ! detailed error message
    if ( present(err_detail) ) write(x_stdout,100) err_detail
    write(x_stdout,100) 'Please report to: '//trim(PACKAGE_BUGREPORT)
    stop
100 format('!! ERROR: ',A,I6)
  end subroutine xsub_err_int

  !> \brief error stop with some additional real information
  !> \author Bin Gao
  !> \date 2009-10-24
  !> \param err_brief is the brief error message
  !> \param err_real is the additional real information
  !> \param err_detail is the detailed error message
  subroutine xsub_err_real( err_brief, err_real, err_detail )
#include <config.h>
    character*(*), intent(in) :: err_brief
    real(xp), intent(in) :: err_real(:)
    character*(*), optional, intent(in) :: err_detail
    ! incremental recorder over real information
    integer ierr
    ! brief error message
    do ierr = 1, size(err_real)
      write(x_stdout,100) err_brief, err_real(ierr)
    end do
    ! dumps the subroutines in stack
    call xsub_track
    ! detailed error message
    if ( present(err_detail) ) write(x_stdout,100) err_detail
    write(x_stdout,100) 'Please report to: '//trim(PACKAGE_BUGREPORT)
    stop
100 format('!! ERROR: ',A,F14.6)
  end subroutine xsub_err_real

  !> \brief dumps the subroutines in stack
  !> \author Bin Gao
  !> \date 2009-10-24
  subroutine xsub_track
    ! incremental recorder over subs.
    integer isub
    ! dumps now
    write(x_stdout,'()')
    write(x_stdout,100) '                      Subs. in Stack                        '
    write(x_stdout,100) '============================================================'
    write(x_stdout,100) 'Name                                              |   Level '
    write(x_stdout,100) '--------------------------------------------------|---------'
    do isub = x_top_subs, 1, -1
      write(x_stdout,100) x_sub_stack(isub)//'|',isub
    end do
    write(x_stdout,100) '============================================================'
    write(x_stdout,'()')
100 format('SUBS >> ',A,I8)
  end subroutine xsub_track

  !> \brief pushes a new subroutine into the stack of subroutines
  !> \author Bin Gao
  !> \date 2009-07-23
  !> \param sub_name is the name of subroutine
  subroutine xsub_enter( sub_name )
    character*(*), intent(in) :: sub_name
    ! new top of the stack
    x_top_subs = x_top_subs + 1
    ! checks if the stack leaks
    if ( x_top_subs > x_substack_size ) then
      write(x_stdout,100) 'Stack leaks! Increse x_substack_size in xtrack.F90!'
      stop
    end if
    ! pushes the new subroutine
    x_sub_stack(x_top_subs) = trim(sub_name)
100 format('!! ERROR: ',A)
  end subroutine xsub_enter

  !> \brief pops the top of stack of subroutines
  !> \author Bin Gao
  !> \date 2009-07-23
  subroutine xsub_leave
    ! checks if the stack is empty
    if ( x_top_subs <= 0 ) then
      write(x_stdout,100) 'Stack is empty! No subs. will be popped!'
      return
    end if
    ! clears the top
    x_sub_stack(x_top_subs) = ''
    ! new top of the stack
    x_top_subs = x_top_subs - 1
100 format('!! WARNING: ',A)
  end subroutine xsub_leave

end module xtrack

