! Copyright 2012 Magnus Ringholm
!
! This source code form is subject to the terms of the
! GNU Lesser General Public License, version 2.1.
! If a copy of the GNU LGPL v2.1 was not distributed with this
! code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.

! Contains routines and functions for selecting the best rule for 
! the response property calculation.

module rsp_choose_rule

  use rsp_field_tuple, only: p_tuple

  implicit none

  public get_bestkn
  public kn_prod
  public rm_dim

  private

  contains
  
  function get_bestkn(pert)

    implicit none

    type(p_tuple) :: pert
    integer, dimension(2) :: get_bestkn
    integer :: min_n, n, csize, minsize

    ! Do the case n = 1 first to establish one value

    minsize = kn_prod(pert%npert - 1, 1, pert%pdim(2:pert%npert), &
              pert%npert - 1 - 1, pert%pdim(1)) 


    minsize = minsize + kn_prod(pert%npert - 1, &
              0, pert%pdim(2:pert%npert), 1, 1)

    min_n = 1

    ! Get the products for the pert dimensions as dictated by k and n
    ! Assume that integer division rounds down to nearest lower integer

    do n = (pert%npert/2), pert%npert - 1

       csize = kn_prod(pert%npert - 1, 1, pert%pdim(2:pert%npert), &
                 pert%npert - n - 1, pert%pdim(1))


       csize = csize + kn_prod(pert%npert - 1, &
                 0, pert%pdim(2:pert%npert), n, 1)

       ! If the products are smaller than the previous min, update 
       ! index identifer and minsize

       if (csize < minsize) then

          min_n = n
          minsize = csize

       end if

    end do

    get_bestkn(2) = min_n
    get_bestkn(1) = pert%npert - min_n - 1

  end function


  recursive function kn_prod(np, lvl, dims, k_or_n, sofar) result (kn_p)

    implicit none

    integer :: np, lvl, k_or_n, sofar, kn_p, i
    integer, dimension(np) :: dims

    kn_p = 0

    if (lvl < k_or_n) then

       do i = 1, np

          kn_p = kn_p + kn_prod(np - 1, lvl + 1, &
                                rm_dim(np, dims, i), k_or_n, dims(i)*sofar)

       end do

    else

       kn_p = sofar

    end if

  end function


  function rm_dim(np, dims, skip)

    implicit none

    integer :: skip, i, j, np
    integer, dimension(np) :: dims
    integer, dimension(np - 1) :: rm_dim

    j = 1

    do i = 1, np

       if (i == skip) then

       else

          rm_dim(j) = dims(i)
  
          j = j + 1   

       end if

    end do

  end function


end module
